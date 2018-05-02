/* FreeCT_ICD is MBIR CT reconstruction Software */
/* Copyright (C) 2018  John Hoffman, Stefano Young, Frederic Noo */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

/* Questions and comments should be directed to */
/* jmhoffman@mednet.ucla.edu with "FreeCT_ICD" in the subject line*/

#include <stdio.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include <chrono>

#include <omp.h>

#include "spinner.h"
#include "recon_structs.h"
#include "icd_iteration.h"
#include "penalties.h"

#include "iteration_kernels.cuh"
#include "icd_iteration_gpu.h"
#include "gpu_helpers.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
        {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
        }
}


#define OMP_N_THREADS 11

void icd_iteration_gpu(const struct recon_params * rp, struct ct_data * data){
    // ==================================================
    // ALLOCATE SOME KEY ARRAYS
    // ==================================================    
    size_t data_size = rp->Readings*rp->n_channels*rp->Nrows_projection;
    
    // Allocate sinogram estimate (all zeros)
    float * sinogram_estimate = new float[rp->Readings*rp->n_channels*rp->Nrows_projection]();
    float * reconstructed_image= new float[rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z];

    // Copy the float recon volume into the vector array (if uninitialized, will just copy zeros);
    for (int i=0; i<rp->num_voxels_x; i++){
        for (int j=0; j<rp->num_voxels_y; j++){
            for (int k=0; k<rp->num_voxels_z; k++){
                size_t idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
                reconstructed_image[idx]=data->recon_volume[idx];
            }
        }
    }

    std::ifstream file(rp->matrix_path, std::ios_base::binary);

    // Initialize iterative parameters
    // Current implementation limited to 2D (hard coded)
    struct iterative_params ip;
    initialize_2d_weights(&ip);
    ip.lambda = rp->lambda;
    ip.Nx=rp->num_voxels_x;
    ip.Ny=rp->num_voxels_y;
    ip.Nz=rp->num_voxels_z;    
    ip.delta  = rp->delta;
    ip.n_elem_sinogram = rp->n_channels*rp->Nrows*rp->Readings;
    
    // ==================================================
    // GPU STUFF: Start 
    // ==================================================
    // Do this manually for the time being
    cudaSetDevice(1);
    cudaDeviceReset();
    query_resources(rp);

    // Allocate our GPU array
    float * d_sinogram;
    float * d_sinogram_estimate;
    float * d_recon_array;
    float * d_values;
    float * d_indices;

    cudaMalloc(&d_sinogram,rp->n_channels*rp->Nrows*rp->Readings* sizeof(float));
    cudaMalloc(&d_sinogram_estimate,rp->n_channels*rp->Nrows*rp->Readings* sizeof(float));    
    cudaMalloc(&d_recon_array,rp->nx*rp->ny*rp->num_voxels_z*sizeof(float));

    // Copy data to GPU
    cudaMemcpyToSymbol(d_ip,&ip,sizeof(struct iterative_params),0,cudaMemcpyHostToDevice);
    cudaMemcpy(d_sinogram,data->raw,rp->n_channels*rp->Nrows*rp->Readings*sizeof(float),cudaMemcpyHostToDevice);
    //cudaMemcpy(d_sinogram_estimate,sinogram_estimate,rp->n_channels*rp->Nrows*rp->Readings*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_recon_array,reconstructed_image,ip.Nx*ip.Ny*ip.Nz*sizeof(float),cudaMemcpyHostToDevice);

    query_resources(rp);

    // ==================================================
    // GPU STUFF: End 
    // ==================================================
    
    // If WFBP was used to inialize the reconstructions, we need to initialize our sinogram estimate.
    if (rp->wfbp_initialize){
        std::cout << "Initializing sinogram estimate..." << std::endl;       

        // run a forward projection to initialize
        init_spinner();
        for (int j=0; j<rp->num_voxels_y; j++){
            update_spinner(j,rp->num_voxels_x);
            for (int i=0; i<rp->num_voxels_x; i++){        
                // Extract column of projection matrix
                size_t nnz;
                file.read((char*)&nnz, sizeof(nnz));

                int num_nonzeros = (int)nnz; // cast to int to avoid potential issues
                
                struct pair * nonzeros = new struct pair[num_nonzeros];
                if (num_nonzeros > 0)
                    file.read((char*)&nonzeros[0], num_nonzeros*sizeof(pair));

                // We need to break the pair array into a vector of indices and values
                int32_t * indices = new int32_t[nnz];
                float * values    = new float[nnz];

                for (int ii=0;i<nnz;i++){
                    indices[ii]=nonzeros[ii].index;
                    values[ii]=nonzeros[ii].value;
                }

                cudaMalloc(&d_indices,nnz*sizeof(int32_t));
                cudaMalloc(&d_values,nnz*sizeof(float));
                gpuErrchk(cudaMemcpy(d_values,values,nnz*sizeof(float),cudaMemcpyHostToDevice));
                gpuErrchk(cudaMemcpy(d_indices,indices,nnz*sizeof(int32_t),cudaMemcpyHostToDevice));
                
                // Loop over all slices for current x,y
                for (int k=0; k<rp->num_voxels_z; k++){
                    ///////////////////////////////////////////////////////////////////////////
                    // Largest *possible* number of nonzeros would be the size of the sinogram
                    // We'll never use this much, however we can read memory up to that size
                    // without segfaulting.
                    ///////////////////////////////////////////////////////////////////////////
                    dim3 allowed_block_size=get_max_threads(1);// max threads for device 1

                    size_t total_threads=32*ceil((float)nnz/32.0f);
                    //total_threads=pow(2, ceil(log(total_threads)/log(2))); // Round up to power of 2 for efficient reduction
                    dim3 threads(min((int)total_threads,(int)allowed_block_size.x),1);
                    dim3 blocks(ceil((float)total_threads/(float)allowed_block_size.x),1);
        
                    size_t central_idx=data->slice_indices[k];
                    // If using flying focal spots, need to guarantee that offset corresponds to the first projection in a FFS "stack"
                    if (rp->Zffs || rp->Phiffs){
                        int n_focal_spots=pow(2.0,rp->Zffs)*pow(2.0,rp->Phiffs);
                        int mod=central_idx%n_focal_spots;
                        central_idx=central_idx-mod;
                    }
                    
                    int offset = (central_idx - rp->num_views_for_system_matrix/2)*rp->n_channels*rp->Nrows_projection; // offset is a projection index
                    
#pragma omp parallel num_threads(OMP_N_THREADS) 
                    {
#pragma omp for
                    for (int m = 0; m<num_nonzeros; m++){
                        int index = nonzeros[m].index + offset; // Raw data index
                        if ((index > -1) && (index < data_size)){
                            size_t voxel_idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
                            sinogram_estimate[index] = sinogram_estimate[index] + reconstructed_image[voxel_idx]*nonzeros[m].value;
                        }
                    }
                    }
                }

                // Free all the things
                cudaFree(d_indices);
                cudaFree(d_values);                
                delete[] nonzeros;
                delete[] indices;
                delete[] values;                
            }            
        }
        destroy_spinner();
        file.clear();
        file.seekg(0, std::ios_base::beg);
    }


    // Write pre-iteration reconstruction to disk 
    std::ostringstream recon_path;       
    recon_path << rp->output_dir << "/reconstructions/iteration0.rcn";
    std::ofstream recon_file(recon_path.str(), std::ios_base::binary);
    recon_file.write((char*)&reconstructed_image[0], rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(reconstructed_image[0]));
    recon_file.close();
    std::cout << "Wrote initial image to disk." << std::endl;
        
    // Write reconstruction to disk
    std::ostringstream sino_est_path;       
    sino_est_path << rp->output_dir << "/reconstructions/sino_estimation" << ".rcn";
    std::ofstream sino_file(sino_est_path.str(), std::ios_base::binary);
    sino_file.write((char*)sinogram_estimate, rp->Readings*rp->Nrows_projection*rp->n_channels*sizeof(float));
    sino_file.close();
    std::cout << "Wrote initial sinogram to disk." << std::endl;

    for (int n = 0; n < rp->num_iterations; n++){
        
        std::cout << "Iteration #" << n+1 << std::endl;
        std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();

        double fov_limit=(rp->acquisition_fov/2.0)*(rp->acquisition_fov/2.0);

        init_spinner();
        for (int j = 0; j < rp->num_voxels_y; j++){
            update_spinner(j,rp->num_voxels_y);
            double y = (j - rp->center_voxel_y)*rp->voxel_size_y;
            for (int i = 0; i < rp->num_voxels_x; i++){

                double x = (i - rp->center_voxel_x)*rp->voxel_size_x;

                // Load column of the system matrix from disk
                size_t nnz;
                file.read((char*)&nnz, sizeof(nnz));                
                int num_nonzeros = (int)nnz; // cast to int to avoid potential issues
                struct pair * nonzeros = new struct pair[num_nonzeros];
                if (num_nonzeros > 0)
                    file.read((char*)&nonzeros[0], num_nonzeros*sizeof(pair));
                else
                    continue;

                // ==================================================
                // GPU STUFF: Start 
                // ==================================================
                // Copy column of system matrix over to GPU

                
                // ==================================================
                // GPU STUFF: end
                // ==================================================

                if ((x*x + y*y) < fov_limit){

                    int q0 = i + rp->num_voxels_x*j;

                    for (int k = 0; k < (rp->num_voxels_z); k++){ 

                        int q = q0 + rp->num_voxels_x*rp->num_voxels_y*k;
                        
                        /// Grab the Z slice locations (spatial+idx)
                        //double curr_slice_location=data->slice_locations[k];
                        size_t central_idx=data->slice_indices[k];

                        // This is the key spot to select slice location (done via the "central_idx" variable)                        
                        // If using flying focal spots, need to guarantee that offset corresponds to the first projection in a FFS "stack"
                        int start_idx=(central_idx - rp->num_views_for_system_matrix/2);
                        int offset = start_idx*rp->n_channels*rp->Nrows_projection;
                        
                        float alpha = 0.0;
                        float beta  = 0.0;

                        // ==================================================
                        // GPU STUFF: Start 
                        // ==================================================

                        // Total number of threads will be equal to
                        // the total number of nonzero elements.
                        // We'll round up to the nearest multiple of
                        // 32.

                        // Largest *possible* number of nonzeros would be the size of the sinogram
                        // We'll never use this much, however we can read memory up to that size
                        // without segfaulting.
                        cudaEvent_t start, stop;
                        cudaEventCreate(&start);
                        cudaEventCreate(&stop);

                        dim3 allowed_block_size=get_max_threads(1);// max threads for device 1

                        size_t total_threads=32*ceil((float)nnz/32.0f);
                        total_threads=pow(2, ceil(log(total_threads)/log(2))); // Round up to power of 2 for efficient reduction
                        dim3 threads(min((int)total_threads,(int)allowed_block_size.x),1);
                        dim3 blocks(ceil((float)total_threads/(float)allowed_block_size.x),1);

                        //printf("Allowed: %d x %d x %d\n",allowed_block_size.x,allowed_block_size.y,allowed_block_size.z);
                        printf("Threads: %d x %d x %d\n",threads.x,threads.y,threads.z);
                        printf(" Blocks: %d x %d x %d\n",blocks.x,blocks.y,blocks.z);

                        // Precompute alpha and beta for all non-zero elements
                        float * d_alpha;
                        float * d_beta;                        
                        cudaMalloc(&d_alpha ,   total_threads*sizeof(float));
                        cudaMalloc(&d_beta  ,   total_threads*sizeof(float));
                        cudaMemset(d_alpha  ,0 ,total_threads*sizeof(float));                        
                        cudaMemset(d_beta   ,0 ,total_threads*sizeof(float));

                        cudaEventRecord(start);                        
                        // Compute alpha and beta, then reduce (i.e. kernel calls)
                        compute_alpha_beta<<<blocks,threads>>>(d_sinogram, d_sinogram_estimate,d_indices,d_values,nnz,offset,d_alpha,d_beta);
                        gpuErrchk(cudaPeekAtLastError());
                        
                        reduction<<<blocks,threads,total_threads*sizeof(float)>>>(nnz,d_alpha);
                        gpuErrchk(cudaPeekAtLastError()) ;
                        
                        reduction<<<blocks,threads,total_threads*sizeof(float)>>>(nnz,d_beta);
                        gpuErrchk(cudaPeekAtLastError());                        

                        // Copy final value back to host
                        cudaMemcpy(&alpha,&d_alpha[0],sizeof(float),cudaMemcpyDeviceToHost);
                        cudaMemcpy(&beta,&d_beta[0],sizeof(float),cudaMemcpyDeviceToHost);                        

                        cudaFree(d_alpha);
                        cudaFree(d_beta);
                        
                        cudaEventRecord(stop);
                        cudaEventSynchronize(stop);
                        float milliseconds = 0;
                        cudaEventElapsedTime(&milliseconds, start, stop);
                        printf("GPU: alpha: %.5f    beta: %.5f\n",alpha,beta);
                        printf("GPU execution time: %.5f\n",milliseconds);
                        
                        alpha=0.0f;
                        beta=0.0f;                        
                        
                        // ==================================================
                        // GPU STUFF: End
                        // ==================================================
                        cudaEventRecord(start);
#pragma omp parallel num_threads(OMP_N_THREADS)
                        {
#pragma omp for reduction(+:alpha,beta)
                            for (int m = 0; m<num_nonzeros; m++){
                                int index = nonzeros[m].index + offset;
                                
                                if ((index > -1) && (index < data_size)){
                                    alpha += nonzeros[m].value * nonzeros[m].value;
                                    beta  += nonzeros[m].value * (data->raw[index] - sinogram_estimate[index]);
                                }                                
                            }
                        }

                        cudaEventRecord(stop);
                        cudaEventSynchronize(stop);
                        cudaEventElapsedTime(&milliseconds, start, stop);
                        printf("CPU: alpha: %.5f    beta: %.5f\n",alpha,beta);
                        printf("CPU execution time: %.5f\n",milliseconds);
                        
                        ip.alpha = alpha;
                        ip.beta  = beta;

                        // Apply selected penalty functions
                        double pixel_update=0.0;
                        /* Quadratic */
                        if (rp->penalty.compare("quadratic")==0){
                            pixel_update=quadratic(q,&ip,reconstructed_image);
                        }
                        /* Edge Preserving*/
                        else if(rp->penalty.compare("edge-preserving")==0){
                            pixel_update=edge_preserving(q,&ip,reconstructed_image);
                        }
                        else{
                            std::cout << "Unrecognized penalty selected. Exiting." << std::endl;
                            exit(1);
                        }

                        //Enforce positivity
                        if (pixel_update+reconstructed_image[q]<0)
                            pixel_update = -reconstructed_image[q];
                        
                        //Update image
                        reconstructed_image[q] += pixel_update;

                        //Update the forward-projection data
#pragma omp parallel num_threads(OMP_N_THREADS)
                        {
#pragma omp for
                            for (int m = 0; m<num_nonzeros; m++){
                                int index = nonzeros[m].index + offset;

                                if ((index > -1) && (index < data_size))
                                    sinogram_estimate[index] += pixel_update*nonzeros[m].value;
                            }
                        }

                    }
                }
                delete[] nonzeros;
            }
        }

        destroy_spinner();

        // Close up our timer
        std::chrono::high_resolution_clock::time_point end=std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
        std::cout << duration << " s" << std::endl;

        // Write reconstruction to disk
        std::ostringstream recon_path;       
        recon_path << rp->output_dir << "/reconstructions/iteration" << n+1 << ".rcn";
        std::ofstream recon_file(recon_path.str(), std::ios_base::binary);
        recon_file.write((char*)&reconstructed_image[0], rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(reconstructed_image[0]));
        recon_file.close();

        // "Rewind" the matrix file for the next iteration
        file.clear();
        file.seekg(0, std::ios_base::beg);        
    }

                    
    // ==================================================
    // GPU STUFF: Start 
    // ==================================================

    // Copy data back from GPU
    //cudaMemcpyToSymbol(d_ip,&ip,sizeof(struct iterative_params),0,cudaMemcpyHostToDevice);
    //cudaMemcpy(d_sinogram,data->raw,rp->n_channels*rp->Nrows*rp->Readings*sizeof(float),cudaMemcpyHostToDevice);
    //cudaMemcpy(d_sinogram_estimate,sinogram_estimate,rp->n_channels*rp->Nrows*rp->Readings*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(reconstructed_image,d_recon_array,ip.Nx*ip.Ny*ip.Nz*sizeof(float),cudaMemcpyDeviceToHost);
    
    cudaFree(d_sinogram);
    cudaFree(d_sinogram_estimate);
    cudaFree(d_recon_array);    

    query_resources(rp);

    // ==================================================
    // GPU STUFF: End
    // ==================================================


    // Copy the final reconstructed volume back into our data structure
    for (int i=0; i<rp->num_voxels_x; i++){
        for (int j=0; j<rp->num_voxels_y; j++){
            for (int k=0; k<rp->num_voxels_z; k++){
                size_t idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
                data->recon_volume[idx]=(float)reconstructed_image[idx];
            }            
        }        
    }

}


