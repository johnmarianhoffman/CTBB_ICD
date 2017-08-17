#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <chrono>

#include <omp.h>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

#include "spinner.h"
#include "recon_structs.h"
#include "icd_iteration.h"
#include "penalties.h"

#define OMP_N_THREADS 8

void icd_iteration(const struct recon_params * rp, struct ct_data * data){

    size_t data_size = rp->Readings*rp->n_channels*rp->Nrows_projection;
    
    // Allocate sinogram estimate (all zeros)
    double * sinogram_estimate = new double[rp->Readings*rp->n_channels*rp->Nrows_projection]();
    std::vector<double> reconstructed_image(rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z, 0.0);

    // Copy the float recon volume into the vector array (if uninitialized, will just copy zeros);
    for (int i=0; i<rp->num_voxels_x; i++){
        for (int j=0; j<rp->num_voxels_y; j++){
            for (int k=0; k<rp->num_voxels_z; k++){
                size_t idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
                reconstructed_image[idx]=(double)data->recon_volume[idx];
            }
        }
    }

    std::ifstream file(rp->matrix_path, std::ios_base::binary);
    
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
                struct pair{
                    int index;
                    float value;
                };
                std::vector<pair> nonzeros(num_nonzeros);
                if (num_nonzeros > 0)
                    file.read((char*)&nonzeros[0], num_nonzeros*sizeof(pair));

                // Loop over all slices for current x,y
                for (int k=0; k<rp->num_voxels_z; k++){
        
                    size_t central_idx=data->slice_indices[k];
        
                    int offset = (central_idx - rp->num_views_for_system_matrix/2)*rp->n_channels*rp->Nrows_projection;
        
                    for (int m = 0; m<num_nonzeros; m++){
                        int index = nonzeros[m].index + offset; // Raw data index
                        if ((index > -1) && (index < data_size)){
                            size_t voxel_idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
                            sinogram_estimate[index] = sinogram_estimate[index] + reconstructed_image[voxel_idx]*nonzeros[m].value;
                        }
                    }
                }
            }            
        }
        destroy_spinner();
        file.clear();
        file.seekg(0, std::ios_base::beg);
    }

    /// tk debugging:
    // Write reconstruction to disk
    std::ostringstream recon_path;       
    recon_path << rp->output_dir << "/reconstructions/iteration_init" << ".rcn";        
    std::ofstream recon_file(recon_path.str(), std::ios_base::binary);
    recon_file.write((char*)&reconstructed_image[0], rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(reconstructed_image[0]));
    recon_file.close();
    std::cout << "Wrote initial image to disk." << std::endl;
        
    // Write reconstruction to disk
    std::ostringstream sino_est_path;       
    sino_est_path << rp->output_dir << "/reconstructions/sino_estimation" << ".rcn";
    std::ofstream sino_file(sino_est_path.str(), std::ios_base::binary);
    sino_file.write((char*)sinogram_estimate, rp->Readings*rp->Nrows_projection*rp->n_channels*sizeof(double));
    sino_file.close();
    std::cout << "Wrote initial sinogram to disk." << std::endl;    

    //tk end debugging
    

    ublas::compressed_vector<float> col(rp->Readings*rp->n_channels*rp->Nrows_projection);

    // Initialize iterative parameters
    // Current implementation limited to 2D (hard coded)
    struct iterative_params ip;
    initialize_2d_weights(&ip);
    ip.lambda = rp->lambda;
    //ip.delta  = rp->delta; // What to do about this???

    for (int n = 0; n < rp->num_iterations; n++){
        
        std::cout << "Iteration #" << n << std::endl;
        std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();

        init_spinner();
        for (int j = 0; j < rp->num_voxels_y; j++){
            update_spinner(j,rp->num_voxels_y);
            double y = (j - rp->center_voxel_y)*rp->voxel_size_y;
            for (int i = 0; i < rp->num_voxels_x; i++){

                double x = (i - rp->center_voxel_x)*rp->voxel_size_x;

                size_t nnz;
                file.read((char*)&nnz, sizeof(nnz));
                
                int num_nonzeros = (int)nnz; // cast to int to avoid potential issues

                struct pair{
                    int index;
                    float value;
                };

                std::vector<pair> nonzeros(num_nonzeros);

                if (num_nonzeros > 0)
                    file.read((char*)&nonzeros[0], num_nonzeros*sizeof(pair));

                if ((x*x + y*y) < (rp->fov_radius*rp->fov_radius)){

                    int q0 = i + rp->num_voxels_x*j;

                    for (int k = 0; k < (rp->num_voxels_z); k++){ 

                        /// Grab the Z slice locations (spatial+idx)
                        double curr_slice_location=data->slice_locations[k];
                        size_t central_idx=data->slice_indices[k];
                        
                        int q = q0 + rp->num_voxels_x*rp->num_voxels_y*k;

                        // This is the key spot to select slice location (done via the "central_idx" variable)
                        int offset = (central_idx - rp->num_views_for_system_matrix/2)*rp->n_channels*rp->Nrows_projection;
                        
                        double alpha = 0.0;
                        double beta  = 0.0;

#pragma omp parallel num_threads(OMP_N_THREADS)
                        {
#pragma omp for reduction(+:alpha,beta)
                            for (int m = 0; m<num_nonzeros; m++){
                                int index = nonzeros[m].index + offset;
                                
                                if ((index > -1) && (index < data_size)){
                                    alpha += nonzeros[m].value * nonzeros[m].value;
                                    beta  += nonzeros[m].value * ((double)data->raw[index] - sinogram_estimate[index]);
                                }                                
                            }
                        }

                        ip.alpha = alpha;
                        ip.beta  = beta;
                        
                        // Apply selected penalty functions
                        double pixel_update=0.0;
                        if (true/* Quadratic */){
                            pixel_update=quadratic(q,&ip,reconstructed_image);
                        }
                        else if(false/* Edge Preserving*/){
                            pixel_update=edge_preserving(q,&ip,reconstructed_image);
                        }
                        else{
                            std::cout << "Unrecognized penalty selected. Exiting." << std::endl;
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

            }
        }

        destroy_spinner();

        // Close up our timer
        std::chrono::high_resolution_clock::time_point end=std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
        std::cout << duration << " s" << std::endl;

        // Write reconstruction to disk
        std::ostringstream recon_path;       
        recon_path << rp->output_dir << "/reconstructions/iteration" << n << ".rcn";        
        std::ofstream recon_file(recon_path.str(), std::ios_base::binary);
        recon_file.write((char*)&reconstructed_image[0], rp->num_voxels_x*rp->num_voxels_y*rp->num_voxels_z*sizeof(reconstructed_image[0]));
        recon_file.close();

        // "Rewind" the matrix file for the next iteration
        file.clear();
        file.seekg(0, std::ios_base::beg);        
    }

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


