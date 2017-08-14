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
    
    //std::vector<float> data(rp->Readings*rp->n_channels*rp->Nrows_projection, 0.0f);
    //size_t data_size = data.size();
    //std::cout << "Theoretical Size: "  << rp->Readings*rp->n_channels*rp->Nrows << std::endl;
    //std::cout << "Actual size: " << data_size << std::endl;
    //std::cout << "Exiting." << std::endl;
    
    //load_all_frames(data); // ?????? Loads data (we have already done this...)

    // Allocate sinogram estimate (all zeros)
    std::vector<double> sinogram_estimate(rp->Readings*rp->n_channels*rp->Nrows_projection, 0.0);
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

    ublas::compressed_vector<float> col(rp->Readings*rp->n_channels*rp->Nrows_projection);

    std::ifstream file(rp->matrix_path, std::ios_base::binary);

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

                    //double alpha = inner_prod(col, col);
                    // tk This is just a constant. Should it be?
                    int k_off = (int)ceil(0.5*rp->beam_width_at_roi_edge / rp->voxel_size_z) + 1;
                    
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
                        if (pixel_update < -reconstructed_image[q])
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

        //// debug
        //std::ofstream debug_outfile("/home/john/Desktop/debug_file.bin",std::ios_base::binary | std::ios::out);
        //std::cout << "Number of elements: " << rp->n_channels*rp->Nrows_projection*rp->Readings << std::endl;
        //std::cout << "N channels: " << rp->n_channels << std::endl;
        //std::cout << "N_rows: " << rp->Nrows_projection << std::endl;
        //std::cout << "N_projections: " << rp->Readings  <<  std::endl;
        //debug_outfile.write((char*)&sinogram_estimate[0],rp->n_channels*rp->Nrows_projection*rp->Readings*sizeof(double));
        //debug_outfile.close();
        //std::cout << "Debug file written to desktop." << std::endl;

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


