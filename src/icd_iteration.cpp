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

    const int num_neighbors = 8;
    std::vector<double> weights(num_neighbors);

    weights[0] = 1 / sqrt(2.0);
    weights[1] = 1.0;
    weights[2] = 1 / sqrt(2.0);
    weights[3] = 1.0;
    weights[4] = 1.0;
    weights[5] = 1 / sqrt(2.0);
    weights[6] = 1.0;
    weights[7] = 1 / sqrt(2.0);

    double weights_scale = 0.0;
    for (int i = 0; i < num_neighbors; i++){
        weights_scale += weights[i];
    }
    
    std::vector<int> neighbor_indices(num_neighbors);

    std::ifstream file(rp->matrix_path, std::ios_base::binary);

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

                    //int k_off=0;

                    // tk There doesn't seem to be any implicit
                    // tk information about slice location contained in
                    // tk this step.  It also seems like while the slice
                    // tk thickness agrees well with WFBP, slice
                    // tk locations are slightly different ultimately
                    // tk causing some disparity between the two
                    // tk reconstructions that we will need to resolve.
                    // tk It's unclear whether this will need to be done
                    // tk in the system matrix step or here in the
                    // tk iteration step.

                    //                   for (int k = 1; k < (rp->num_voxels_z - 1); k++){ // I am skipping edges, but those are still contained in the system matrix
                    for (int k = 0; k < (rp->num_voxels_z); k++){ // I am skipping edges, but those are still contained in the system matrix                        

                        /// Grab the Z slice locations (spatial+idx)
                        double curr_slice_location=data->slice_locations[k];
                        size_t central_idx=data->slice_indices[k];
                        
                        int q = q0 + rp->num_voxels_x*rp->num_voxels_y*k;

                        //////// This is (one of) the key spot(s) to change location of reconstructed slice. Still missing some small detail though. tk
                        //int offset = rp->views_per_slice*(k - k_off)*rp->n_channels*rp->Nrows_projection;                        
                        //int offset = (central_idx - (rp->views_per_slice/2) - rp->views_per_slice*k_off)*rp->n_channels*rp->Nrows_projection;
                        //int offset = central_idx*rp->n_channels*rp->Nrows_projection - (rp->views_per_slice/2)     - rp->views_per_slice*k_off*rp->n_channels*rp->Nrows_projection;
                        int offset = (central_idx - rp->num_views_for_system_matrix/2)*rp->n_channels*rp->Nrows_projection;
                        
                        double alpha = 0.0;
                        double beta = 0.0;

#pragma omp parallel num_threads(OMP_N_THREADS)
                        {
#pragma omp for reduction(+:alpha,beta)
                            for (int m = 0; m<num_nonzeros; m++){
                                int index = nonzeros[m].index + offset;
                                
                                if ((index > -1) && (index < data_size)){
                                    alpha += nonzeros[m].value*nonzeros[m].value;
                                    beta += nonzeros[m].value * ((double)data->raw[index] - sinogram_estimate[index]);
                                }                                
                            }
                        }

                        //First find the indices and weights of the neighbors for pixel q
                        //q=i+NX*j+NX*NY*k

                        //in-plane
                        neighbor_indices[0] = q - rp->num_voxels_x - 1;
                        neighbor_indices[1] = q - rp->num_voxels_x;
                        neighbor_indices[2] = q - rp->num_voxels_x + 1;
                        neighbor_indices[3] = q - 1;
                        neighbor_indices[4] = q + 1;
                        neighbor_indices[5] = q + rp->num_voxels_x - 1;
                        neighbor_indices[6] = q + rp->num_voxels_x;
                        neighbor_indices[7] = q + rp->num_voxels_x + 1;

                        double sum1 = 0.0;
                        for (int n = 0; n < num_neighbors; n++){
                            sum1 += weights[n] * (reconstructed_image[neighbor_indices[n]] - reconstructed_image[q]);
                        }

                        double pixel_update = (beta + rp->lambda*sum1) / (alpha + rp->lambda*weights_scale);

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


