#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "generate_system_matrix.h"
#include "recon_structs.h"

#include <vector>

//#include <Eigen/Dense>
//#include <Eigen/SparseCore>

#include "spinner.h"

#include "boost/numeric/ublas/vector_sparse.hpp" 
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
namespace ublas = boost::numeric::ublas;

#define PI 3.14159265359

#define debug_disp(VARIABLE_NAME) std::cout << #VARIABLE_NAME": " << VARIABLE_NAME << std::endl; 

void save_system_matrix(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix);

void generate_system_matrix(const struct recon_params * rp, struct ct_data * data){
    std::cout << "Generating system matrix..." << std::endl;

    ublas::vector<double> source_position(3);
    ublas::vector<double> direction_cosine(3);
    ublas::vector<double> e_w(3);
    ublas::vector<double> e_u(3);
    ublas::vector<double> e_z(3);

    std::cout << rp->num_voxels_x << " : " << rp->num_voxels_y << std::endl;
    
    std::vector<ublas::compressed_vector<float>> system_matrix(rp->num_voxels_x*rp->num_voxels_y);
    for (int i = 0; i < system_matrix.size(); i++){
        system_matrix[i].resize((size_t)(rp->Readings*rp->n_channels*rp->Nrows_projection));
    }

    // Report initial information about system matrix
    std::cout << "System matrix size: " << system_matrix.size() << std::endl;
    std::cout << "Allocating system matrix..." << std::endl;

    // Attempts with BOOST
    if ((!rp->Zffs)&(!rp->Phiffs)){
        std::cout << "Matrix generation for SciFFSNone configuration."  << std::endl;

        init_spinner();

        for (int i = 0; i < rp->num_views_for_system_matrix; i++){

            update_spinner(i,rp->num_views_for_system_matrix);

            // tk this could be one place where we take the initial tube angle into account
            double tube_angle = (i - 0.5*rp->num_views_for_system_matrix)*rp->tube_angle_increment;
 
            //Define local coordinate system for the current view
            e_w(0) = cos(tube_angle);
            e_w(1) = sin(tube_angle);
            e_w(2) = 0.0;

            e_u(0) = -sin(tube_angle);
            e_u(1) = cos(tube_angle);
            e_u(2) = 0.0;

            e_z(0) = 0.0;
            e_z(1) = 0.0;
            e_z(2) = 1.0;

            source_position = rp->focal_spot_radius*e_w + rp->table_direction*((double)i - 0.5*(double)rp->num_views_for_system_matrix)*rp->tube_z_increment*e_z;

            for (int j = 0; j < rp->n_channels; j++){
                double transaxial_position = (j - rp->center_channel_non_ffs)*rp->transaxial_detector_spacing;
                double transaxial_angle = transaxial_position / rp->source_detector_distance;
                double cos_transaxial_angle = cos(transaxial_angle);
                double sin_transaxial_angle = sin(transaxial_angle); 

                for (int k = 0; k < rp->Nrows_projection; k++){
                    double axial_position = (k - rp->center_row)*rp->axial_detector_spacing;

                    direction_cosine = (-rp->source_detector_distance*cos_transaxial_angle*e_w  \
                                        -rp->source_detector_distance*sin_transaxial_angle*e_u  \
                                        -axial_position*e_z)/sqrt(axial_position*axial_position \
                                        +rp->source_detector_distance*rp->source_detector_distance);

                    int q = k + rp->Nrows_projection*j + rp->Nrows_projection*rp->n_channels*i;
                    
                    //Compute the contribution of the ray to the system matrix
                    double x, y,
                        x_hat, y_hat, z_hat,
                        i_hat, j_hat, k_hat;

                    double abs_alpha_x = fabs(direction_cosine(0));
                    double abs_alpha_y = fabs(direction_cosine(1));

                    if (abs_alpha_x > abs_alpha_y){

                        //std::cout << "x focused" << std::endl;
                        
                        for (int i = 0; i < rp->num_voxels_x; i++){
                            x = (i - rp->center_voxel_x)*rp->voxel_size_x;

                            y_hat = source_position(1) + direction_cosine(1) / direction_cosine(0)*(x - source_position(0));
                            j_hat = y_hat / rp->voxel_size_y + rp->center_voxel_y;
                            int j = (int)(floor(j_hat));

                            z_hat = source_position(2) + direction_cosine(2) / direction_cosine(0)*(x - source_position(0));
                            k_hat = z_hat / rp->voxel_size_z + rp->center_voxel_z;
                            int k = (int)(floor(k_hat)); // for negative values of k, casting to int gives the wrong answer so I use floor

                            if ((j >= 0) && (j <= (int)(rp->num_voxels_y - 2))){
                                double scale = rp->voxel_size_x / abs_alpha_x;
                                
                                if (k == -1){
                                    system_matrix[i + rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(1 - (j_hat - j))));
                                    system_matrix[i + rp->num_voxels_x*(j + 1)].push_back(q, (float)(scale*(k_hat - k)*(j_hat - j)));                                    
                                }
                                if (k == 0){
                                    system_matrix[i + rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(1 - (j_hat - j))));
                                    system_matrix[i + rp->num_voxels_x*(j + 1)].push_back(q, (float)(scale*(1 - (k_hat - k))*(j_hat - j)));
                                }
                            }
                        }
                    }
                    else{

                        //std::cout << "y focused" << std::endl;
                        
                        for (int j = 0; j < rp->num_voxels_y; j++){                            
                            y = (j - rp->center_voxel_y)*rp->voxel_size_y;

                            x_hat = source_position(0) + direction_cosine(0) / direction_cosine(1)*(y - source_position(1));
                            i_hat = x_hat / rp->voxel_size_x + rp->center_voxel_x;
                            int i = (int)(floor(i_hat));

                            z_hat = source_position(2) + direction_cosine(2) / direction_cosine(1)*(y - source_position(1));
                            k_hat = z_hat / rp->voxel_size_z + rp->center_voxel_z;
                            int k = (int)(floor(k_hat));

                            if ((i >= 0) && (i <= (int)(rp->num_voxels_x - 2))){
                                double scale = rp->voxel_size_y / abs_alpha_y;

                                if (k == -1){
                                    system_matrix[i +       rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(1 - (i_hat - i))));
                                    system_matrix[(i + 1) + rp->num_voxels_x*j].push_back(q, (float)(scale*(k_hat - k)*(i_hat - i)));                                    
                                }
                                if (k == 0){
                                    system_matrix[i +       rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(1 - (i_hat - i))));
                                    system_matrix[(i + 1) + rp->num_voxels_x*j].push_back(q, (float)(scale*(1 - (k_hat - k))*(i_hat - i)));
                                }
                            }
                        }
                    }
                }                
            }
        }

        // Clear our user feeback
        destroy_spinner();

        // Save matrix to disk
        std::ofstream output_filepath(rp->output_dir + "/matrix.bin",std::ios_base::binary);
        for (int i=0; i<system_matrix.size(); i++){
            size_t num_nonzeros=(size_t)system_matrix[i].filled();
            output_filepath.write((char*)&num_nonzeros, sizeof(num_nonzeros));
            if (num_nonzeros > 0){
                struct pair{
                    int index;
                    float value;
                } current_pair;

                std::vector<pair> nonzeros;

                for (auto pos = system_matrix[i].begin(); pos != system_matrix[i].end(); ++pos){
                    current_pair.index = (int)pos.index();
                    current_pair.value = (*pos);

                    nonzeros.push_back(current_pair);
                }
                output_filepath.write((char*)&nonzeros[0], num_nonzeros*sizeof(pair));
            }
            ublas::compressed_vector<float>(system_matrix[i].size(), 0).swap(system_matrix[i]);
        }
        output_filepath.close();

        std::cout << "Done!" << std::endl;
    }
}

void save_system_matrix(struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix){
}
