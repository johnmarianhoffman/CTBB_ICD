
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "generate_system_matrix_ffs.h"
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

struct pair{
    int index;
    float value;
};

ublas::vector<double> generate_ffs_offset(int proj_idx,double da, double dr, double anode_angle, int ZFFS, int PHIFFS);

// Defined in generate_system_matrix.cpp (code available at bottom of this file however)
extern void save_system_matrix_first_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block);
extern void       save_system_matrix_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block);

void generate_system_matrix_ffs(const struct recon_params * rp, struct ct_data * data){

    std::cout << "Generating system matrix (with flying focal spots)..." << std::endl;
    
    // Configure some data for our flying focal spots
    // The focal spot deflection amounts
    double dr=rp->source_detector_distance*rp->CollSlicewidth/(4.0*(rp->source_detector_distance-rp->focal_spot_radius)*tan(rp->anode_angle*PI/180.0));
    double da=rp->source_detector_distance*rp->focal_spot_radius*sin(rp->fan_angle_increment)/(4.0*(rp->source_detector_distance-rp->focal_spot_radius));

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

    if (rp->Phiffs && !rp->Zffs)
        std::cout << "Matrix generation for SciFFSPhi configuration."  << std::endl;
    else if (!rp->Phiffs && rp->Zffs)
        std::cout << "Matrix generation for SciFFSZ configuration."  << std::endl;
    else if (rp->Phiffs && rp->Zffs)
        std::cout << "Matrix generation for SciFFSPhiZ configuration."  << std::endl;
    else{
        std::cout << "You shouldn't be here. Something's gone wrong. Exiting." << std::endl;
        exit(1);
    }
    
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

        // Determine base source position (if we had no FFS) 
        source_position = rp->focal_spot_radius*e_w + rp->table_direction*((double)i - 0.5*(double)rp->num_views_for_system_matrix)*rp->tube_z_increment*e_z;
        // Determine the source position deflection caused by the current FFS
        ublas::vector<double> ffs_offset(3);
        ffs_offset=generate_ffs_offset(i,da,dr,rp->anode_angle,rp->Zffs,rp->Phiffs);

        // Tranform the deflection out of the rotating coordinate frame of the x-ray source, into the gantry coordinate frame
        ublas::vector<double> ffs_offset_rot(3);
        ffs_offset_rot(0) = ublas::inner_prod(ffs_offset,e_w);
        ffs_offset_rot(1) = ublas::inner_prod(ffs_offset,e_u);
        ffs_offset_rot(2) = ffs_offset(2);
        
        // Create a "true" source position with the current deflected focal spot
        ublas::vector<double> ffs_source_position = source_position+ffs_offset_rot; // +inc_source_position;
                
        for (int j = 0; j < rp->n_channels; j++){
            double transaxial_position = (j - rp->center_channel_non_ffs)*rp->transaxial_detector_spacing;
            double transaxial_angle = transaxial_position / rp->source_detector_distance;
            double cos_transaxial_angle = cos(transaxial_angle);
            double sin_transaxial_angle = sin(transaxial_angle); 

            for (int k = 0; k < rp->Nrows_projection; k++){
                double axial_position = (k - rp->center_row)*rp->axial_detector_spacing;

                //b-a = - D cos g ew - D sin g eu - z ez 
                //direction_cosine = (b-a-inc)/||b-a-inc||
                direction_cosine = (-rp->source_detector_distance*cos_transaxial_angle*e_w \
                                    -rp->source_detector_distance*sin_transaxial_angle*e_u  \
                                    -axial_position*e_z-ffs_offset_rot);
                direction_cosine = direction_cosine/ublas::norm_2(direction_cosine);

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

                        y_hat = ffs_source_position(1) + direction_cosine(1) / direction_cosine(0)*(x - ffs_source_position(0));
                        j_hat = y_hat / rp->voxel_size_y + rp->center_voxel_y;
                        int j = (int)(floor(j_hat));

                        z_hat = ffs_source_position(2) + direction_cosine(2) / direction_cosine(0)*(x - ffs_source_position(0));
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

                        x_hat = ffs_source_position(0) + direction_cosine(0) / direction_cosine(1)*(y - ffs_source_position(1));
                        i_hat = x_hat / rp->voxel_size_x + rp->center_voxel_x;
                        int i = (int)(floor(i_hat));

                        z_hat = ffs_source_position(2) + direction_cosine(2) / direction_cosine(1)*(y - ffs_source_position(1));
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
            // Save system matrix chunk to disk (should happen ~5 times throughout generation)
            if (i==300){
                save_system_matrix_first_block(rp,system_matrix);
            }
            else if (i%300==0 && i!=0){
                save_system_matrix_block(rp,system_matrix);
            }
    }

    // Save the final chunk to disk
    save_system_matrix_block(rp,system_matrix);

    // Clear our user feeback
    destroy_spinner();
    std::cout << "Done!" << std::endl;
    
}

ublas::vector<double> generate_ffs_offset(int proj_idx,double da, double dr, double anode_angle, int ZFFS, int PHIFFS){
    ublas::vector<double> ffs_offset(3);
    // Phi-only
    if (PHIFFS && !ZFFS){
        // There are more clever ways to do this, but we implement this way for clarity
        int rho = proj_idx%2;
        if (rho==0){
            ffs_offset(0)=0.0;
            ffs_offset(1)=da;
            ffs_offset(2)=0.0;                
        }
        else{
            ffs_offset(0)=0.0;
            ffs_offset(1)=-da;
            ffs_offset(2)=0.0;                
        }
    }
    // Z-only
    else if (!PHIFFS && ZFFS){
        int rho = proj_idx%2;
        if (rho==0){
            ffs_offset(0) = -dr;
            ffs_offset(1) = 0.0;
            ffs_offset(2) = -dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }
        else{
            ffs_offset(0) = dr;
            ffs_offset(1) = 0.0;
            ffs_offset(2) = dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }
    }
    // Z & Phi ffs
    else if (PHIFFS && ZFFS){
        int rho = proj_idx%4;
        if (rho==0){
            ffs_offset(0) = -dr; 
            ffs_offset(1) = da; 
            ffs_offset(2) = -dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }                    
        else if (rho==1){
            ffs_offset(0) = -dr; 
            ffs_offset(1) = -da; 
            ffs_offset(2) = -dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }
        else if (rho==2){
            ffs_offset(0) = dr; 
            ffs_offset(1) = da; 
            ffs_offset(2) = dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }
        else if (rho==3){
            ffs_offset(0) = dr; 
            ffs_offset(1) = -da; 
            ffs_offset(2) = dr*tan(anode_angle); // tk we may have gotten the direction wrong here
        }
        else {}
    }        
    else{
        std::cout << "Something went wrong (you should never have ended up here). Exiting." << std::endl;
        exit(1);                
    }

    return ffs_offset;
}

//void save_system_matrix_first_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block){
//
//    std::cout << "Saving first block to disk..." << std::endl;
//
//    // Open the output file
//    std::string matrix_filepath=rp->output_dir+"/matrix.bin";    
//    std::ofstream output_file(matrix_filepath.c_str(),std::ios_base::binary);
//
//    // Loop over all columns of our system matrix
//    for (int i=0; i<rp->num_voxels_x*rp->num_voxels_y; i++){
//
//        // Get the information about the column chunk in memory
//        size_t num_nonzeros_memory=(size_t)system_matrix_block[i].filled();
//        struct pair * col_memory=new struct pair[num_nonzeros_memory];
//        
//        size_t nonzero_idx=0;
//        for (auto pos = system_matrix_block[i].begin(); pos != system_matrix_block[i].end(); ++pos){
//            col_memory[nonzero_idx].index=(int)pos.index();
//            col_memory[nonzero_idx].value=(*pos);
//            nonzero_idx=nonzero_idx + 1;
//        }
//
//        // Flush the data to disk (more lines of code than needed for clarity);
//        size_t num_nonzeros=num_nonzeros_memory;
//        output_file.write((char*)&num_nonzeros,sizeof(num_nonzeros));
//        output_file.write((char*)col_memory,num_nonzeros_memory*sizeof(struct pair));
//
//        // Empty contents, and force reallocation of matrix column
//        ublas::compressed_vector<float>(system_matrix_block[i].size(), 0).swap(system_matrix_block[i]);
//
//        delete[] col_memory;
//    }
//
//    output_file.close();
//}
//
//
//void save_system_matrix_block(const struct recon_params * rp,std::vector<ublas::compressed_vector<float>> & system_matrix_block){
//
//    std::cout << "Saving next block to disk..." << std::endl;
//
//    // Filepaths
//    std::string tmp_matrix_filepath=rp->output_dir + "/matrix_tmp.bin";
//    std::string matrix_filepath=rp->output_dir+"/matrix.bin";
//
//    // Open our input and output files
//    std::ofstream   output_file(tmp_matrix_filepath.c_str(),std::ios_base::binary);
//    std::ifstream existing_file(matrix_filepath.c_str(),std::ios_base::binary);
//
//    for (int i=0; i<rp->num_voxels_x*rp->num_voxels_y; i++){
//
//        // Read the already saved chunk of the current column from disk
//        size_t num_nonzeros_existing;
//        existing_file.read((char*)&num_nonzeros_existing,sizeof(num_nonzeros_existing));
//
//        struct pair * col_disk = new struct pair[num_nonzeros_existing];
//        existing_file.read((char*)col_disk,num_nonzeros_existing*sizeof(struct pair));
//
//        // Get the information about the chunk in memory
//        size_t num_nonzeros_memory=(size_t)system_matrix_block[i].filled();
//        struct pair * col_memory=new struct pair[num_nonzeros_memory];
//        
//        size_t nonzero_idx=0;
//        for (auto pos = system_matrix_block[i].begin(); pos != system_matrix_block[i].end(); ++pos){            
//            col_memory[nonzero_idx].index=(int)pos.index();
//            col_memory[nonzero_idx].value=(*pos);
//            nonzero_idx=nonzero_idx + 1;
//        }
//
//        // Flush the data to disk (more lines of code than needed for clarity);
//        size_t num_nonzeros=num_nonzeros_existing+num_nonzeros_memory;
//        output_file.write((char*)&num_nonzeros,sizeof(num_nonzeros));
//        output_file.write((char*)col_disk,num_nonzeros_existing*sizeof(struct pair));
//        output_file.write((char*)col_memory,num_nonzeros_memory*sizeof(struct pair));
//
//        // Empty contents, and force reallocation of matrix column
//        ublas::compressed_vector<float>(system_matrix_block[i].size(), 0).swap(system_matrix_block[i]);
//
//        delete[] col_memory;
//        delete[] col_disk;
//    }
//
//    // Overwrite the existing matrix file with the updated temporary file
//    int success=rename(tmp_matrix_filepath.c_str(),matrix_filepath.c_str()); // 0 return means success
//    if (success!=0){
//        std::cout << "There was an error joining matrix chunks. Exiting."  << std::endl;
//        exit(1);
//    }
//
//    output_file.close();
//    existing_file.close();
//}
