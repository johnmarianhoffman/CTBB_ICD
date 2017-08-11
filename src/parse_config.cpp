#include <iostream>
#include "recon_structs.h"
#include "parse_config.h"
#include <yaml-cpp/yaml.h>

// Macro that actually reads our YAML and handles keyword detection
#define parse_item(TAG_NAME,TYPE) try{\
    std::cout << #TAG_NAME": " ;                                      \
    for (int i=0; i< 35-sizeof(#TAG_NAME);i++){std::cout << " ";}     \
    rp->TAG_NAME=config[#TAG_NAME].as<TYPE>();                        \
    std::cout << "FOUND "  << std::endl;                              \
    }catch(YAML::RepresentationException& e){std::cout << "NOT FOUND" << std::endl;}

void parse_config(std::string config_file, struct recon_params * rp){
    // Load our YAML config file
    std::cout << "Config File: " << config_file << std::endl;    
    YAML::Node config = YAML::LoadFile(config_file);
    std::cout << std::endl;

    // Pertinent paths block
    parse_item(initial_recon_path,std::string);
    parse_item(matrix_path,std::string);
    parse_item(output_dir,std::string);
    parse_item(output_file,std::string);
    parse_item(potential,std::string);
    parse_item(recon_path,std::string);
    parse_item(sinogram_path,std::string);

    // Scanner Geometry
    parse_item(acquisition_fov,double);
    parse_item(anode_angle,double);
    parse_item(axial_detector_spacing,double);
    parse_item(axial_focal_spot_shift,double);
    parse_item(center_channel_non_ffs,double);
    parse_item(center_row,double);
    parse_item(focal_spot_radius,double);
    parse_item(n_channels,size_t);
    parse_item(num_views_per_turn_without_ffs,size_t);
    parse_item(source_detector_distance,double);
    parse_item(transaxial_detector_spacing,double);
    parse_item(transaxial_focal_spot_shift,double);

    // Iterative Recon parameters
    parse_item(lambda,double);
    parse_item(num_iterations,size_t);
    parse_item(num_views_for_system_matrix,size_t);

    // Scan specific parameters
    parse_item(first_view_angle,double);
    parse_item(num_z_ffs,size_t);
    parse_item(num_phi_ffs,size_t);
    parse_item(num_views,size_t);
    parse_item(table_feed_per_rotation,double);
    parse_item(tube_angle_increment,double);
    parse_item(z_end,double);
    parse_item(z_first_view,double);
    parse_item(z_good_slice,double);
    parse_item(z_matrix_start,double);
    parse_item(z_start,double);

    // Recon Geometry
    parse_item(recon_fov,double);
    parse_item(nx,size_t);
    parse_item(ny,size_t);
    parse_item(slice_thickness,double);

    // John's Additions
    parse_item(FileType,size_t);
    parse_item(FileSubType,size_t);
    parse_item(RawOffset,size_t);
    parse_item(Readings,size_t);
    parse_item(Zffs,size_t);
    parse_item(Phiffs,size_t);
    parse_item(CollSlicewidth,double);
    parse_item(Nrows,size_t);
    parse_item(wfbp_initialize,size_t);

    // Deprecated
    // ========================================
    // parse_item(dx,double);
    // parse_item(dy,double);
    // parse_item(dz,double);
    // parse_item(fov_radius,double);
    // parse_item(center_voxel_x,double);
    // parse_item(center_voxel_y,double);
    // parse_item(center_voxel_z,double);
    // parse_item(nz,size_t);
}
