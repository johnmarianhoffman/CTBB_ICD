#include <iostream>
#include <fstream>
#include <string>

#include "recon_structs.h"
#include "setup.h"
#include "initialize_recon_volume.h"
#include "generate_system_matrix.h"
#include "rotate_slices.h"
#include "icd_iteration.h"

struct flags {
    bool testing;
    bool verbose;
    bool timing;
};

void usage(){
    printf("\n");
    printf("usage: ctbb_icd [options] input_prm_file\n\n");
    printf("    Options:\n");
    printf("          -v: verbose.\n");
    printf("    --timing: Display timing information for each step of the recon process\n");
    printf("\n");
    printf("Copyright Stefano Young, John Hoffman 2017\n\n");
    exit(0);
}

inline bool exists(const std::string& name){
    std::ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char ** argv){

    struct ct_data data={};
    struct recon_params rp={}; // initialize parameter structure to zero
    struct flags flags={};
    
    /*--- Parse command line options ---*/ 
    if (argc<2)
        usage();

    for (int i=1; i<argc-1; i++){
        std::string curr_arg = argv[i];

        if (curr_arg.compare("-t")==0)
            flags.testing=true;
        else if (curr_arg.compare("-v")==0)
            flags.verbose=true;
        else if (curr_arg.compare("--timing")==0)
            flags.timing=true;
        else
            usage();
    }
    
    /*--- Parse our configuration file (parse_config.cpp)---*/
    std::string parameter_file = argv[argc-1];
    rp=configure_recon_params(parameter_file);

    /*--- Get raw data (setup.cpp) ---*/
    load_raw(&rp,&data);

    /*--- All parameters now set to final values, map to CONSTANT ---*/
    const struct recon_params rp_const=rp;
    
    /*--- Initialize reconstruction volume (setup.cpp) ---*/
    // Perform wFBP reconstruction if using as input to ICD
    initialize_recon_volume(&rp_const,&data);

    /*--- Generate system matrix (generate_system_matrix.cpp) ---*/
    // If matrix file does not exist, generate
    if (!exists(rp_const.matrix_path)){
        std::cout << "No existing matrix file found for reconstruction." << std::endl;
        generate_system_matrix(&rp_const,&data);
    }
    else{
        std::cout << "Existing matrix file FOUND." << std::endl;
        std::cout << "Using matrix file: " << rp_const.matrix_path << std::endl;
    }

    std::ofstream debug_outfile("/home/john/Desktop/debug_file.bin",std::ios_base::binary | std::ios::out);
    std::cout << "Number of elements: " << rp.nx*rp.ny*rp.num_voxels_z << std::endl;
    debug_outfile.write((char*)&data.recon_volume[0],rp.nx*rp.ny*rp.num_voxels_z*sizeof(float));
    debug_outfile.close();
    std::cout << "Debug file written to desktop." << std::endl;

    /*--- Perform ICD iterations (icd_iteration.cpp) ---*/
    icd_iteration(&rp_const,&data);

    /*--- De-rotate our ICD slices (rotate_slices.cpp) ---*/
    rotate_slices_rotating2fixed(&rp_const,&data);

    /*--- Write data to disk and clean up (setup.cpp) ---*/
    clean_up(&rp,&data);
    
    return 0;    
}


