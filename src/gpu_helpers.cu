#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "recon_structs.h"
#include "gpu_helpers.h"

#define MB 1000000

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
        {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
        }
}

size_t filesize(const char* filename){
        std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
        return (size_t)in.tellg(); 
}

void query_resources(const struct recon_params * rp){

    // Sizes of different things for the current geometry
    size_t projection_size = rp->Nrows*rp->n_channels * sizeof(float);    
    size_t slice_size = rp->nx * rp->ny * sizeof(float);
    size_t matrix_size = filesize(rp->matrix_path.c_str());
    
    std::cout << "Querying GPU memory..." << std::endl; 
    size_t free;
    size_t total;    
    gpuErrchk( cudaMemGetInfo(&free,&total) );

    std::cout << "========================================" << std::endl;
    std::cout << "             GPU memory report:         " << std::endl;
    std::cout << "========================================" << std::endl;    
    
    std::cout << "Memory available:         " << free/MB <<  " of " << total/MB <<  " Mbytes"  << std::endl;
    std::cout << "No. reconstructed slices: " << (int)((float)free/(float)slice_size) << std::endl;
    std::cout << "No. projections:          " << (int)((float)free/(float)projection_size) << std::endl;
    std::cout << "System matrices:          " << (int)((float)free/(float)(matrix_size)) << " ("<< (matrix_size/MB) << " MB)" << std::endl;

}
