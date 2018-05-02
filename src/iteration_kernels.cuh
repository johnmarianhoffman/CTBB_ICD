#include "recon_structs.h"

// Allocate structs
__constant__ struct iterative_params d_ip;

__global__ void test_kernel(float * raw, float * sinogram, float * recon_volume, pair * sys_matrix_column, size_t nnz, size_t central_idx){

}

__global__ void compute_alpha_beta(float *raw, float * sinogram, float * indices,float * values, size_t nnz, int offset,float * alpha, float * beta){
    int m = threadIdx.x + blockIdx.x*blockDim.x;

    if (m<nnz){
        
        int idx     = indices[m];
        float value = values[m];

        size_t data_size=d_ip.n_elem_sinogram;

        if ((m < nnz) && ( idx>-1 ) && (idx < data_size)){
            alpha[m]=value*value;
            beta[m]=value * (raw[idx]-sinogram[idx]); // This will hurt

            atomicAdd(&alpha[m],4);
            
        }
        else{
            alpha[m]=0.0f;
            beta[m]=0.0f; // This will hurt
        }
        
    }
}

__global__ void reduction(int nnz, float * array){
    extern __shared__ float f[];

    // Each thread loads one elem from global to shared
    int tid=threadIdx.x;
    int idx=threadIdx.x+blockIdx.x*blockDim.x;
    f[tid]=array[idx];
    __syncthreads();

    // Do the reduction
    //for (unsigned int s=1; s < blockDim.x; s *= 2) {
    //    if (tid % (2*s) == 0) {
    //        f[tid] += f[tid + s];
    //    }
    //    __syncthreads();
    //}
    
    for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            f[tid] += f[tid + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0)
        array[blockIdx.x] = f[0];
}

__global__ void forward_project(float * sinogram_estimate, float * recon_volume, pair * sys_matrix_column, int offset){
//    size_t central_idx=data->slice_indices[k];
//    // If using flying focal spots, need to guarantee that offset corresponds to the first projection in a FFS "stack"
//    if (rp->Zffs || rp->Phiffs){
//        int n_focal_spots=pow(2.0,rp->Zffs)*pow(2.0,rp->Phiffs);
//        int mod=central_idx%n_focal_spots;
//        central_idx=central_idx-mod;
//    }
//                    
//    int offset = (central_idx - rp->num_views_for_system_matrix/2)*rp->n_channels*rp->Nrows_projection; // offset is a projection index
//                    
//    for (int m = 0; m<num_nonzeros; m++){
//        int index = nonzeros[m].index + offset; // Raw data index
//        if ((index > -1) && (index < data_size)){
//            size_t voxel_idx=i+j*rp->num_voxels_x+k*rp->num_voxels_x*rp->num_voxels_y;
//            sinogram_estimate[index] = sinogram_estimate[index] + reconstructed_image[voxel_idx]*nonzeros[m].value;
//        }
//    }
//

}
