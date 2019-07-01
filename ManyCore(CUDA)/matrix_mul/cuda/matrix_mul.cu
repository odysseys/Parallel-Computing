/*
    Copyright (C) 2011  Abhinav Jauhri (abhinav.jauhri@gmail.com), Carnegie Mellon University - Silicon Valley

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <cuda.h>
#include <cuda_runtime.h>
#include "matrix_mul.h"
#include <stdio.h>
#define TILE_WIDTH 2
#define BLOCK_SIZE 32


namespace cuda
{
  __global__
  void
  matrix_mul_kernel(float *sq_matrix_1, float *sq_matrix_2, float *sq_matrix_result, int sq_dimension)
  {
    //For repeatly access blocks, use shared memory to speed up
    __shared__ float local_mat_1[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ float local_mat_2[BLOCK_SIZE][BLOCK_SIZE];

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int block_offsetx = blockIdx.x * BLOCK_SIZE;
    int block_offsety = blockIdx.y * BLOCK_SIZE;
    float sum = 0.0f;
    
    #pragma unroll
    for(int i = 0; i < sq_dimension; i += BLOCK_SIZE)
    {
      //Transfer to 2-D matrix to avoid memory bank conflict
      //local_mat_1 is row-based moving for each 32*32 block in 1st input matrix
      if(ty + block_offsety < sq_dimension && tx + i < sq_dimension)
        local_mat_1[ty][tx] = sq_matrix_1[(ty + block_offsety) * sq_dimension + tx + i            ];
      else
        local_mat_1[ty][tx] = 0;
      
      //local_mat_2 is column-based moving for each 32*32 block in 2nd input matrix
      if(tx + block_offsetx < sq_dimension && ty + i < sq_dimension)
        local_mat_2[ty][tx] = sq_matrix_2[(ty + i            ) * sq_dimension + tx + block_offsetx];
      else
        local_mat_2[ty][tx] = 0;
      
      //must wait all threads finishing moving data into shared memory
      __syncthreads();
      
      #pragma unroll
      for(int k = 0; k < BLOCK_SIZE; k++)
      {
        sum += local_mat_1[ty][k] * local_mat_2[k][tx];
      }
      
      __syncthreads();//must wait all threads sum up
    }

    if(tx + block_offsetx < sq_dimension && ty + block_offsety < sq_dimension)
    {
      sq_matrix_result[(ty + block_offsety) * sq_dimension + tx + block_offsetx] = sum;
      //calculate the correct position of the product
    }
  }

  void
  matrix_multiplication(float *sq_matrix_1, float *sq_matrix_2, float *sq_matrix_result, unsigned int sq_dimension)
  {
    int size = sq_dimension * sq_dimension * sizeof(float);
    float *sq_matrix_1_d, *sq_matrix_2_d, *sq_matrix_result_d;

    /***************************************************
  1st Part: Allocation of memory on device memory
    ****************************************************/

    /* copy sq_matrix_1 and sq_matrix_2 to device memory */
    cudaMalloc((void**) &sq_matrix_1_d, size);
    cudaMemcpy(sq_matrix_1_d, sq_matrix_1, size, cudaMemcpyHostToDevice);
    cudaMalloc((void**) &sq_matrix_2_d, size);
    cudaMemcpy(sq_matrix_2_d, sq_matrix_2, size, cudaMemcpyHostToDevice);

    /*allocate sq_matrix_result on host */
    cudaMalloc((void**) &sq_matrix_result_d, size);

    /***************************************************
   2nd Part: Inovke kernel
    ****************************************************/
    //fix the size to maximum number 1024
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    //if sq_dimension*sq_dimension>1024, use other blocks to calculate
    int gridx = int (sq_dimension + BLOCK_SIZE - 1) / int(BLOCK_SIZE);
    dim3 dimGrid(gridx, gridx);
    matrix_mul_kernel<<<dimGrid, dimBlock>>>(sq_matrix_1_d, sq_matrix_2_d, sq_matrix_result_d, sq_dimension);
    //, dimBlock.x * dimBlock.x * sizeof(float)
    /***************************************************
   3rd Part: Transfer result from device to host
    ****************************************************/
    cudaMemcpy(sq_matrix_result, sq_matrix_result_d, size, cudaMemcpyDeviceToHost);
    cudaFree(sq_matrix_1_d);
    cudaFree(sq_matrix_2_d);
    cudaFree(sq_matrix_result_d);
  }
} // namespace cuda
