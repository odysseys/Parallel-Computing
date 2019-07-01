#include <omp.h>
#include <stdlib.h>
#include "matrix_mul.h"
#include <memory.h>
#define min(a,b)  ((a < b) ? a : b)

namespace omp
{
  /*Transpose matrix so that the product would be symmetric*/
  void transpose(float *sq_matrix_1, float *sq_matrix_2, unsigned int sq_dimension)
  {
    /*Parallel process the matrix transpose*/
    #pragma omp parallel for
    for(unsigned int i = 0; i < sq_dimension; i++)
    {
         for(unsigned int j = 0; j < sq_dimension; j++)
         { 
            sq_matrix_2[j*sq_dimension+i] = sq_matrix_1[i*sq_dimension+j];
         }
     }
  }

  /*Naive matrix multiply*/
  void matrix_multiplication(float *sq_matrix_1, float *sq_matrix_2, float *sq_matrix_result, unsigned int sq_dimension)
  {
    unsigned int block_size = 64;
    /*As matrix result would be variable size, use memset which copy 0 into matrix to initialize*/
   unsigned int result_size = sq_dimension * sq_dimension * sizeof(float);
    memset(sq_matrix_result, 0, result_size);
   
    /*Allocate the memory for transpose matrix*/ 
    float *M = (float*)malloc(sizeof(float)*sq_dimension*sq_dimension);
    transpose(sq_matrix_2, M, sq_dimension);

    #pragma omp parallel for 
    for (unsigned int i = 0; i < sq_dimension; i+=block_size)
    {
       for(unsigned int j = 0; j < sq_dimension; j+=block_size)
       {
          for(unsigned int k = 0; k < sq_dimension; k++)
          {
             /*Partition large array into smaller block*/
             for(unsigned int a = i; a < min(sq_dimension,(i+block_size)); a++)
              {  
                /*Initial block_sum to record the result of block*/
                float block= 0;
		 for(unsigned int b = j; b < min(sq_dimension,(j+block_size)); b++)
                 {
		  block+= sq_matrix_1[k*sq_dimension+b] * M[a*sq_dimension+b];
                 }
                 /*Assign the result back after sum of block result to reduce the cumputation */
                  sq_matrix_result[k*sq_dimension+a] += block;
              }
             }
         }
      }
  }
} //namespace omp
