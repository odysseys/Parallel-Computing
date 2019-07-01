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

#include <omp.h>
#include "matrix_mul.h"
//#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#include <iostream>

namespace omp
{
  void transpose (float *m1, float *m1T, unsigned int dim)
  {
    if(dim<64)
    {// normal loop if number of dimension < 64
      unsigned int i,j;
      #pragma omp parallel for schedule(dynamic)
      for (i = 0; i < dim; i++)
      {
          for (j = 0; j < dim; j++)
          {
            m1T[j * dim + i] = m1[i * dim + j];
          }
      }
    }   
    else
    {// loop tilting if number of dimension >= 64
      unsigned int i1,i2,j, B;
      B=64;
      {
        for (i1 = 0; i1 < dim; i1+=B)
        {
          for (j = 0; j < dim; j++)
          {
            for (i2=0; i2<std::min(B,dim-i1);i2++)
              m1T[j * dim + i1+i2] = m1[(i1+i2) * dim + j];
          }
        }
      }   
    }
  }

  void
  matrix_multiplication(float *sq_matrix_1, float *sq_matrix_2, float *sq_matrix_result, unsigned int sq_dimension )
  {
    float *sq_matrix_2T; //Transpose matrix of sq_matrix_2
    sq_matrix_2T = new float[sq_dimension * sq_dimension];
    transpose(sq_matrix_2, sq_matrix_2T, sq_dimension);
    float *sq_matrix_ex1, *sq_matrix_ex2T; //Insert zeros to expand sq_matrix to make number of dimension be multiple of 4
    unsigned sq_dimension_ex = sq_dimension; //new dimension after insert zeros
    __m128 vresult = _mm_setzero_ps(); //SIMD items for multiplication
    __m128 vsub_matrix1, vsub_matrix2;
    float read_out[4];
    //Insert zeros to expand matrix to make dimension be multiple of 4 for SIMD
    if(sq_dimension % 4 != 0)
    {
        sq_dimension_ex = sq_dimension - (sq_dimension % 4) + 4;//new number of dimension after extension
        sq_matrix_ex1 = (float *) calloc(sq_dimension * sq_dimension_ex, sizeof(float));
        sq_matrix_ex2T = (float *) calloc(sq_dimension * sq_dimension_ex, sizeof(float));
        for(unsigned int d = 0; d < sq_dimension; d++)//copy the original sq_matrix content into extended one
        {
            memcpy(&sq_matrix_ex1[d * sq_dimension_ex], &sq_matrix_1[d * sq_dimension], sq_dimension * sizeof(float));
            memcpy(&sq_matrix_ex2T[d * sq_dimension_ex], &sq_matrix_2T[d * sq_dimension], sq_dimension * sizeof(float));
        }
        sq_matrix_1 = sq_matrix_ex1;//override the input
        sq_matrix_2T = sq_matrix_ex2T;//override the input
    }
//Matrix multiplication by SIMD
#pragma omp parallel for reduction( +: vresult) schedule(dynamic) private(vsub_matrix1, vsub_matrix2, read_out)
    for (unsigned int i = 0; i < sq_dimension; i++) 
    {	
        for(unsigned int j = 0; j < sq_dimension; j++) 
	{
            vresult = _mm_setzero_ps();//after each point, accumulated container has to be set 0s.
	    for (unsigned int k = 0; k < sq_dimension_ex; k+=4)
            {//SIMD multiplication
                vsub_matrix1 = _mm_load_ps(&sq_matrix_1[i * sq_dimension_ex + k]);//read 4 points in sq_matrix_1
                vsub_matrix2 = _mm_load_ps(&sq_matrix_2T[j * sq_dimension_ex + k]);//read 4 points in sq_matrix2
                vresult = _mm_add_ps(vresult, _mm_mul_ps(vsub_matrix1, vsub_matrix2));//each product accumulates
            }
            _mm_store_ps(read_out, vresult);//take out the accumulated products
            sq_matrix_result[i * sq_dimension + j] = read_out[0] + read_out[1] + read_out[2] + read_out[3];//sum up to get matrix multiplication of point (i, j)
	}
    }// End of parallel region
  }
  
} //namespace omp 
