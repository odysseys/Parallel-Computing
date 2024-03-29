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

#include "matrix_mul.h"
#include<iostream>
#include <stdio.h>

namespace sequential
{

	void matrix_multiplication(	float *sq_matrix_1, float *sq_matrix_2,
					float *sq_matrix_result, unsigned int sq_dimension	)
	{
		for( unsigned int i = 0; i < sq_dimension; i++ )
		{
			for( unsigned int j = 0; j < sq_dimension; j++ )
			{
        //printf("x1 = %d, y1 = %d, value = %d\n", j, i, sq_matrix_1[i * sq_dimension + j] );
        //printf("x2 = %d, y2 = %d, value = %d\n", j, i, sq_matrix_2[i * sq_dimension + j] );
				float sum = 0.0f;
				for( unsigned int k = 0; k < sq_dimension; k++ )
					sum += sq_matrix_1[i*sq_dimension + k] * sq_matrix_2[k*sq_dimension + j];
				sq_matrix_result[i*sq_dimension + j] = sum;
			}
		}
	}

} //namespace sequential


