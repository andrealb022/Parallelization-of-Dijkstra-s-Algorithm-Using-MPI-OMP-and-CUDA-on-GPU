/* 
 * Course: High Performance Computing 2023/2024 
 * 
 * Lecturer: Francesco Moscato	fmoscato@unisa.it 
 *
 * Student : 
 * Alberti Andrea	0622702370	a.alberti2@studenti.unisa.it
 *
 * 
 * Copyright (C) 2023 - All Rights Reserved 
 *
 * This file is part of DijkstraFinalProjectHPC 
 *
 * DijkstraFinalProjectHPC   is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 *
 * DijkstraFinalProjectHPC   is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 *
 * You should have received a copy of the GNU General Public License 
 * along with DijkstraFinalProjectHPC . If not, see <http://www.gnu.org/licenses/>. 
/**
 *                                  REQUIREMENTS OF THE ASSIGNMENT
 *
 * Student shall provide a parallel version of Dijkstra algorithm with both "OpenMP + MPI" and "OpenMP + Cuda" approaches,
 * comparing results with a known solution on single-processing node. 
 * Results and differences shall be discussed for different inputs (type and size).
 * The parallel algorithm used in "OpenMP + MPI" solution could not be the same of the "OpenMP + CUDA" approach.
 *
 * @file utility.h
 * @copyright Copyright (c) 2023
 */

#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <limits.h>

#define MAXINT INT_MAX    //INFINITY
#define SOURCE 0          //SOURCE PROCESS IN MPI
#define SOURCE_VERTEX 0   //SOURCE VERTEX
#define FALSE 0
#define TRUE 1

/**
 * @brief Prints a matrix to the console.
 * 
 * @param N Number of rows/columns in the matrix.
 * @param adj_matrix The matrix to be printed.
 */
void print_matrix(int N, int* adj_matrix);

/**
 * @brief Prints a vector to the console.
 * 
 * @param N Length of the vector.
 * @param vec The vector to be printed.
 */
void print_vector(int N, int* vec);

/**
 * @brief Prints a matrix to a file.
 * 
 * @param fp File pointer.
 * @param N Number of rows/columns in the matrix.
 * @param adj_matrix The matrix to be printed.
 */
void print_matrix_on_file(FILE* fp, int N, int* adj_matrix);

/**
 * @brief Prints a vector to a file.
 * 
 * @param fp File pointer.
 * @param N Length of the vector.
 * @param vec The vector to be printed.
 */
void print_vector_on_file(FILE* fp, int N, int* vec);

/**
 * @brief Initializes an array with a specified value.
 *
 * @param array The array to be initialized.
 * @param size  The size of the array.
 * @param value The value to set for each element.
 */
void initializeArray(int* array, int size, int value);


#endif // UTILITY_H

