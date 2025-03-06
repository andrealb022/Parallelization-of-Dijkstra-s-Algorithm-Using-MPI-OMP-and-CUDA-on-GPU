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
 * @file graph.h
 * @copyright Copyright (c) 2023
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define MAXINT INT_MAX

/**
 * @brief Creates an adjacency matrix for a disconnected graph.
 * 
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @return Number of edges in the graph.
 */
int createDisconnectedGraph(int N, int* adj_matrix);

/**
 * @brief Creates an adjacency matrix for a sparse graph.
 * 
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createSparseGraph(int N, int* adj_matrix, int seed);

/**
 * @brief Creates an adjacency matrix for a dense graph.
 * 
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createDenseGraph(int N, int* adj_matrix, int seed);

/**
 * @brief Creates an adjacency matrix for a fully connected graph.
 * 
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createFullConnectedGraph(int N, int* adj_matrix, int seed);

/**
 * @brief Allocates and returns a vector of length N*N.
 * 
 * @param N Number of vertices.
 * @return Pointer to the allocated vector.
 */
int* createGraph(int N);

/**
 * @brief Frees the memory allocated for the adjacency matrix.
 * 
 * @param adj_matrix The adjacency matrix to be freed.
 */
void freeGraph(int* adj_matrix);

/**
 * @brief Creates a graph based on the specified type.
 * 
 * @param type Type of the graph (0: Fully Connected, 1: Dense, 2: Intermediate Density, 3: Sparse).
 * @param numberOfVertices Number of vertices in the graph.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createGraphByType(int type, int numberOfVertices, int* adj_matrix, int seed);

#endif // GRAPH_H


