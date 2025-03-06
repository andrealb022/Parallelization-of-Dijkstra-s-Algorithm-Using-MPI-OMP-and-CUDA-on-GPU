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
 * @file graph.c
 * @copyright Copyright (c) 2023
 */

#include "../Headers/graph.h"

/**
 * @brief Creates an adjacency matrix for a disconnected graph.
 *
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @return Number of edges in the graph.
 */
int createDisconnectedGraph(int N, int* adj_matrix) {
    // Initialize the adjacency matrix for a disconnected graph
    int i, j, Nedges = 0;
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (i == j) {
                adj_matrix[i*N + j] = 0;
            } else {
                adj_matrix[i*N + j] = MAXINT;
                adj_matrix[j*N + i] = MAXINT; // Add this line to maintain symmetry
            }
        }
    }
    return Nedges;
}

/**
 * @brief Creates an adjacency matrix for a graph with intermediate density.
 *
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createIntermediateDensityGraph(int N, int* adj_matrix, int seed) {
    // Initialize the adjacency matrix for a graph with few edges
    srand(seed);

    int i, j, Nedges = 0;
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (i == j) {
                adj_matrix[i*N + j] = 0;
            } else {
                int weight = rand() % 100 + 1;
                if (weight % 3 == 0) { // 30%
                    adj_matrix[i*N + j] = MAXINT;
                    adj_matrix[j*N + i] = MAXINT;
                } else {
                    adj_matrix[i*N + j] = weight;
                    adj_matrix[j*N + i] = weight; // Add this line to maintain symmetry
                    Nedges += 1;
                }
            }
        }
    }
    return Nedges;
}

/**
 * @brief Creates an adjacency matrix for a sparse graph.
 *
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createSparseGraph(int N, int* adj_matrix, int seed) {
    // Initialize the adjacency matrix for a graph with few edges
    srand(seed);

    int i, j, Nedges = 0;
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (i == j) {
                adj_matrix[i*N + j] = 0;
            } else {
                int weight = rand() % 100 + 1;
                if (weight % 2 == 0) { // 50%
                    adj_matrix[i*N + j] = MAXINT;
                    adj_matrix[j*N + i] = MAXINT;
                } else {
                    adj_matrix[i*N + j] = weight;
                    adj_matrix[j*N + i] = weight; // Add this line to maintain symmetry
                    Nedges += 1;
                }
            }
        }
    }
    return Nedges;
}

/**
 * @brief Creates an adjacency matrix for a dense graph.
 *
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createDenseGraph(int N, int* adj_matrix, int seed) {
    // Initialize the adjacency matrix for a graph with many edges
    srand(seed);

    int i, j, Nedges = 0;
    for (i = 0; i < N; ++i) {
        for (j = i; j < N; ++j) {
            if (i == j) {
                adj_matrix[i*N + j] = 0;
            } else {
                int weight = rand() % 100 + 1;
                if (weight % 10 == 0) { // 10%
                    adj_matrix[i*N + j] = MAXINT;
                    adj_matrix[j*N + i] = MAXINT;
                } else {
                    adj_matrix[i*N + j] = weight;
                    adj_matrix[j*N + i] = weight; // Add this line to maintain symmetry
                    Nedges += 1;
                }
            }
        }
    }
    return Nedges;
}

/**
 * @brief Creates an adjacency matrix for a fully connected graph.
 *
 * @param N Number of vertices.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createFullConnectedGraph(int N, int* adj_matrix, int seed) {
    // Create a fully connected graph
    srand(seed);

    int i, j, Nedges = 0;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            if (i == j) {
                adj_matrix[i*N + j] = 0;
            } else {
                adj_matrix[i*N + j] = rand() % 100 + 1;
                adj_matrix[j*N + i] = adj_matrix[i*N + j]; // Add this line to maintain symmetry
                Nedges += 1;
            }
        }
    }
    return Nedges;
}

/**
 * @brief Allocates and returns a vector of length N*N.
 *
 * @param N Number of vertices.
 * @return Pointer to the allocated vector.
 */
int* createGraph(int N) {
    int* graph = (int*)malloc(N * N * sizeof(int));
    if (graph == NULL) {
        // Gestisci l'errore di allocazione di memoria
        fprintf(stderr, "Error: impossible to allocate memory.\n");
        exit(EXIT_FAILURE);
    }
    return graph;
}


/**
 * @brief Frees the memory allocated for the adjacency matrix.
 *
 * @param adj_matrix The adjacency matrix to be freed.
 */
void freeGraph(int* adj_matrix) {
    // Free the memory allocated for the vector
    free(adj_matrix);
}

/**
 * @brief Creates a graph based on the specified type.
 *
 * @param type Type of the graph (0: Fully Connected, 1: Dense, 2: Intermediate Density, 3: Sparse).
 * @param numberOfVertices Number of vertices in the graph.
 * @param adj_matrix The adjacency matrix to be filled.
 * @param seed Random seed for reproducibility.
 * @return Number of edges in the graph.
 */
int createGraphByType(int type, int numberOfVertices, int* adj_matrix, int seed) {
    switch (type) {
        case 0:
            return createFullConnectedGraph(numberOfVertices, adj_matrix, seed);
            break;
        case 1:
            return createDenseGraph(numberOfVertices, adj_matrix, seed);
            break;
        case 2:
            return createIntermediateDensityGraph(numberOfVertices, adj_matrix, seed);
            break;
        case 3:
            return createSparseGraph(numberOfVertices, adj_matrix, seed);
            break;
        default:
            return -1;
            break;
    }
}
