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
 * @file dijkstra_cuda.cu
 * @copyright Copyright (c) 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include "graph.c"
#include "utility.c"

#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "CUDA_SAFE_CALL: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

/**
 * @brief Finds the index of the minimum distance in the distances array that has not been visited.
 *
 * @param distances   Array of distances.
 * @param marker      Array indicating whether a node has been visited.
 * @param array_size  Size of the arrays.
 * @return            Index of the minimum distance.
 */
int searchMinIndex(int* distances, int* marker, int array_size) {
    int j;
    int local_min_pair[2];
    local_min_pair[0]=MAXINT;
    local_min_pair[1]=-1;

    #pragma omp parallel for private(j) shared(local_min_pair)
        for (j = 0; j < array_size; j++) {
            if (!marker[j] && distances[j] < local_min_pair[0] && distances[j] != MAXINT) {
                #pragma omp critical
                {
                    if (distances[j] < local_min_pair[0]) {
                        local_min_pair[0] = distances[j];
                        local_min_pair[1] = j;
                    }
                }
            }
        }
    marker[local_min_pair[1]] = 1;
    return local_min_pair[1];
}

/**
 * @brief CUDA kernel to update distances based on the Dijkstra algorithm.
 *
 * @param graph         Adjacency matrix of the graph.
 * @param node_dist     Array of distances from the source vertex.
 * @param path          Array representing the predecessor shortest path from the source node.
 * @param visited       Array indicating whether a node has been visited.
 * @param source        Source node to update distances for Dijkstra's algorithm.
 * @param num_vertices  Number of vertices in the graph.
 */
__global__ void cuda_update_distance(int* graph, int* node_dist, int* path, int* visited, int source, int num_vertices) {
    int node = blockIdx.x * blockDim.x + threadIdx.x;
    if(node < num_vertices){
        visited[source] = 1;
        int edge = graph[source * num_vertices + node];
        int new_dist = node_dist[source] + edge;
        if((edge != MAXINT) && (edge != 0)){
            if ((visited[node] != 1) && (new_dist < node_dist[node])) {
                node_dist[node] = new_dist;
                path[node] = source;
            }
        }
    }
}

/**
 * @brief Main function implementing the OMP/CUDA Dijkstra algorithm.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments.
 * @return     0 on successful execution, 1 on failure.
 */
int main(int argc, char* argv[]) {

    if (argc != 5) {
    printf("Usage: %s <number_of_vertices> <graph_type> <seed> <optimization_type>\n", argv[0]);
    exit(1);
    }

    struct timeval all_start, all_end, start, end, total_time, graph_creation_time, dijkstra_time, allocation_time;
    gettimeofday(&all_start, NULL); // the variable all_start is used to measure the run time of the program

    int N,n_edges;  //number of vertices,edges
    N = atoi(argv[1]);

    int numThreadsPerBlock;      // blockSize
    int minGridSize;    // grid min size
    int numBlocks;      //gridSize

    if (N <= 0){ // check if the number of vertices in input is coherent
        printf("Number of vertices is wrong!\n");
        exit(1);
    }
    if(atoi(argv[2])<0 || atoi(argv[2])>3){ // check if the type of graph is coherent
        printf("Type of graph is wrong!\n");
        exit(1);
    }
    if(atoi(argv[4])<0 || atoi(argv[4])>3){ // check if the optimization is coherent
        printf("Optimization choice is wrong!\n");
        exit(1);
    }

    //calculate blocksize optimally
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &numThreadsPerBlock, (void*)cuda_update_distance, 0, N);

    // Calculate gridSize & blockSize based on the size of the problem
    numBlocks = (N + numThreadsPerBlock - 1) / numThreadsPerBlock;
    dim3 gridSize(numBlocks, 1, 1);
    dim3 blockSize(numThreadsPerBlock, 1, 1);

    int type = atoi(argv[2]); // type is the type of graph to be created (e.g. sparse, dense, etc.)
    int seed = atoi(argv[3]); // seed is the random seed generator used to guarantee the same graph in the tests

    //host vectors
    int* graph;
    int* host_distances;
    int* host_path;
    int* host_visited;

    //device vectors
    int* device_graph;
    int* device_distances;
    int* device_path;
    int* device_visited;

    //start the timer to calculate the graph creation time
    gettimeofday(&start, NULL);

    //host allocations
    graph = createGraph(N);
    n_edges = createGraphByType(type, N, graph, seed);
    assert(graph != NULL); // check if the graph is created correctly
    gettimeofday(&end, NULL); // stop the timer and calculate the graph creation time
    timersub(&end, &start, &graph_creation_time);

    host_distances = (int*)malloc(N * sizeof(int));
    host_path = (int*)malloc(N * sizeof(int));
    host_visited = (int*)malloc(N * sizeof(int));

    // arrays inizializations
    initializeArray(host_distances, N, MAXINT);              //all node distances are infinity
    initializeArray(host_path, N, SOURCE_VERTEX);    //parent nodes are SOURCE (no parents yet)
    initializeArray(host_visited, N, 0);              //no nodes have been visited
    host_distances[SOURCE_VERTEX] = 0;                     //start distance for SOURCE_VERTEX is 0;
    host_path[SOURCE_VERTEX]= -1;                  //start partent for SOURCE_VERTEX is -1;

    //start the timer to calculate gpu allocation time
    gettimeofday(&start, NULL);

    //device allocations
    CUDA_SAFE_CALL(cudaMalloc((void**)&device_graph, N * N * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&device_distances, N * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&device_path, N * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&device_visited, N * sizeof(int)));
    //gpu source        cpu source      memory size     HtD or DtH
    CUDA_SAFE_CALL(cudaMemcpy(device_graph, graph, N * N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(device_distances, host_distances, N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(device_path, host_path, N * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(device_visited, host_visited, N * sizeof(int), cudaMemcpyHostToDevice));

    gettimeofday(&end, NULL); // stop the timer and calculate device allocation time
    timersub(&end, &start, &allocation_time);

    //start the timer to calculate dijkstra execution time
    gettimeofday(&start, NULL);

    for (int i = 0; i < N; i++) {
        //step 1: search min
        int min = searchMinIndex(host_distances, host_visited, N);

        //step 2: update distance
        cuda_update_distance<<<gridSize,blockSize>>>(device_graph, device_distances, device_path, device_visited, min, N);

        // transfers the updated distance vector from device to host
        cudaDeviceSynchronize();
        CUDA_SAFE_CALL(cudaMemcpy(host_distances, device_distances, N * sizeof(int), cudaMemcpyDeviceToHost));
    }
    // transfers the path vector from device to host
    cudaDeviceSynchronize();
    CUDA_SAFE_CALL(cudaMemcpy(host_path, device_path, N * sizeof(int), cudaMemcpyDeviceToHost));

    gettimeofday(&end, NULL); // stop the timer and calculate the Dijkstra execution time
    timersub(&end, &start, &dijkstra_time);

    //calculate all time program time
    gettimeofday(&all_end, NULL);
    timersub(&all_end, &all_start, &total_time);

    //print results
    FILE *fp;
    char filepath[200];
    sprintf(filepath, "ResultsOMP+CUDA/ResultsDijkstra/opt%d/type%d/OpenMP+CUDA_%d_%d_%d_%d.txt", atoi(argv[4]), type, N, n_edges, omp_get_max_threads(), numThreadsPerBlock); // define the path where to store the results
    char *filename = filepath;
    fp = fopen(filename, "w"); // open a file to write the results founded by the algorithm on it
    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", filename);
        return 1; // Return an error code
    }
    fprintf(fp, "Distance Vector: \n");
    print_vector_on_file(fp,N,host_distances);
    fprintf(fp, "------------------------------------------------------\n");
    fprintf(fp, "Path Vector: \n");
    print_vector_on_file(fp,N,host_path);
    fprintf(fp, "------------------------------------------------------\n");
    fclose(fp);

    // Open a file to write al the calculated times founded by the algorithm
    FILE *fp2;
    char filepath2[200];
    sprintf(filepath2, "ResultsOMP+CUDA/InfoTimeDijkstra/opt%d/type%d/%d_%d.csv", atoi(argv[4]), type, N, n_edges); // define the path where to store the times calculated for each phase.
    char *filename2 = filepath2;
    fp2 = fopen(filename2, "a+");
    if (fp2 == NULL) {
        fprintf(stderr, "Error opening file %s for writing\n", filename2);
        return 1; // Return an error code
    }
    fprintf(fp2, "OpenMp+CUDA;%d;%d; %ld.%06ld; %ld.%06ld; %ld.%06ld; %ld.%06ld;\n", omp_get_max_threads(),numThreadsPerBlock,(long int)total_time.tv_sec, (long int)total_time.tv_usec, (long int)graph_creation_time.tv_sec, (long int)graph_creation_time.tv_usec, (long int)allocation_time.tv_sec, (long int)allocation_time.tv_usec,(long int)dijkstra_time.tv_sec, (long int)dijkstra_time.tv_usec);
    fclose(fp2);

    // Free device memory
    CUDA_SAFE_CALL(cudaFree(device_graph));
    CUDA_SAFE_CALL(cudaFree(device_distances));
    CUDA_SAFE_CALL(cudaFree(device_path));
    CUDA_SAFE_CALL(cudaFree(device_visited));

    // Free host memory
    freeGraph(graph);
    free(host_distances);
    free(host_path);
    free(host_visited);
    return 0;
}
