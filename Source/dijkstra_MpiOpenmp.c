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
 * @file dijkstra_MpiOpenmp.c
 * @copyright Copyright (c) 2023
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include <omp.h>
#include <assert.h>
#include "../Headers/graph.h"
#include "../Headers/utility.h"
/**
 * @brief Checks if the given number of MPI processes is valid for the specified number of vertices.
 *
 * @param n_vertices Number of vertices in the graph.
 * @param nprocess Number of MPI processes.
 * @return TRUE if the number of MPI processes is valid, FALSE otherwise.
 */
int checkNprocess(int n_vertices, int nprocess) {
    if (nprocess > 0)
        if (n_vertices % nprocess == 0)
            return TRUE;
    return FALSE;
}

/**
 * @brief Performs Dijkstra's algorithm to find the shortest paths from a source vertex to all other vertices in a graph.
 *
 * @param n_vertices Number of vertices in the graph.
 * @param source Source vertex for the shortest paths.
 * @param local_adj_matrix Weighted adjacency matrix representing the local portion of the graph.
 * @param distance Array to store the calculated shortest distances.
 * @param path Array to store the shortest paths.
 * @param communicator MPI communicator.
 */
void Dijkstra(int n_vertices, int source, int *local_adj_matrix, int *distance, int *path, MPI_Comm communicator) {
    int i, j, current_vertex;
    int n_vertices_local;   /* The number of vertices stored locally */
    int *vertex_marker;  /* Used to mark the vertices belonging to Vo */
    int first_vertex; /* The index number of the first vertex that is stored locally */
    int last_vertex;  /* The index number of the last vertex that is stored locally */
    int minimum_vertex, minimum_distance;
    int local_min_pair[2], global_min_pair[2];

    int num_processes, rank;
    MPI_Comm_size(communicator, &num_processes);
    MPI_Comm_rank(communicator, &rank);

    n_vertices_local = n_vertices / num_processes;
    first_vertex = rank * n_vertices_local;
    last_vertex = first_vertex + n_vertices_local - 1;
    vertex_marker = (int *)malloc(n_vertices_local * sizeof(int));

    /* Set the initial distances from source to all the other vertices */
    for (j = 0; j < n_vertices_local; j++) {
        distance[j] = local_adj_matrix[j * n_vertices + source];
        vertex_marker[j] = FALSE;
        path[j] = source;
    }

    /* The process that stores the source vertex marks it as being seen */
    if (source >= first_vertex && source <= last_vertex) {
        vertex_marker[source - first_vertex] = TRUE;
        path[source - first_vertex] = -1;
    }

    /* The main loop of Dijkstra's algorithm */
    for (i = 1; i < n_vertices; i++) {

        /* Step 1: Find the local vertex that is at the smallest distance from source */
        local_min_pair[0] = MAXINT; /*large number */
        local_min_pair[1] = -1;

        #pragma omp parallel for private(j) shared(local_min_pair)
        for (j = 0; j < n_vertices_local; j++) {
            if (!vertex_marker[j] && distance[j] < local_min_pair[0] && distance[j] != MAXINT) {
                #pragma omp critical
                {
                    if (distance[j] < local_min_pair[0]) {
                        local_min_pair[0] = distance[j];
                        local_min_pair[1] = first_vertex + j;
                    }
                }
            }
        }

        /* Step 2: Compute the global minimum vertex*/
        MPI_Allreduce(local_min_pair, global_min_pair, 1, MPI_2INT, MPI_MINLOC, communicator);
        minimum_distance = global_min_pair[0];
        minimum_vertex = global_min_pair[1];
        if (minimum_vertex == -1)
            return;
        else {
            /* The process that stores the minimum vertex marks it as being seen */
            if (minimum_vertex == local_min_pair[1]) {
                vertex_marker[minimum_vertex - first_vertex] = TRUE;
            }
        }

        /* Step 3: Update the distances given that u got inserted */
        #pragma omp parallel for
        for (j = 0; j < n_vertices_local; j++) {
            if (local_adj_matrix[j * n_vertices + minimum_vertex] != MAXINT) {
                if (!vertex_marker[j] && ((minimum_distance + local_adj_matrix[j * n_vertices + minimum_vertex]) < distance[j])) {
                    distance[j] = minimum_distance + local_adj_matrix[j * n_vertices + minimum_vertex];
                    path[j] = minimum_vertex;
                }
            }
        }
    }
    free(vertex_marker);
}

/**
 * @brief Main function to execute the MPI+OpenMP version of Dijkstra's algorithm.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return 0 on successful execution, 1 on failure.
 */
int main(int argc, char *argv[]){

    if (argc != 5) {
    printf("Usage: %s <number_of_vertices> <graph_type> <seed> <optimization_type>\n", argv[0]);
    exit(1);
    }

    struct timeval all_start, all_end, start, end, total_time, graph_creation_time, dijkstra_time, communication_time = {0, 0};
    gettimeofday(&all_start, NULL); // the variable all_start is used to measure the run time of the program

    int nprocess,myrank,N;
    N = atoi(argv[1]);

    /* Initialize MPI and get system information */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (!checkNprocess(N,nprocess)) // check if the MPI_processes number is coherent
    {
        printf("Number of MPI process is wrong!\n");
        exit(1);
    }
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

    int type = atoi(argv[2]); // type is the type of graph to be created (e.g. sparse, dense, etc.)
    int seed = atoi(argv[3]); // seed is the random seed generator used to guarantee the same graph in the tests
    int* adj_matrix; /*adjacency matrix*/
    int* distance;    /*distance vector*/
    int* path;        /*path vector*/
    int* local_adj_matrix;   /*local weight array*/
    int* localDistance; /*local distance vector*/
    int* localPath;     /*path distance vector*/
    int i, j, k, nlocal, n_edges;

    distance = (int *)malloc(N * sizeof(int));
    path = (int *)malloc(N * sizeof(int));
    /* Compute the number of elements to be stored locally. */
    nlocal = N / nprocess;

    /*allocate local adj matrix, local distance and local path arrays for each prosess*/
    local_adj_matrix = (int *)malloc(nlocal * N * sizeof(int));
    localDistance = (int *)malloc(nlocal * sizeof(int));
    localPath = (int *)malloc(nlocal * sizeof(int));

    //read matrix
    if (myrank == SOURCE)
    {
        gettimeofday(&start, NULL); // start the timer to measure the graph creation time

        adj_matrix = createGraph(N);
        n_edges = createGraphByType(type, N, adj_matrix, seed);
        assert(adj_matrix != NULL); // check if the graph is created correctly

        gettimeofday(&end, NULL); // stop the timer to measure the graph creation time
        timersub(&end, &start, &graph_creation_time);
    }

    if(nprocess == 1){
        gettimeofday(&start, NULL);

        /*Implement the single source dijkstra's algorithm*/
        Dijkstra(N, SOURCE_VERTEX, adj_matrix, distance, path, MPI_COMM_WORLD);
    }
    else{
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&start, NULL);
        /*distribute data dividing for each process nlocal rows of the adj_matrix*/
        MPI_Scatter(adj_matrix, nlocal * N, MPI_INT, local_adj_matrix, nlocal * N, MPI_INT, SOURCE,
                    MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        gettimeofday(&end, NULL); // stop the timer to measure the communication time
        timersub(&end, &start, &communication_time);

        gettimeofday(&start, NULL);

        /*Implement the single source dijkstra's algorithm*/
        Dijkstra(N, SOURCE_VERTEX, local_adj_matrix, localDistance, localPath, MPI_COMM_WORLD);

        /*collect local distance vector at the source process*/
        MPI_Gather(localDistance, nlocal, MPI_INT, distance, nlocal, MPI_INT, SOURCE,
                MPI_COMM_WORLD);
        MPI_Gather(localPath, nlocal, MPI_INT, path, nlocal, MPI_INT, SOURCE,
                MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    gettimeofday(&end, NULL); // stop the timer to measure dijkstra time and total_time
    timersub(&end, &start, &dijkstra_time);

    //calculate total program execution time
    gettimeofday(&all_end, NULL);
    timersub(&all_end, &all_start, &total_time);

    //print results
    if (myrank == SOURCE)
    {
        FILE *fp;
        char filepath[200];
        sprintf(filepath, "ResultsOMP+MPI/ResultsDijkstra/opt%d/type%d/Mpi+OpenMP_%d_%d_%d_%d.txt", atoi(argv[4]), type, N, n_edges, nprocess, omp_get_max_threads()); // define the path where to store the results
        char *filename = filepath;
        fp = fopen(filename, "w"); // open a file to write the results founded by the algorithm on it
        if (fp == NULL) {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            return 1; // Return an error code
        }

        //print results
        fprintf(fp, "Distance Vector: \n");
        print_vector_on_file(fp,N,distance);
        fprintf(fp, "------------------------------------------------------\n");
        fprintf(fp, "Path Vector: \n");
        print_vector_on_file(fp,N,path);
        fprintf(fp, "------------------------------------------------------\n");
        fclose(fp);
    }

    // open a file to write al the calculated times founded by the algorithm
    if(myrank == SOURCE){
        freeGraph(adj_matrix);
        FILE *fp2;
        char filepath2[200];
        sprintf(filepath2, "ResultsOMP+MPI/InfoTimeDijkstra/opt%d/type%d/%d_%d.csv", atoi(argv[4]), type, N, n_edges);
        char *filename2 = filepath2;
        fp2 = fopen(filename2, "a+");
        if (fp2 == NULL)
        {
            fprintf(stderr, "Error opening file %s for writing\n", filename2);
            return 1; // Return an error code
        }
        fprintf(fp2, "OpenMP+MPI;%d;%d; %ld.%06ld; %ld.%06ld; %ld.%06ld; %ld.%06ld;\n", omp_get_max_threads(), nprocess, (long int)total_time.tv_sec, (long int)total_time.tv_usec, (long int)graph_creation_time.tv_sec, (long int)graph_creation_time.tv_usec, (long int)communication_time.tv_sec, (long int)communication_time.tv_usec, (long int)dijkstra_time.tv_sec, (long int)dijkstra_time.tv_usec);
        fclose(fp2);
    }

    //free memory
    free(local_adj_matrix);
    free(localDistance);
    free(localPath);
    free(path);
    free(distance);

    MPI_Finalize();
    return 0;
}
