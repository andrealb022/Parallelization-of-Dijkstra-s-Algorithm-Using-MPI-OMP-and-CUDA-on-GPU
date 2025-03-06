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
 * @file dijkstra_serial.c
 * @copyright Copyright (c) 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include "../Headers/graph.h"
#include "../Headers/utility.h"


/**
 * @brief Implements Dijkstra's algorithm for finding the shortest paths in a graph.
 *
 * @param source The source vertex for the shortest paths.
 * @param num_vertices The total number of vertices in the graph.
 * @param graph The adjacency matrix representing the graph.
 * @param distances The array to store the calculated distances from the source.
 * @param visited Array to mark visited vertices during the algorithm.
 * @param path Array to store the predecessors in the shortest paths.
 */
void dijkstra(int source, int num_vertices, int *graph, int *distances, int *visited, int *path) {
    int i, j;
    int tmp, current_vertex;

    // Initialization: Set the initial distances from source to all the other vertices
    for (i = 0; i < num_vertices; i++) {
        distances[i] = graph[source * num_vertices + i];
        visited[i] = FALSE;
        path[i] = source;
    }
    visited[source] = TRUE;
    path[source] = -1;

    // There are N-1 steps in the greedy algorithm:
    for (j = 1; j < num_vertices; j++) {

        /* Step 1: Find the local vertex that is at the smallest distance from source */
        tmp = MAXINT;
        current_vertex = -1;

        for (i = 0; i < num_vertices; i++) {
            if (!visited[i] && distances[i] < tmp && distances[i] != MAXINT) {
                current_vertex = i;
                tmp = distances[i];
            }
        }

        /* Step 2: mark the minimum vertex as being seen */
        if (current_vertex == -1)
            return;
        else
            visited[current_vertex] = TRUE;

        /* Step 3: Update the distances given that u got inserted */
        for (i = 0; i < num_vertices; i++) {
            if (graph[current_vertex * num_vertices + i] != MAXINT) {
                if (!visited[i] && ((distances[current_vertex] + graph[current_vertex * num_vertices + i]) < distances[i])) {
                    distances[i] = distances[current_vertex] + graph[current_vertex * num_vertices + i];
                    path[i] = current_vertex;
                }
            }
        }
    }
}

/**
 * @brief Main function for executing Dijkstra's algorithm sequentially.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return 0 on successful execution, 1 on failure.
 */
int main(int argc, char *argv[]) {
    struct timeval all_start, all_end, start, end, total_time, graph_creation_time, dijkstra_time;
    gettimeofday(&all_start, NULL); // Variable used to measure the program's runtime

    if (argc != 5) {
        printf("Usage: %s <number_of_vertices> <graph_type> <seed> <optimization_type>\n", argv[0]);
        exit(1);
    }
    if (atoi(argv[1]) <= 0) {
        printf("Number of vertices is wrong!\n");
        exit(1);
    }
    if (atoi(argv[2]) < 0 || atoi(argv[2]) > 3) {
        printf("Type of graph is wrong!\n");
        exit(1);
    }
    if (atoi(argv[4]) < 0 || atoi(argv[4]) > 3) {
        printf("Optimization choice is wrong!\n");
        exit(1);
    }

    FILE *fp;
    char filepath[200];
    int num_vertices = atoi(argv[1]);
    int type = atoi(argv[2]);
    int seed = atoi(argv[3]);
    int i, num_edges;
    int *adj_matrix; /* Adjacency matrix */
    int *distances;
    int *visited;
    int *path; // path[] is indicating where from we came to that node, the predecessor
    // allocate  visited, distance and path arrays
    distances = (int *)malloc(num_vertices * sizeof(int));
    visited = (int *)malloc(num_vertices * sizeof(int));
    path = (int *)malloc(num_vertices * sizeof(int));

    gettimeofday(&start, NULL); // Start the timer to measure the graph creation time

    //create  adj matrix
    adj_matrix = createGraph(num_vertices);
    num_edges = createGraphByType(type, num_vertices, adj_matrix, seed);
    assert(adj_matrix != NULL);

    gettimeofday(&end, NULL); // Stop the timer to measure the graph creation time
    timersub(&end, &start, &graph_creation_time);

    gettimeofday(&start, NULL); // Start the timer to measure the Dijkstra's algorithm time

    /*Implement the single source dijkstra's algorithm*/
    dijkstra(SOURCE, num_vertices, adj_matrix, distances, visited, path);

    gettimeofday(&end, NULL); // After all the computation and freeing the allocated memory, stop the timer and calculate the Dijkstra execution time
    timersub(&end, &start, &dijkstra_time);

    //calculate total program execution time
    gettimeofday(&all_end, NULL);
    timersub(&all_end, &all_start, &total_time);

    // PRINT RESULTS FOR OPENMP+MPI & OPENMP+CUDA VERSIONS
    char *filepaths[] = {
        "ResultsOMP+MPI/ResultsDijkstra/opt%d/type%d/Sequential%d_%d.txt",
        "ResultsOMP+CUDA/ResultsDijkstra/opt%d/type%d/Sequential%d_%d.txt",
        "ResultsOMP+MPI/InfoTimeDijkstra/opt%d/type%d/%d_%d.csv",
        "ResultsOMP+CUDA/InfoTimeDijkstra/opt%d/type%d/%d_%d.csv"
    };
    for(i=0;i<2;i++){
        sprintf(filepath, filepaths[i], atoi(argv[4]), type, num_vertices, num_edges); // define the path where to store the results
        char *filename = filepath;
        fp = fopen(filename, "w"); // open a file to write the results founded by the algorithm on it
        if (fp == NULL) {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            return 1; // Return an error code
        }

        fprintf(fp, "Distance Vector: \n");
        print_vector_on_file(fp,num_vertices,distances);
        fprintf(fp, "------------------------------------------------------\n");
        fprintf(fp, "Path Vector: \n");
        print_vector_on_file(fp,num_vertices,path);
        fprintf(fp, "------------------------------------------------------\n");
        fclose(fp);
    }

    // open a file to write al the calculated times founded by the algorithm
    for(i=2;i<4;i++){
        sprintf(filepath, filepaths[i], atoi(argv[4]), type, num_vertices, num_edges); // define the path where to store the times calculated for each phase, also for the communicational that in this case in 0.0 because is the sequential algorithm on just one processor
        char *filename = filepath;
        fp = fopen(filename, "a+");
        if (fp == NULL) {
            fprintf(stderr, "Error opening file %s for writing\n", filename);
            return 1; // Return an error code
        }
        fprintf(fp, "Sequential;0;0; %ld.%06ld; %ld.%06ld; 0.0; %ld.%06ld;\n", (long int)total_time.tv_sec, (long int)total_time.tv_usec, (long int)graph_creation_time.tv_sec, (long int)graph_creation_time.tv_usec, (long int)dijkstra_time.tv_sec, (long int)dijkstra_time.tv_usec);
        fclose(fp);
    }
    //free memory
    freeGraph(adj_matrix);
    free(distances);
    free(visited);
    free(path);

    return 0;
}
