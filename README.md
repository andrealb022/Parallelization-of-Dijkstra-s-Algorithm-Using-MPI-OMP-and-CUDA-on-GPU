# DijkstraFinalProjectHPC

## Dependencies

* CMake
* MPICH
* OPENMP
* NVCC
* Python3

## How to run
Code for generating directories, tables and plots requires python interpreter and matplotlib library, the last to be installed with the command pip3 install matplotlib.

For running project programs, the following steps should be executed:
1.	Navigate to the folder containing the makefile
2.	To clear previously obtained achievements and previous builds, enter the command: make clean
3.	To generate the necessary directories and compile and linking the various source codes, enter the command: make all
4.	To run the algorithm for making tests, producing results, measurements, graphs and tables, enter the command: make test
5.	To clear previously obtained achievements and previous builds, enter the command: make clean

Results are found in the folders "OPENMP+MPI" & "OPENMP+CUDA", organized by version.
Results of the algorithms can be viewed in the "ResultDijkstra" folder, organized by optimization and type of graph and size.
Execution times of the algorithms and their average values can be viewed respectively in the "InfoTimeDijkstra" and "FinalAnalysisDijkstra" folders, organized by optimization and type of graph and size.
Results in graphical and tabular format can be viewed respectively in the "Plots" and "Tables" folders, organized by optimization and type of graph and size.

NB: Due to high execution time, the number of trials has been reduced by setting the corresponding variable "it"=5 in the makefile.


