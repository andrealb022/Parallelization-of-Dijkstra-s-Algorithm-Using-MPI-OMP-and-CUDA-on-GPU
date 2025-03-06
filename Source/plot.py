 # Course: High Performance Computing 2023/2024 
 # 
 # Lecturer: Francesco Moscato	fmoscato@unisa.it 
 #
 # Student : 
 # Alberti Andrea	0622702370	a.alberti2@studenti.unisa.it
 #
 # 
 # Copyright (C) 2023 - All Rights Reserved 
 #
 # This file is part of DijkstraFinalProjectHPC 
 #
 # DijkstraFinalProjectHPC   is free software: you can redistribute it and/or modify 
 # it under the terms of the GNU General Public License as published by 
 # the Free Software Foundation, either version 3 of the License, or 
 # (at your option) any later version. 
 #
 # DijkstraFinalProjectHPC   is distributed in the hope that it will be useful, 
 # but WITHOUT ANY WARRANTY; without even the implied warranty of 
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 # GNU General Public License for more details. 
 #
 # You should have received a copy of the GNU General Public License 
 # along with DijkstraFinalProjectHPC . If not, see <http://www.gnu.org/licenses/>. 
 #
 #                                  REQUIREMENTS OF THE ASSIGNMENT
 #
 # Student shall provide a parallel version of Dijkstra algorithm with both "OpenMP + MPI" and "OpenMP + Cuda" approaches,
 # comparing results with a known solution on single-processing node. 
 # Results and differences shall be discussed for different inputs (type and size).
 # The parallel algorithm used in "OpenMP + MPI" solution could not be the same of the "OpenMP + CUDA" approach.
 #
 # @file plot.py
 # @copyright Copyright (c) 2023

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
def makePlot(sourceFile,outputPath,outputName):
    data = read_data_from_file(sourceFile)
    num_processes = [0, 1, 2, 4]
    omp_threads = [1, 2, 4, 8]

    # Creazione del grafico
    plt.figure(figsize=(10, 6))

    # Aggiungi curve per ciascuna configurazione parallela
    for i in range(len(omp_threads)):
        plt.plot(num_processes, data[i], label=f'{omp_threads[i]} OMP Threads')

    plt.plot(num_processes, num_processes, linestyle='--', color='black', label='Ideal Speedup')
    # Personalizza il grafico
    plt.xlabel('Number of Processes')
    plt.ylabel('Speedup')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(outputPath, outputName)+".png", bbox_inches='tight')
    plt.close()
    
def read_data_from_file(sourceFile):
    file=open(sourceFile,"r")
    lines = file.readlines()

    data = []
    for line in lines[1:]:
        values = line.strip().split(';')
        if len(values) == 9:  # Ensure that the line has the expected number of values
            data.append({
                float(values[7]), 
            })
        else:
            print(f"Error processing line: {line}")

    data.remove(data[0]) 
    values = [list(d)[0] for d in data] 
    while len(values) % 4 != 0:
        values.append(0.0)
    tmp = [values[i:i+3] for i in range(0,len(data),3)] 
    speedup = [[0] + sublist for sublist in tmp] 
    file.close() 
    return speedup
   
if os.path.exists("ResultsOMP+MPI/Plots"):#if Plots directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI/Plots')

os.mkdir("ResultsOMP+MPI/Plots")#create the plots directory

for path, currentDirectory, files in os.walk("ResultsOMP+MPI/FinalAnalysisDijkstra"):
    for file in files:
        resultPath=path.replace("ResultsOMP+MPI/FinalAnalysisDijkstra","ResultsOMP+MPI/Plots")
        if not(os.path.exists(resultPath)):
            os.makedirs(resultPath)
        makePlot(os.path.join(path, file),resultPath,file.replace(".csv",""))

