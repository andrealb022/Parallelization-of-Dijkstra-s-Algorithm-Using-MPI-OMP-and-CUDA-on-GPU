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
 # @file plotCUDA.py
 # @copyright Copyright (c) 2023

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
def makePlot(sourceFile,outputPath,outputName):
    data = read_data_from_file(sourceFile)
    omp_threads = [0, 1, 2, 4, 8, 16]
    ax = [0, 1, 2, 4, 6, 8, 10, 12, 14, 16]
    #ax = range(0,17)
    # Creazione del grafico
    plt.figure(figsize=(10, 6))

    # Aggiungi curve per ciascuna configurazione parallela
    plt.plot(omp_threads, data, marker='o', label='768 Threads per block')


    plt.plot(omp_threads, omp_threads, linestyle='--', color='black', label='Ideal Speedup')
    
    # Personalizza il grafico
    plt.xlabel('OMP Threads')
    plt.ylabel('Speedup')
    # Imposta gli assi x in modo esplicito
    plt.xticks(ax)
    plt.yticks(ax)
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(outputPath, outputName)+".png", bbox_inches='tight')
    plt.close()
    
def read_data_from_file(sourceFile):
    file = open(sourceFile, "r")
    lines = file.readlines()

    data = []
    for line in lines[1:]:
        values = line.strip().split(';')
        if len(values) == 9:
            data.append(float(values[7]))
        else:
            print(f"Error processing line: {line}")

    # Rimuovi il primo elemento solo se la lista non è vuota
    if data:
        data.pop(0)
    data.insert(0,0.0)
    return data
   
if os.path.exists("ResultsOMP+CUDA/Plots"):#if Plots directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+CUDA/Plots')

os.mkdir("ResultsOMP+CUDA/Plots")#create the plots directory

for path, currentDirectory, files in os.walk("ResultsOMP+CUDA/FinalAnalysisDijkstra"):
    for file in files:
        resultPath=path.replace("ResultsOMP+CUDA/FinalAnalysisDijkstra","ResultsOMP+CUDA/Plots")
        if not(os.path.exists(resultPath)):
            os.makedirs(resultPath)
        makePlot(os.path.join(path, file),resultPath,file.replace(".csv",""))

