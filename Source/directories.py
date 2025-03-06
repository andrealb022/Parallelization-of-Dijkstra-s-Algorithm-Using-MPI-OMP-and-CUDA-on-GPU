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
 # @file directories.py
 # @copyright Copyright (c) 2023

import os
import shutil

if os.path.exists("ResultsOMP+MPI"):#if Informations directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI')

os.mkdir("ResultsOMP+MPI")#create the Informations directory

if os.path.exists("ResultsOMP+MPI/InfoTimeDijkstra"):#if Informations directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI/InfoTimeDijkstra')

os.mkdir("ResultsOMP+MPI/InfoTimeDijkstra")#create the Informations directory

if os.path.exists("ResultsOMP+MPI/FinalAnalysisDijkstra"):#if Results directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI/FinalAnalysisDijkstra')
 
os.mkdir("ResultsOMP+MPI/FinalAnalysisDijkstra")#create the Results directory

if os.path.exists("ResultsOMP+MPI/ResultsDijkstra"):#if ResultSCC directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI/ResultsDijkstra')
 
os.mkdir("ResultsOMP+MPI/ResultsDijkstra")#create the ResultsSCC directory


os.chdir("ResultsOMP+MPI/InfoTimeDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)

    os.chdir("..")

os.chdir("../FinalAnalysisDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)
    os.chdir("..")

os.chdir("../ResultsDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)
    os.chdir("..")

os.chdir("../..")
if os.path.exists("ResultsOMP+CUDA"):#if Informations directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+CUDA')

os.mkdir("ResultsOMP+CUDA")#create the Informations directory

if os.path.exists("ResultsOMP+CUDA/InfoTimeDijkstra"):#if Informations directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+CUDA/InfoTimeDijkstra')

os.mkdir("ResultsOMP+CUDA/InfoTimeDijkstra")#create the Informations directory

if os.path.exists("ResultsOMP+CUDA/FinalAnalysisDijkstra"):#if Results directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+CUDA/FinalAnalysisDijkstra')
 
os.mkdir("ResultsOMP+CUDA/FinalAnalysisDijkstra")#create the Results directory

if os.path.exists("ResultsOMP+CUDA/ResultsDijkstra"):#if ResultSCC directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+CUDA/ResultsDijkstra')
 
os.mkdir("ResultsOMP+CUDA/ResultsDijkstra")#create the ResultsSCC directory


os.chdir("ResultsOMP+CUDA/InfoTimeDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)

    os.chdir("..")

os.chdir("../FinalAnalysisDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)
    os.chdir("..")

os.chdir("../ResultsDijkstra") 

for i in range(4):
    opt="opt"+str(i)
    os.mkdir(opt)
    os.chdir(opt)
    for j in range(4):
        type="type"+str(j)
        os.mkdir(type)
    os.chdir("..")
