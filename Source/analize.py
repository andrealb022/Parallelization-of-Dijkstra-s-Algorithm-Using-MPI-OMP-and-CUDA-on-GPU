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
 # @file analize.py
 # @copyright Copyright (c) 2023

import os
import shutil

def analizeFile(sourceFile,resultFile):#analize the source file and save the information in the result file
    f=open(sourceFile,"r")#open file and read all its lines
    lines=f.readlines()
    data={}
    num={}
    for line in lines:#for each line get the informations
        type,omp,mpi,totalTime,creationTime,comunicationTime,executionTime=line.strip()[:-1].split(";")
        if type not in data:
            data[type]={}
            num[type]={}
        if omp not in data[type]:
            data[type][omp]={}
            num[type][omp]={}
        if mpi not in data[type][omp]:
            data[type][omp][mpi]=[float(totalTime),float(creationTime),float(comunicationTime),float(executionTime)]
            num[type][omp][mpi]=1
        else:
            data[type][omp][mpi][0]+=float(totalTime)
            data[type][omp][mpi][1]+=float(creationTime)
            data[type][omp][mpi][2]+=float(comunicationTime)
            data[type][omp][mpi][3]+=float(executionTime)
            num[type][omp][mpi]+=1
    f.close()#close the source file and open the result file
    f=open(resultFile,"w")
    f.write("Modality;OMP;MPI;total time;graph creation time;comunication time;Dijkstra execution time;Speedup;Efficiency \n")
    for typeKey in data.keys():#for each element in the dictionary calculate the mean and save it in the result file
        ompList=list(data[typeKey].keys())
        ompList.sort(key=lambda x: (int(x), x))
        for ompKey in ompList:
            mpiList=list(data[typeKey][ompKey].keys())
            mpiList.sort(key=lambda x: (int(x), x))
            for mpiKey in mpiList:
                it=num[typeKey][ompKey][mpiKey]
                data[typeKey][ompKey][mpiKey][0]/=it
                data[typeKey][ompKey][mpiKey][1]/=it
                data[typeKey][ompKey][mpiKey][2]/=it
                data[typeKey][ompKey][mpiKey][3]/=it
                
                if typeKey =="Sequential":
                    # speedup
                    seq_time = data[typeKey][ompKey][mpiKey][0]
                    data[typeKey][ompKey][mpiKey].extend([1,1])
                else:
                    par_time = seq_time / data[typeKey][ompKey][mpiKey][0]
                    data[typeKey][ompKey][mpiKey].extend([par_time,par_time / ( int(mpiKey) * int(ompKey) )])
                
                
                f.write(typeKey+";"+str(ompKey)+";"+str(mpiKey)+";"+str(data[typeKey][ompKey][mpiKey][0])+
                ";"+str(data[typeKey][ompKey][mpiKey][1])+";"+str(data[typeKey][ompKey][mpiKey][2])+
                ";"+str(data[typeKey][ompKey][mpiKey][3])+";"+str(data[typeKey][ompKey][mpiKey][4]) + ";"+str(data[typeKey][ompKey][mpiKey][5]) + "\n")
    f.close()#close the result file

if os.path.exists("ResultsOMP+MPI/FinalAnalysisDijkstra"):#if Result directory exist, delete it and all its elements
    shutil.rmtree('ResultsOMP+MPI/FinalAnalysisDijkstra')

os.mkdir("ResultsOMP+MPI/FinalAnalysisDijkstra")#create the result directory

for path, currentDirectory, files in os.walk("ResultsOMP+MPI/InfoTimeDijkstra"):#for each file in the directory Informations, analize them and save the results in the relative path in Results directory 
    for file in files:
        resultPath=path.replace("ResultsOMP+MPI/InfoTimeDijkstra","ResultsOMP+MPI/FinalAnalysisDijkstra")
        if not(os.path.exists(resultPath)):
            os.makedirs(resultPath)
        analizeFile(os.path.join(path, file),os.path.join(resultPath,file))
