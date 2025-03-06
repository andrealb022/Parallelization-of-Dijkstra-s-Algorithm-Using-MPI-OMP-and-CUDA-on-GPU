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
 # @file tableCreator.py
 # @copyright Copyright (c) 2023

import matplotlib.pyplot as plt
import os
import shutil

def createTable(sourceFile, outputPath, outputName):
    with open(sourceFile, "r") as f:
        # Read the header of the file to get the column names
        col = f.readline().strip().split(";")
        col.pop(0)  # Remove the "Modality" column

        # Initialize lists for table rows and cells
        row = []
        cell = []

        # Read the rest of the lines from the file
        lines = f.readlines()

        # For each line in the file
        for line in lines:
            # Split the line into elements separated by semicolons
            elem = line.strip().split(";")
            # Extract the modality (first element of the line)
            modality = elem.pop(0)
            
            # Add the modality to the list of rows
            row.append(modality)

            # Add the values of the row to the list of cells
            new_row = [elem.pop(0), elem.pop(0)]
            for val in elem:
                new_row.append(str(format(float(val), ".7f")))
            cell.append(new_row)

    # Create the table using Matplotlib
    fig, ax = plt.subplots()
    ax.set_axis_off()
    table = ax.table(
        cellText=cell,
        rowLabels=row,
        colLabels=col,
        rowColours=["#d0d0d0"] * len(row),
        colColours=["#d0d0d0"] * len(col),
        cellLoc='center',
        loc='upper center',
        rowLoc="center",
        colLoc="center"
    )

    # Set the dimensions of the table
    table.scale(0.95,0.95)
    table.auto_set_font_size(False)
    table.set_fontsize(3.6)

    # Set the title of the table
    name = outputName.split("_")
    title = f"{name[0]} nodes & {name[1]} edges"
    ax.set_title(title, fontweight="bold")

    # Save the image of the table
    plt.savefig(os.path.join(outputPath, outputName) + ".png", dpi=300, bbox_inches='tight')
    plt.close()

if os.path.exists("ResultsOMP+MPI/Tables"):
    shutil.rmtree('ResultsOMP+MPI/Tables')

os.mkdir("ResultsOMP+MPI/Tables")

for path, currentDirectory, files in os.walk("ResultsOMP+MPI/FinalAnalysisDijkstra"):
    for file in files:
        resultPath = path.replace("ResultsOMP+MPI/FinalAnalysisDijkstra", "ResultsOMP+MPI/Tables")
        if not (os.path.exists(resultPath)):
            os.makedirs(resultPath)
        createTable(os.path.join(path, file), resultPath, file.replace(".csv", ""))


if os.path.exists("ResultsOMP+CUDA/Tables"):
    shutil.rmtree('ResultsOMP+CUDA/Tables')

os.mkdir("ResultsOMP+CUDA/Tables")

for path, currentDirectory, files in os.walk("ResultsOMP+CUDA/FinalAnalysisDijkstra"):
    for file in files:
        resultPath = path.replace("ResultsOMP+CUDA/FinalAnalysisDijkstra", "ResultsOMP+CUDA/Tables")
        if not (os.path.exists(resultPath)):
            os.makedirs(resultPath)
        createTable(os.path.join(path, file), resultPath, file.replace(".csv", ""))
