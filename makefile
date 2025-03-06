.PHONY: all clean test compile0 compile1 compile2 compile3 test0 test1 test2 test3 testCuda0 testCuda1 testCuda2 testCuda3

all: directories compile0 compile1 compile2 compile3

clean:
	rm -rf ./Build/*
	rm -rf ./ResultsOMP+MPI
	rm -rf ./ResultsOMP+CUDA

test: test0 test1 test2 test3 testCuda0 testCuda1 testCuda2 testCuda3 statistics

types=3
it=5
mpiMax=4
ompMax=8
ompMaxCUDA=16
dimensions=3
seed=4632

directories:
	python3 ./Source/directories.py

compile0:
	gcc -c -o ./Build/utility.o ./Source/utility.c -O0
	gcc -c -o ./Build/graph.o ./Source/graph.c -O0
	gcc -c -o ./Build/dijkstra_serial.o ./Source/dijkstra_serial.c -O0
	mpicc -c -o ./Build/dijkstra_MpiOpenmp.o ./Source/dijkstra_MpiOpenmp.c -fopenmp -O0
	nvcc -c -o ./Build/dijkstra_cuda.o ./Source/dijkstra_cuda.cu -Xcompiler -fopenmp -O0
	
	gcc -o ./Build/dijkstra_serial0.exe ./Build/dijkstra_serial.o ./Build/utility.o ./Build/graph.o -O0
	mpicc -o ./Build/dijkstra_MpiOpenmp0.exe ./Build/dijkstra_MpiOpenmp.o ./Build/utility.o ./Build/graph.o -fopenmp -O0
	nvcc -o ./Build/dijkstra_cuda0.exe ./Build/dijkstra_cuda.o ./Build/utility.o ./Build/graph.o -Xcompiler -fopenmp -O0

compile1:
	gcc -c -o ./Build/utility.o ./Source/utility.c -O1
	gcc -c -o ./Build/graph.o ./Source/graph.c -O1
	gcc -c -o ./Build/dijkstra_serial.o ./Source/dijkstra_serial.c -O1
	mpicc -c -o ./Build/dijkstra_MpiOpenmp.o ./Source/dijkstra_MpiOpenmp.c -fopenmp -O1
	nvcc -c -o ./Build/dijkstra_cuda.o ./Source/dijkstra_cuda.cu -Xcompiler -fopenmp -O1

	gcc -o ./Build/dijkstra_serial1.exe ./Build/dijkstra_serial.o ./Build/utility.o ./Build/graph.o -O1
	mpicc -o ./Build/dijkstra_MpiOpenmp1.exe ./Build/dijkstra_MpiOpenmp.o ./Build/utility.o ./Build/graph.o -fopenmp -O1
	nvcc -o ./Build/dijkstra_cuda1.exe ./Build/dijkstra_cuda.o ./Build/utility.o ./Build/graph.o -Xcompiler -fopenmp -O1

compile2:
	gcc -c -o ./Build/utility.o ./Source/utility.c -O2
	gcc -c -o ./Build/graph.o ./Source/graph.c -O2
	gcc -c -o ./Build/dijkstra_serial.o ./Source/dijkstra_serial.c -O2
	mpicc -c -o ./Build/dijkstra_MpiOpenmp.o ./Source/dijkstra_MpiOpenmp.c -fopenmp -O2
	nvcc -c -o ./Build/dijkstra_cuda.o ./Source/dijkstra_cuda.cu -Xcompiler -fopenmp -O2
	
	gcc -o ./Build/dijkstra_serial2.exe ./Build/dijkstra_serial.o ./Build/utility.o ./Build/graph.o -O2
	mpicc -o ./Build/dijkstra_MpiOpenmp2.exe ./Build/dijkstra_MpiOpenmp.o ./Build/utility.o ./Build/graph.o -fopenmp -O2
	nvcc -o ./Build/dijkstra_cuda2.exe ./Build/dijkstra_cuda.o ./Build/utility.o ./Build/graph.o -Xcompiler -fopenmp -O2

compile3:
	gcc -c -o ./Build/utility.o ./Source/utility.c -O3
	gcc -c -o ./Build/graph.o ./Source/graph.c -O3
	gcc -c -o ./Build/dijkstra_serial.o ./Source/dijkstra_serial.c -O3
	mpicc -c -o ./Build/dijkstra_MpiOpenmp.o ./Source/dijkstra_MpiOpenmp.c -fopenmp -O3
	nvcc -c -o ./Build/dijkstra_cuda.o ./Source/dijkstra_cuda.cu -Xcompiler -fopenmp -O3
	
	gcc -o ./Build/dijkstra_serial3.exe ./Build/dijkstra_serial.o ./Build/utility.o ./Build/graph.o -O3
	mpicc -o ./Build/dijkstra_MpiOpenmp3.exe ./Build/dijkstra_MpiOpenmp.o ./Build/utility.o ./Build/graph.o -fopenmp -O3
	nvcc -o ./Build/dijkstra_cuda3.exe ./Build/dijkstra_cuda.o ./Build/utility.o ./Build/graph.o -Xcompiler -fopenmp -O3

statistics:
	python3 ./Source/analize.py
	python3 ./Source/analizeCUDA.py
	python3 ./Source/tableCreator.py
	python3 ./Source/plot.py
	python3 ./Source/plotCUDA.py

test0:
	opt=0 ; \
	nodes=1536 ; \
	max=10; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; mpi=1; i=1 ; while [ $$i -le $(it) ] ; do \
				./Build/dijkstra_serial0.exe $$nodes $$type $(seed) $$opt ; \
				while [ $$mpi -le $$(($(mpiMax))) ] ; do \
					while [ $$omp -le $(ompMax) ] ; do \
						export OMP_NUM_THREADS=$$omp && \
						if [ $$(($$mpi + $$omp)) -le $$max ]; then \
							mpirun -np $$mpi ./Build/dijkstra_MpiOpenmp0.exe $$nodes $$type $(seed) $$opt; \
						fi; \
						omp=$$((omp*2)); \
					done; \
					omp=1; \
					mpi=$$((mpi*2)); \
				done; \
				i=$$((i+1)); \
				mpi=1; \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done; 
test1:
	opt=1 ; \
	nodes=1536 ; \
	max=10 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; mpi=1; i=1 ; while [ $$i -le $(it) ] ; do \
				./Build/dijkstra_serial1.exe $$nodes $$type $(seed) $$opt ; \
				while [ $$mpi -le $$(($(mpiMax))) ] ; do \
					while [ $$omp -le $(ompMax) ] ; do \
						export OMP_NUM_THREADS=$$omp && \
						if [ $$(($$mpi + $$omp)) -le $$max ]; then \
							mpirun -np $$mpi ./Build/dijkstra_MpiOpenmp1.exe $$nodes $$type $(seed) $$opt; \
						fi; \
						omp=$$((omp*2)); \
					done; \
					omp=1; \
					mpi=$$((mpi*2)); \
				done; \
				i=$$((i+1)); \
				mpi=1; \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done; 
test2:
	opt=2 ; \
	nodes=1536 ; \
	max=10 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; mpi=1; i=1 ; while [ $$i -le $(it) ] ; do \
				./Build/dijkstra_serial2.exe $$nodes $$type $(seed) $$opt ; \
				while [ $$mpi -le $$(($(mpiMax))) ] ; do \
					while [ $$omp -le $(ompMax) ] ; do \
						export OMP_NUM_THREADS=$$omp && \
						if [ $$(($$mpi + $$omp)) -le $$max ]; then \
							mpirun -np $$mpi ./Build/dijkstra_MpiOpenmp2.exe $$nodes $$type $(seed) $$opt; \
						fi; \
						omp=$$((omp*2)); \
					done; \
					omp=1; \
					mpi=$$((mpi*2)); \
				done; \
				i=$$((i+1)); \
				mpi=1; \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done; 
test3:
	opt=3 ; \
	nodes=1536 ; \
	max=10 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; mpi=1; i=1 ; while [ $$i -le $(it) ] ; do \
				./Build/dijkstra_serial3.exe $$nodes $$type $(seed) $$opt ; \
				while [ $$mpi -le $$(($(mpiMax))) ] ; do \
					while [ $$omp -le $(ompMax) ] ; do \
						export OMP_NUM_THREADS=$$omp && \
						if [ $$(($$mpi + $$omp)) -le $$max ]; then \
							mpirun -np $$mpi ./Build/dijkstra_MpiOpenmp3.exe $$nodes $$type $(seed) $$opt; \
						fi; \
						omp=$$((omp*2)); \
					done; \
					omp=1; \
					mpi=$$((mpi*2)); \
				done; \
				i=$$((i+1)); \
				mpi=1; \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done; 

testCuda0:
	opt=0 ; \
	nodes=1536 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; i=1 ; while [ $$i -le $(it) ] ; do \
				while [ $$omp -le $(ompMaxCUDA) ] ; do \
					export OMP_NUM_THREADS=$$omp && \
					./Build/dijkstra_cuda0.exe $$nodes $$type $(seed) $$opt; \
					omp=$$((omp*2)); \
				done; \
				omp=1; \
				i=$$((i+1)); \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done;
	
testCuda1:
	opt=1 ; \
	nodes=1536 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; i=1 ; while [ $$i -le $(it) ] ; do \
				while [ $$omp -le $(ompMaxCUDA) ] ; do \
					export OMP_NUM_THREADS=$$omp && \
					./Build/dijkstra_cuda1.exe $$nodes $$type $(seed) $$opt; \
					omp=$$((omp*2)); \
				done; \
				omp=1; \
				i=$$((i+1)); \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done; 
	
testCuda2:
	opt=2 ; \
	nodes=1536 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; i=1 ; while [ $$i -le $(it) ] ; do \
				while [ $$omp -le $(ompMaxCUDA) ] ; do \
					export OMP_NUM_THREADS=$$omp && \
					./Build/dijkstra_cuda2.exe $$nodes $$type $(seed) $$opt; \
					omp=$$((omp*2)); \
				done; \
				omp=1; \
				i=$$((i+1)); \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done;  
	
testCuda3:
	opt=3 ; \
	nodes=1536 ; \
	dimension=1; while [ $$dimension -le $(dimensions) ] ; do \
		type=0; while [ $$type -le $$(($(types))) ] ; do \
			omp=1; i=1 ; while [ $$i -le $(it) ] ; do \
				while [ $$omp -le $(ompMaxCUDA) ] ; do \
					export OMP_NUM_THREADS=$$omp && \
					./Build/dijkstra_cuda3.exe $$nodes $$type $(seed) $$opt; \
					omp=$$((omp*2)); \
				done; \
				omp=1; \
				i=$$((i+1)); \
				omp=1; \
			done; \
			type=$$((type+1)); \
		done; \
		nodes=$$((nodes+6016)); \
		dimension=$$((dimension+1)); \
	done;


