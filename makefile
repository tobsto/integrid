CC=mpiCC
#CC=mpiCC -O3
euo.out: grid.o multigrid.o mesh.o main.o
	$(CC) main.o grid.o multigrid.o mesh.o -o euo.out -lm -lboost_mpi -lboost_serialization
grid.o: grid.cpp grid.h
	$(CC) -c grid.cpp 
multigrid.o: multigrid.cpp multigrid.h
	$(CC) -c multigrid.cpp 
mesh.o: mesh.cpp mesh.h
	$(CC) -c mesh.cpp 
main.o: main.cpp multigrid.o grid.o mesh.o
	$(CC) -c main.cpp 
clean:
	rm *.o *.out
