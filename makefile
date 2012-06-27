INSTALL_DIR=/usr/local/
CC=g++

libmultigrid.so: grid.o multigrid.o mesh.o
	$(CC) -c -fPIC grid.cpp multigrid.cpp mesh.cpp
	$(CC) -shared -o libmultigrid.so grid.o multigrid.o mesh.o

install: libmultigrid.so
	cp grid.h multigrid.h mesh.h $(INSTALL_DIR)/include/
	cp libmultigrid.so $(INSTALL_DIR)/lib/


example: grid.o multigrid.o mesh.o main.o
	$(CC) main.o grid.o multigrid.o mesh.o -o example.out

example_lib: main.cpp  
	$(CC) -I $(INSTALL_DIR)/include/ -L $(INSTALL_DIR)/lib/ -o example.out main.cpp -lmultigrid

grid.o: grid.cpp grid.h
	$(CC) -c grid.cpp 
multigrid.o: multigrid.cpp multigrid.h
	$(CC) -c multigrid.cpp 
mesh.o: mesh.cpp mesh.h
	$(CC) -c mesh.cpp 
main.o: main.cpp grid.o mesh.o grid.o
	$(CC) -c main.cpp 
clean:
	rm -rf *.o *.out output/ *.so
