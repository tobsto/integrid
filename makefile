# ****************************************************************************
# This file is part of Multigrid.
#
# Multigrid is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Multigrid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Multigrid.  If not, see <http:#www.gnu.org/licenses/>.
#
# Copyright 2012 Tobias Stollenwerk
# ****************************************************************************

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
