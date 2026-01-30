CC = g++

CFLAGS = -std=gnu++17 #-Wall -Wextra -Iinclude

kalman: main.cpp kalman.hpp
	$(CC) $(CFLAGS) -Iinclude $^ -o kalman


bonus: main.cpp kalman.hpp
	$(CC) $(CFLAGS) -O3 -Iinclude  -DBONUS $^ -o kalman


speed: main.cpp kalman.hpp
	$(CC) $(CFLAGS) -O3  -Iinclude  $^ -o kalman


graph:
	cd graph_vtk && cmake . && make
	cp graph_vtk/graph_vtk ./graph


clean:
	rm -f kalman graph *.o *.txt log.txt variances.txt core*
	rm -rf graph_vtk/CMakeFiles graph_vtk/CMakeCache.txt graph_vtk/Makefile graph_vtk/cmake_install.cmake
	cd graph_vtk && make clean


re: clean kalman

.PHONY: 
