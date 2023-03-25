MPICXX = mpicxx
CC = g++
CFLAGS = -g -O3
SRC = tsp.cpp tsp-omp.cpp tsp-mpi.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o # Add algorithms.o as a dependency
INCLUDES = -I/usr/include/x86_64-linux-gnu/mpich/

all: tsp tsp-omp tsp-mpi

tsp: tsp.o algorithms.o # Add algorithms.o as a dependency
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp.o

tsp-omp: tsp-omp.o algorithms.o # Add algorithms.o as a dependency
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-omp.o

tsp-mpi: tsp-mpi.o algorithms.o # Add algorithms.o as a dependency
	$(MPICXX) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-mpi.o

tsp.o: tsp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

tsp-omp.o: tsp-omp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

tsp-mpi.o: tsp-mpi.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp # Add a rule to compile algorithms.cpp to algorithms.o
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

.PHONY: all

clean:
	rm -f tsp tsp-omp tsp-mpi $(OBJ) algorithms.o