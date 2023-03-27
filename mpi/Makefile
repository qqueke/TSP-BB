CC = g++
MPICXX = mpicxx
CFLAGS = -g -O3
SRC = tsp-mpi.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o 
INCLUDES = -I/usr/include/x86_64-linux-gnu/mpich/

tsp-mpi: tsp-mpi.o algorithms.o 
	$(MPICXX) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-mpi.o -lmpi

tsp-mpi.o: tsp-mpi.cpp
	$(MPICXX) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp 
	$(MPICXX) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

clean:
	rm -f tsp-mpi $(OBJ) algorithms.o
