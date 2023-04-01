CC = mpicxx
CFLAGS = -O3
SRC = tsp-mpi.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o 
INCLUDES = -I/usr/include/x86_64-linux-gnu/mpich/

tsp-mpi: tsp-mpi.o algorithms.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-mpi.o -lmpi

tsp-mpi.o: tsp-mpi.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

clean:
	rm -f tsp-mpi $(OBJ) algorithms.o
