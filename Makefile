CC = g++
CFLAGS = -O3
SRC = tsp.cpp tsp-omp.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o # Add algorithms.o as a dependency

all: tsp tsp-omp

tsp: tsp.o algorithms.o # Add algorithms.o as a dependency
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp.o

tsp-omp: tsp-omp.o algorithms.o # Add algorithms.o as a dependency
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-omp.o

tsp.o: tsp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

tsp-omp.o: tsp-omp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp # Add a rule to compile algorithms.cpp to algorithms.o
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

.PHONY: all

clean:
	rm -f tsp tsp-omp $(OBJ) algorithms.o