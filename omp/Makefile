CC = g++
CFLAGS = -O3
SRC = tsp-omp.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o 

all: tsp-omp 

tsp-omp: tsp-omp.o algorithms.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp-omp.o

tsp-omp.o: tsp-omp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

clean:
	rm -f tsp-omp $(OBJ) algorithms.o
