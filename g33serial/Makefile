CC = g++
CFLAGS = -O3
SRC = tsp.cpp
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o # Add algorithms.o as a dependency

tsp: tsp.o algorithms.o # Add algorithms.o as a dependency
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp.o

tsp.o: tsp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp # Add a rule to compile algorithms.cpp to algorithms.o
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f tsp $(OBJ) algorithms.o