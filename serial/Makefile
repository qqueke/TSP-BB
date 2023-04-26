CC = g++
CFLAGS = -O3
SRC = tsp.cpp 
OBJ = $(SRC:.cpp=.o)
LIBS = algorithms.o 

tsp: tsp.o algorithms.o 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -fopenmp -o $@ tsp.o

tsp.o: tsp.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

algorithms.o: algorithms.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -fopenmp -c -o $@ $<

clean:
	rm -f tsp $(OBJ) algorithms.o
