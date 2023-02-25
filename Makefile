CC = g++
CFLAGS = -O3
SRC = tsp.cpp
OBJ = $(SRC:.cpp=.o)
LIBS =

tsp: $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f tsp $(OBJ)
