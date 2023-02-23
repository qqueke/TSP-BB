CC = gcc
CFLAGS = -O3
SRC = tsp.c queue.c
OBJ = $(SRC:.c=.o)
INCLUDES = -I/home/qqueke/PDC/PDC---Projeto/queue
LIBS =

tsp: $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f tsp $(OBJ)
