CC=icpc

CFLAGS= -Ofast -Wall -lgsl

all: project

project: main.o  read.o
	$(CC) $(CFLAGS)  main.o read.o
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
read.o: read.cpp
	$(CC) $(CFLAGS) -c read.cpp

clean:
	rm *.o a.out
