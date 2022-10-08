all: clean main

test: clean utest

utest: utest.o CLab.o utilities.o
	gcc utest.o CLab.o utilities.o -o utest.out -lcblas

utest.o: utest.c CLab.h utilities.h
	gcc -c utest.c

main: main.o CLab.o peerMethods.o utilities.o
	gcc -Wall -ansi -pedantic main.o CLab.o utilities.o peerMethods.o -o main.out -lcblas -lm

main.o: main.c CLab.h peerMethods.h utilities.h
	gcc -c main.c

peerMethods.o: peerMethods.h peerMethods.c
	gcc -c peerMethods.c

CLab.o: CLab.h CLab.c
	gcc -c CLab.c

utilities.o: utilities.h utilities.c
	gcc -c utilities.c

clean:
	rm -f *.o
	rm -f *.out