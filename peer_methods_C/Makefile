all: clean main cleanObj

test: clean utest

main: main.o CLab.o peerMethods.o utilities.o
	gcc -Wall -pedantic -std=c11 -mfpmath=387 main.o CLab.o utilities.o peerMethods.o -o main.out -lcblas -lm

main.o: main.c CLab.h peerMethods.h utilities.h
	gcc -c main.c

peerMethods.o: peerMethods.h peerMethods.c
	gcc -c peerMethods.c

CLab.o: CLab.h CLab.c
	gcc -c CLab.c

utilities.o: utilities.h utilities.c
	gcc -c utilities.c

cleanObj:
	rm -f *.o

clean:
	rm -f *.o
	rm -f *.out