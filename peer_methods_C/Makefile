all: clean main

main: main.o CLab.o utilities.o
	gcc main.o CLab.o utilities.o -o main.out

main.o: main.c CLab.h utilities.h
	gcc -c main.c

CLab.o: CLab.h CLab.c
	gcc -c CLab.c

utilities.o: utilities.h utilities.c
	gcc -c utilities.c

clean:
	rm -f *.o
	rm -f *.out