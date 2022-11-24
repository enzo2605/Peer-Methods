CC = gcc
INCLUDE_DIRS = -I external_libs/CBLAS/include
LINKED_LIB = -L external_libs/CBLAS/lib
FLAGS = -Wall -pedantic -std=c11 -mfpmath=387 -Wextra
LIB_FLAGS = -lcblas -lm -lblas
OBJECTS = main.o CLab.o utilities.o peerMethods.o

.PHONY: all

all: cleaner main clean

main: ${OBJECTS}
	${CC} ${INCLUDE_DIRS} ${FLAGS} ${OBJECTS} -o main.out ${LINKED_LIB} ${LIB_FLAGS}

main.o: main.c CLab.h peerMethods.h utilities.h
	${CC} ${FLAGS} -c main.c ${LIB_FLAGS}

peerMethods.o: peerMethods.h peerMethods.c
	${CC} ${FLAGS} -c peerMethods.c ${LIB_FLAGS}

CLab.o: CLab.h CLab.c
	${CC} ${FLAGS} -c CLab.c

utilities.o: utilities.h utilities.c
	${CC} ${FLAGS} -c utilities.c

.PHONY: clean

clean:
	rm -f *.o

cleaner: clean
	rm -f *.out