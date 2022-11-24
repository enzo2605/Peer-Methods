CC = gcc
INCLUDE_DIRS = -I peerMethods/include -I peerMethods/external_libs/CBLAS/include
LINKED_LIB = -L peerMethods/lib -L pedantic/external_libs/CBLAS/lib
FLAGS = -Wall -pedantic -std=c11 -mfpmath=387 -Wextra
LIB_FLAGS = -lcblas -lm -lblas -lpeerMethods
OBJECTS = main.o helper.o

.PHONY: all

all: cleaner main clean

main: ${OBJECTS}
	${CC} ${INCLUDE_DIRS} ${FLAGS} ${OBJECTS} -o main.out ${LINKED_LIB} ${LIB_FLAGS}

main.o: main.c helper.h
	${CC} ${INCLUDE_DIRS} ${FLAGS} -c main.c ${LIB_FLAGS}

helper.o: helper.h helper.c
	${CC} ${FLAGS} -c helper.c

.PHONY: clean

clean:
	rm -f *.o

cleaner: clean
	rm -f *.out