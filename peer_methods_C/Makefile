CC = gcc
INCLUDE_DIRS = -I peerMethods/include
LINKED_LIB = -L peerMethods/lib
CFLAGS = -Wall -pedantic -std=c11 -mfpmath=387 -Wextra
LIB_FLAGS = -lpeerMethods -lm -lcblas
OBJECTS = main.o helper.o

.PHONY: all

all: cleaner main clean

main: ${OBJECTS}
	${CC} ${INCLUDE_DIRS} ${CFLAGS} ${OBJECTS} -o $@.out ${LINKED_LIB} ${LIB_FLAGS}

main.o: main.c helper.h
	${CC} ${INCLUDE_DIRS} ${CFLAGS} -c $< ${LIB_FLAGS}

helper.o: helper.c helper.h
	${CC} ${CFLAGS} -c $<

.PHONY: clean

clean:
	rm -f *.o

cleaner: clean
	rm -f *.out