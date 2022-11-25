CC = gcc
INCLUDE_DIRS = -I ../external_libs/CBLAS/include
LINKED_LIB = -L ../external_libs/CBLAS/lib
FLAGS = -Wall -pedantic -std=c11 -mfpmath=387 -Wextra -fPIC
ARCHFLAGS = -Wl,--whole-archive -shared
LIB_FLAGS = -lcblas -lm -lblas
OBJECTS = CLab.o utilities.o peerMethods.o
CLIB = ../lib/libpeerMethods.so

.PHONY: all

all: cleaner lib clean

lib: ${OBJECTS}
	${CC} ${FLAGS} ${ARCHFLAGS} -o ${CLIB} ${OBJECTS} -Wl,--no-whole-archive

peerMethods.o: peerMethods.c
	${CC} -I../include ${INCLUDE_DIRS} ${FLAGS} -c peerMethods.c ${LINKED_LIB} ${LIB_FLAGS}

CLab.o: CLab.c
	${CC} -I../include ${FLAGS} -c CLab.c

utilities.o: utilities.c
	${CC} -I../include ${FLAGS} -c utilities.c

.PHONY: clean

clean:
	rm -f *.o

rmlib:
	rm -f ${CLIB}

cleaner: clean rmlib
	rm -f *.out