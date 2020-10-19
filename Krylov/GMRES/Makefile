DEBUG=-D__DEBUG
OPTIONS= $(OUTSIDE)
CC=g++
CFLAGS=-Wall -I. -I.. -I../.. $(DEBUG) $(OPTIONS)
HDEPS = ../../matrix.h GMRES.h ../Krylov.h ../SparseMatrix.h
DEPS = Makefile $(HDEPS)

all: GMRES_example

GMRES_example: GMRES_example.cc $(DEPS)
	$(CC) GMRES_example.cc -o $@ $(CFLAGS)

clean:
	rm GMRES_example
