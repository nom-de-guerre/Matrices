# DEBUG=-D__DEBUG
DEBUG=-O3
OPTIONS= $(OUTSIDE)
CC=g++
CFLAGS=-Wall -I. -I.. -I../.. $(DEBUG) $(OPTIONS)
HDEPS = ../../matrix.h ConjugateGradient.h
DEPS = Makefile $(HDEPS)

all: CG_Example

CG_Example: CG_Example.cc $(DEPS)
	$(CC) CG_Example.cc -o $@ $(CFLAGS)

clean:
	rm CG_Example
