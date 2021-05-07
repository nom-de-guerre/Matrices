DEBUG=-D__DEBUG
OPTIONS= $(OUTSIDE)
CC=g++
CFLAGS=-Wall -I. $(DEBUG) $(OPTIONS)
HDEPS = matrix.h
DEPS = Makefile $(HDEPS)

eigen: eigen.cc francis.cc $(DEPS)
	$(CC) eigen.cc francis.cc -o $@ $(CFLAGS)
	./eigen

example: example.cc $(DEPS)
	$(CC) example.cc -o $@ $(CFLAGS)
	./example

regression: regression.cc $(DEPS)
	$(CC) regression.cc -o $@ $(CFLAGS)
	./regression

all: example regression eigen

clean:
	rm regression
	rm example
	rm eigen
