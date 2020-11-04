DJS_MATH=/Users/dsantry/Desktop/Shuttle/Matrix\ V2/
OPTIONS= $(OUTSIDE)
CC=g++
# DEBUG=-O0 -g
DEBUG=-O3 -mavx
CFLAGS=-Wall -I. $(DEBUG) $(OPTIONS) -I$(DJS_MATH) $(DEBUG)
HDEPS = NNLM.h NNLM.tcc data.h
DEPS = *.h $(DJS_MATH)/matrix.h Makefile $(HDEPS)

classify: classify.cc $(DEPS) softmax.h
	$(CC) -o $@ $(CFLAGS) classify.cc

sine: sine.cc $(DEPS) regression.h
	$(CC) -o $@ $(CFLAGS) sine.cc

all: sine classify

clean:
	rm classify
	rm sine
