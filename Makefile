CC      ?= cc
CFLAGS  ?= -std=c99 -Wall -Wextra -pedantic -O2
LDLIBS   = -lm

all: libsimplex.a example

libsimplex.a: src/simplex.o
	$(AR) rcs $@ $^

example: example.o libsimplex.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

tests/test_simplex: tests/test_simplex.o libsimplex.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

test: tests/test_simplex
	./tests/test_simplex

clean:
	rm -f src/*.o *.o tests/*.o libsimplex.a example tests/test_simplex

.PHONY: all test clean
