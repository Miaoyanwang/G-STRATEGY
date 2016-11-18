CC=gcc
LIB=-lc -lm
CFLAGS = -O2 -Wno-unused-result
LDLIBS = -lz -lm



OBJECTS =  mastor_tpr.o read.o calc.o mle.o svdcomp.o brent.o brak.o nrutil.o cholesky.o datacheck.o tpr.o


all: mastor_tpr GSTRATEGY

.PHONY: CLAPACK_routine

mastor_tpr: $(OBJECTS) CLAPACK_routine
	$(CC) $(OBJECTS) $(wildcard clapack_routine/*.o)  $(LDLIBS) -o mastor_tpr

CLAPACK_routine:
	$(MAKE) -C clapack_routine/


GSTRATEGY: G-STRATEGY.o
	$(CC) $(wildcard clapack_routine/*.o) $(LIB) G-STRATEGY.o -o G-STRATEGY

clean:
	@rm -rf *.o clapack_routine/*.o *MASTOR.txt TPR.txt

