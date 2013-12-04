CC= mpicc
CFLAGS= -O3 -g
LFLAGS= -lm -lbsponmpi 
OBJALP= bspalpspbr.o bspedupack.o


all:allpairrb



allpairrb: $(OBJALP)
	$(CC) $(CFLAGS) -o allpairrb $(OBJALP) $(LFLAGS)	


	
bspalpspbr.o: bspalpspbr.c bspedupack.h
	$(CC) $(CFLAGS) -c bspalpspbr.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o  allpairrb
