CC= mpicc
CFLAGS= -O3 -g
LFLAGS= -lm -lbsponmpi 
OBJALP= bspalpspbnt.o bspedupack.o


all:allpairrb_ntn



allpairrb_ntn: $(OBJALP)
	$(CC) $(CFLAGS) -o allpairrb_ntn $(OBJALP) $(LFLAGS)	


	
bspalpspbnt.o: bspalpspbnt.c bspedupack.h
	$(CC) $(CFLAGS) -c bspalpspbnt.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o  allpairrb_ntn
