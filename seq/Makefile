CC= mpicc
CFLAGS= -O3 -g
LFLAGS= -lm -lbsponmpi 
OBJALP= bspapspseq.o bspedupack.o


all:allpairseq



allpairseq: $(OBJALP)
	$(CC) $(CFLAGS) -o allpairseq $(OBJALP) $(LFLAGS)	


	
bspapspseq.o: bspapspseq.c bspedupack.h
	$(CC) $(CFLAGS) -c bspapspseq.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o  allpairseq
