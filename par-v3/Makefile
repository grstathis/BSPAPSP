CC= mpicc
CFLAGS= -O3 -g
LFLAGS= -lm -lbsponmpi 
OBJALP= bspapspbc.o bspedupack.o


all:allpairbc

allpairbc: $(OBJALP)
	$(CC) $(CFLAGS) -o allpairbc $(OBJALP) $(LFLAGS)	


	
bspapspbc.o: bspapspbc.c bspedupack.h
	$(CC) $(CFLAGS) -c bspapspbc.c

bspedupack.o: bspedupack.c bspedupack.h
	$(CC) $(CFLAGS) -c bspedupack.c

clean:
	rm -f *.o  allpairbc
