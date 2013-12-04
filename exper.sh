#!/bin/bash
echo "" > /home/stathis/Dropbox/paralg/APSPc/apsp.out
for j in 5 32 100 200 400 500 1000 2000
do
	echo "--------------------" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	echo " size of matrix  $j" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	sed -i -e "s/.*#define N.*/#define N $j/g" /home/stathis/Dropbox/paralg/APSPc/\seq/bspapspseq.c
	sed -i -e "s/.*#define N.*/#define N $j/g" /home/stathis/Dropbox/paralg/APSPc/par-v1/bspalpspbr.c
	sed -i -e "s/.*#define N.*/#define N $j/g" /home/stathis/Dropbox/paralg/APSPc/par-v2/bspalpspbnt.c
	sed -i -e "s/.*#define N.*/#define N $j/g" /home/stathis/Dropbox/paralg/APSPc/par-v3/bspapspbc.c
	
	cd /home/stathis/Dropbox/paralg/APSPc/seq
	make
	cd /home/stathis/Dropbox/paralg/APSPc/par-v1
	make
	cd /home/stathis/Dropbox/paralg/APSPc/par-v2
	make
	echo "-------------------" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	cd /home/stathis/Dropbox/paralg/APSPc/
	./seq/allpairseq  >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
for i in 1 2 4 8 16 20 32
do
	export MP_PROCS=$i
	echo "" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	echo " num of procs $i" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	echo "************" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	#echo "version 1 copy all matrix Row Distr" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	cd /home/stathis/Dropbox/paralg/APSPc/par-v1
	mpirun -np $i allpairrb	>> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	#echo "version 2 Row Distr need to know" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	cd /home/stathis/Dropbox/paralg/APSPc/par-v2
	mpirun -np $i allpairrb_ntn >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
	echo "" >> /home/stathis/Dropbox/paralg/APSPc/apsp.out
done
done
