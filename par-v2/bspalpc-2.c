#include "bspedupack.h"




#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#define N 5


void main(){
 
 int* gen_graph(int n, double p); 
 int nloc(int p, int s, int n);
 void mainloop();


bsp_init(mainloop, 0, NULL );
mainloop();

}


//Description: generate random graph with integer weights of size n*n, with probability p
int* gen_graph(int n, double p)
{
	int i,j;
	int* l = calloc(n*n, sizeof(int));
	//struct mt19937p state;
	//sgenrand(10302011UL, &state);
	for ( j = 0; j < n; ++j) {
		for (i = 0; i < n; ++i)
		  if(((double) rand() / (RAND_MAX)) > p){
			l[j*n+i] = rand() % 1000;
		  }
			l[j*n+j] = 0;
		}	
	return l;
}



int nloc(int p, int s, int n){
    /* Compute number of rows of processor s for vector
       of length n block distributed over p processors. */
 	double t = ((double)n/(double)p);
	if (s==(p-1)){
    	return  (n - s*(ceil(t))) ;
	}else{
		return  (ceil(t));
	} 
} /* end nloc */



void mainloop(){

int init[N*N] = {0,3,8,1000,-4, 1000,0,1000,1,7,1000,4,0,1000,1000,
2,1000,-5,0,1000,1000,1000,1000,6,0};

int i,j,k,l,v,t,lsize,*lsize_m,llsize,*lrow,*lcol, *linit, *linter, *lresp,*startrow_m;
int li,lj,lk,startrow, endrow,g;

int* init = gen_graph(N, 0.05);  

bsp_begin(bsp_nprocs());




lsize = nloc(bsp_nprocs(),bsp_pid(), N);
lrow = vecalloci(lsize*N);
lcol = vecalloci(N);
startrow_m = vecalloci(bsp_nprocs());
lsize_m = vecalloci(bsp_nprocs());
linter = vecalloci(lsize*N);

bsp_push_reg(startrow_m,bsp_nprocs()*SZINT);
bsp_push_reg(lsize_m,bsp_nprocs()*SZINT);
bsp_push_reg(lrow,lsize*N*SZINT);


if(bsp_pid() == (bsp_nprocs() - 1)){
 startrow = (N - lsize);
 endrow = N;
}else{
 startrow = bsp_pid()*lsize;
 endrow = bsp_pid()*lsize + lsize;
}


//bsp_sync();


li=0;
for ( i= startrow; i < endrow; i++) {
	lj=0;
	 for(j=0; j < N; j++) {	
   		lrow[N*li+lj] = init[N*i+j];
		lj++;
   	 } 
 li++;
}

for ( i=0; i<bsp_nprocs(); i++) {
			startrow_m[i] = 0;
			lsize_m[i] = 0;
}
bsp_sync();
for(g=0; g<bsp_nprocs();g++){
bsp_put(g,&startrow,&startrow_m[0],bsp_pid()*SZINT,SZINT);
bsp_put(g,&lsize,&lsize_m[0],bsp_pid()*SZINT,SZINT);
}
bsp_sync();

j=1;
while ((N-1) > j) {
 
		bsp_sync();		

		for ( i=0; i<N*lsize; i++) {
			linter[i] = 1000;
		}
		for ( i=0; i<N; i++) {
			lcol[i] = 0;
		}

	   		
        	for ( lj=0; lj < N; lj++) {
				bsp_sync();
				t=0;
				for(g=0; g < bsp_nprocs();g++){
				  for(v=0; v<lsize_m[g]; v++){				
					bsp_get(g,&lrow[0],(lj+v*N)*SZINT,&lcol[t],SZINT);
					t++;
				  }
				}
				bsp_sync();

				for ( li = 0; li < lsize; li++){
				//printf("\n I am proc %d and i have row %d and i have lcol %d \n",bsp_pid(),li,lj);
				//printf("\n i have row %d and i have lcol %d \n",li,lj);
				//bsp_sync();
				for ( lk=0; lk < N; lk++) {
					//if(bsp_pid()==0){
					//	printf("row \t %d col \t %d \n",lrow[N*li+lk], lcol[lk]);
					//}
					linter[N*li+lj] = fmin(linter[N*li+lj], lrow[N*li+lk]+lcol[lk]);
				} 
				

        	}
			bsp_sync();
    	}
 



		memcpy(lrow,linter,N*lsize*SZINT);
  		j=2*j;
		
		bsp_sync();
}

bsp_sync();

for(g = 0; g < bsp_nprocs(); g++){
if(bsp_pid()==g){
 printf("\n i am proc %d and i have APSP Mat \n",bsp_pid());
  for(k=0;k<lsize;k++)
     {
	  printf("\n");
		 for(l=0;l<N;l++){
		    printf("\t %d",lrow[N*k+l]);
			  }
			printf("\n \n ");
		}
	}
	bsp_sync();
}


bsp_end();   
}


 
