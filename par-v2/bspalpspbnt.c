#include "bspedupack.h"




#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#define N 10000


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
	int t = ceil(((double)n/(double)p));
	if((s+1)*t > n ){
		if(s*t < n){
			return (n - s*t);
		}else{
			return 0;
		}
	}else{
		return  (ceil(t));
	}
} /* end nloc */



void mainloop(){

//int init[N*N] = {0,3,8,1000,-4, 1000,0,1000,1,7,1000,4,0,1000,1000,
//2,1000,-5,0,1000,1000,1000,1000,6,0};

int i,j,k,l,v,t,lsize,*lsize_m,*lrow,*lcol, *linit, *linter,*startrow_m;
int li,lj,lk,startrow, endrow,g;

int* init = gen_graph(N, 0.05);  

bsp_begin(bsp_nprocs());


/**********Initialization***************/

/*******Comp. Superstep 0******/

lsize = nloc(bsp_nprocs(),bsp_pid(), N); //Get the number of rows of processor s
lrow = vecalloci(lsize*N);				 //The main storing array of processor s
lcol = vecalloci(N);					 //array to hold the column for the matrix squaring
startrow_m = vecalloci(bsp_nprocs());    //array to hold all processors starting global row
lsize_m = vecalloci(bsp_nprocs());		 //array to hold the number of rows of all processors
linter = vecalloci(lsize*N);			 //Intermidiate array used for the matrix "multiplication"

bsp_push_reg(startrow_m,bsp_nprocs()*SZINT);
bsp_push_reg(lsize_m,bsp_nprocs()*SZINT);
bsp_push_reg(lrow,lsize*N*SZINT);

/****Get the first and last global row of processor s***/
if(bsp_pid() == (bsp_nprocs() - 1)){
 startrow = (N - lsize);
 endrow = N;
}else{
 startrow = bsp_pid()*lsize;
 endrow = bsp_pid()*lsize + lsize;
}



//Distribute Data, according row block distribution
li=0;
for ( i= startrow; i < endrow; i++) {
	lj=0;
	 for(j=0; j < N; j++) {	
   		lrow[N*li+lj] = init[N*i+j];
		lj++;
   	 } 
 li++;
}
vecfreei(init); //out of the shared enviroment

//initialize arrays
for ( i=0; i<bsp_nprocs(); i++) {
			startrow_m[i] = 0;
			lsize_m[i] = 0;
}

bsp_sync();
/*******End Comp. Superstep 0******/


/*********Comm. Superstep 1********/
//Communicate the global starting rows of all processors
for(g=0; g<bsp_nprocs();g++){
	bsp_put(g,&startrow,&startrow_m[0],bsp_pid()*SZINT,SZINT);
	bsp_put(g,&lsize,&lsize_m[0],bsp_pid()*SZINT,SZINT);
}
/*********End Comm. Superstep 1*****/
bsp_sync();
/**********End Initialization***************/

double time0= bsp_time();
/*********Repeated Squaring loop start*************/
j=1;
while ((N-1) > j) {
 
		/****Comp. Superstep j0****/ 
		//initialize arrays
		for ( i=0; i<N*lsize; i++) {
			linter[i] = 1000;
		}
		for ( i=0; i<N; i++) {
			lcol[i] = 0;
		}
		bsp_sync();
		/****End Comp. Superstep j0****/ 
	   		
        	for ( lj=0; lj < N; lj++) {
				/***Comm. SuperStep jlj0*******/
				//get global column lj 
				t=0;
				for(g=0; g < bsp_nprocs();g++){
				  for(v=0; v<lsize_m[g]; v++){				
					bsp_get(g,&lrow[0],(lj+v*N)*SZINT,&lcol[t],SZINT);
					t++;
				  }
				}
				bsp_sync();
				/***End Comm. SuperStep jlj0***/
				/***Comp. SuperStep jlj1*******/
				//update the values that use global column lj
				for ( li = 0; li < lsize; li++){
					for ( lk=0; lk < N; lk++) {
						linter[N*li+lj] = fmin(linter[N*li+lj], lrow[N*li+lk]+lcol[lk]);
					} 
        		}
				bsp_sync();
				/***End Comp. SuperStep jlj1***/
    		}
 		/****Comp. Superstep j1****/ 
		memcpy(lrow,linter,N*lsize*SZINT);
  		j=2*j;
		bsp_sync();
		/****End Comp. Superstep j1****/ 
}
/*********Repeated Squaring loop end*************/
double time1= bsp_time();
bsp_sync();
/*********display matrices and time*********/
if(bsp_pid()==0){
	printf( " \n Block Row Distr (need to know basis) calculation of APSP took: %f seconds \n", time1-time0 ); 
}

/*for(g = 0; g < bsp_nprocs(); g++){
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
}*/


//Clean up
bsp_pop_reg(startrow_m);
bsp_pop_reg(lsize_m);
bsp_pop_reg(lrow);


vecfreei(lrow);
vecfreei(lcol);
vecfreei(startrow_m);
vecfreei(lsize_m);
vecfreei(linter);

bsp_end();   
}


 
