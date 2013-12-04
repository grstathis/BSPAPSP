#include "bspedupack.h"




#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>

#define N 500
#define IMAX 1000
void main(){

void mainloop();



//int init[N*N] = {0,3,8,1000,-4, 1000,0,1000,1,7,1000,4,0,1000,1000,
//2,1000,-5,0,1000,1000,1000,1000,6,0};
bsp_init(mainloop,0,NULL);
mainloop();

}

void mainloop(){

int* gen_graph(int n, double p); 
void matrixMultiply(int *res, int *input );
int* init = gen_graph(N, 0.05);  



bsp_begin(1);


  int inter[N*N];	//intermidiate matrix used from the repeated squaring algorithm 
  int path[N*N];	//final APSP matrix
  int i,j,k,l;		//indeces for general usage on matrices	   

/*****Initialization of final APSP matrix****/

  for ( i=0; i<N*N; i++) {
                path[i] = init[i];
  }



/*****Initialization end********************/

double time0= bsp_time();

/*********Repeated Squaring loop start*************/
  i=1;
  while ((N-1) > i) {

    matrixMultiply(inter, path); 	//<----- the updating of values of the APSP is done here
    memcpy(path,inter,N*N*SZINT); 	//<----  keep results and continue loop
   	i=2*i;							//<----  logN repeats
  }
/*********Repeated Squaring loop end*************/
double time1= bsp_time();

 /*********display matrices and time*********/
printf( " \n Sequential calculation of APSP took: %f seconds \n", time1-time0 ); 

/*printf("Init matrix \n");

 for(i=0;i<N;i++)
 {
     for(j=0;j<N;j++){
      printf("\t %d",init[N*i+j]);
    }
    printf("\n");

 }

 printf("\n APSP matrix \n");

 for(i=0;i<N;i++)
 {
    
   for(j=0;j<N;j++){
    printf("\t %d",path[N*i+j]);
   }
    printf("\n");

 }*/
free(init);
bsp_end();

}

//Description: the matrix "funny" doubling multiplication function
//Input: pointer to N*N array A
//Output: pointer to N*N array B = A^2
 void inline matrixMultiply(int *res, int* restrict input) {
    int i,j,k;

     for ( i=0; i<N; i++) {
        for ( j=0; j<N; j++) {
			res[N*i+j] = IMAX;
			for ( k=0; k<N; k++) {
                        res[N*i+j] = fmin(res[N*i+j], input[N*i+k]+input[N*k+j]); 
            }
        }
    }

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




