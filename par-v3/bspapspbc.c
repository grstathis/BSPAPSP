#include "bspedupack.h"




#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#define N 2000

void main(){

 int* gen_graph(int n, double p); 
 int nloc(int p, int s, int n);
int compare_int (const int *a, const int *b);
 void mainloop();


bsp_init(mainloop, 0, NULL );
mainloop();

}

//used for sorting an array
int compare_int (const int *a, const int *b)
{
  return (int) (*a - *b);
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
    /* Compute number of components of processor s for vector
       of length n cyclicly distributed over p processors. */

	return  (n+p-s-1)/p ;
	
} /* end nloc */



void mainloop(){

//int init[N*N] = {0,3,8,1000,-4, 1000,0,1000,1,7,1000,4,0,1000,1000,
//2,1000,-5,0,1000,1000,1000,1000,6,0};
   
int nlr,nlc,s,t,i,j,k,l,li,lsize,tsize0, tsize1,tempp,tempoff,rpos,cpos, 
*lpart,*linter,*gindx,*lcol,*lrow,*lsrow, *lscol, *ltrow, *ltcol, *temp;

int* init = gen_graph(N, 0.05);  

bsp_begin(bsp_nprocs());

/**********Initialization SuperStep 0***************/

//Compute global row and column indeces for each element
int pm = sqrt(bsp_nprocs());
int pn = (bsp_nprocs())/pm;
/* Compute 2D processor numbering from 1D numbering 
 with failsafe if the number of processors are not enough, back to simple 1D cyclic distribution */ 
if ( pn  != pm ){
	pn = bsp_nprocs();
	pm = 1;
	t = bsp_pid();
	s = 0;
  
}else{
	s= bsp_pid()%pm;  /* 0 <= s < pm */
	t= bsp_pid()/pn;  /* 0 <= t < pn */
}

nlr=  nloc(pm,s,N); /* number of local rows */
nlc=  nloc(pn,t,N); /* number of local columns */

lsize = nlr*nlc;						  //interpret 2D size to array size
lpart = vecalloci(lsize);				  //Initialize local part of processor s
linter = vecalloci(lsize);				  //Intermidiate array used for the matrix "multiplication"
gindx = vecalloci(lsize);				  //Array to store the global indeces of the local elements
lcol  = vecalloci(lsize);				  //Array to store the glocal column index
lrow  = vecalloci(lsize);				  //Array to store the glocal row index
bsp_push_reg(lpart,lsize*SZINT);

//Distribute the Data
li=0;
for ( i= 0; i < N; i++){
	for ( j= 0; j < N; j++){
		if ((j % pn) == t){
			lpart[li] = init[N*i+j];
			lrow[li] = i;
			lcol[li] = j;
			gindx[li] = N*i+j;
			li++;	
		}
	}
}


/*for ( i= 0; i < N*N; i++) {

		if(bsp_pid() == (i % bsp_nprocs())){
   			lpart[li] = init[i];
			lrow[li] = i/N;
			lcol[li] = i % N;
			gindx[li] = i;
			li++;	
		}
		

}*/
vecfreei(init);//out of the shared space

tsize0 = tsize1 =lsize;
temp = lrow;

//find unique global rows for processor s
for(i=0;i<tsize0;i++){
    for(j=0;j<tsize0;j++){
         if(i==j){
             continue;
         }
         else if(*(temp+i)==*(temp+j)){
             k=j;
             tsize0--;
             while(k < tsize0){
                 *(temp+k)=*(temp+k+1);
                 k++;
             }
              j=0;
         }
    }
}
temp = lcol;

//find unique global column for processor s
for(i=0;i<tsize1;i++){
    for(j=0;j<tsize1;j++){
         if(i==j){
             continue;
         }
         else if(*(temp+i)==*(temp+j)){
             k=j;
             tsize1--;
             while(k < tsize1){
                 *(temp+k)=*(temp+k+1);
                 k++;
             }
              j=0;
         }
    }
}


//keep unique global rows and columns in arrays
//initialize arrays to hold the elements of those rows and columns(ltcol, ltrow)
lscol  = vecalloci(tsize1); 
lsrow  = vecalloci(tsize0);
ltcol  = vecalloci(N*tsize1);
ltrow  = vecalloci(N*tsize0);

for(i=0;i < tsize0;i++){
    lsrow[i] = lrow[i];
  }
for(i=0;i < tsize1;i++){
    lscol[i] = lcol[i];
  }


vecfreei(lcol);//not needed from this point on
vecfreei(lrow);//we use lscol, lsrow, ltrow, ltcol

//sort arrays
qsort (lsrow, tsize0, sizeof(int), compare_int);
qsort (lscol, tsize1, sizeof(int), compare_int);
bsp_sync();
/**********End Initialization SuperStep 0***************/

double time0= bsp_time();
/*********Repeated Squaring loop start*************/
j=1;
while ((N-1) > j) {

/*************Comm. SuperStep j0*************/
for(i=0;i < tsize1;i++){
	for(k=0; k<N;k++){
		tempp=((N*k+lscol[i]) % bsp_nprocs());
		tempoff = ((double)(N*k+lscol[i])/(double)bsp_nprocs());
		bsp_get(tempp, &lpart[0],tempoff*SZINT, &ltcol[N*i+k],SZINT);
	} 
}

for(i=0;i < tsize0;i++){
	for(k=0; k<N;k++){
		tempp=((N*lsrow[i]+k) % bsp_nprocs());
		tempoff = ((double)(N*lsrow[i]+k)/(double)bsp_nprocs());
		bsp_get(tempp, &lpart[0],tempoff*SZINT, &ltrow[N*i+k],SZINT);
	} 
}
bsp_sync();
/*************End Comm. SuperStep j0*************/

/*************Comp. SuperStep j1*************/
for ( i=0; i<lsize; i++) {
  
	int gcol = gindx[i] % N; //get global col indx of current element
	int grow = gindx[i]/N;	 //get global row indx of current element

    linter[i]=1000;//initiliaze array
	//find appropriate indx of the global rows and columns to perform "multiplication"
	/*for ( l=0; l < tsize0;l++){
		if(grow == lsrow[l]){
			rpos =l;
			break;
		}
	}*/
	int *rp = bsearch (&grow, lsrow, tsize0, sizeof (lsrow),compare_int);
	rpos = rp - lsrow;
	

	int *cp = bsearch (&gcol, lscol, tsize1, sizeof (lscol),compare_int);
	cpos = cp - lscol;
	
	/*for ( l=0; l < tsize1;l++){
		if(gcol == lscol[l]){
			cpos =l;
			break;
		}
	}*/

	//this is where the update is done
	for(k=0;k<N;k++){
		linter[i] = fmin(linter[i], ltrow[N*rpos + k]+ltcol[N*cpos + k]);
	}

}

memcpy(lpart,linter,lsize*SZINT);
j = 2*j;
bsp_sync();
/*************End Comp. SuperStep j1*************/

}
/*********Repeated Squaring loop end*************/
double time1= bsp_time();
bsp_sync();
/*********display matrices and time*********/
if(bsp_pid()==0){
	printf( " \n Block Cyclic Distr  calculation of APSP took: %f seconds \n", time1-time0 ); 
}
/*printf("\n The array is, proc %d \n ", bsp_pid());
  for(i=0;i < lsize;i++){
    	printf(" %d",lpart[i]);
	
}*/
printf("\n ");

//clean up
bsp_pop_reg(lpart);
vecfreei(lpart);
vecfreei(linter);
vecfreei(lscol);
vecfreei(lsrow);
vecfreei(ltcol);
vecfreei(ltrow);
vecfreei(gindx);

bsp_end();   
}


 
