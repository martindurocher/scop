/*=================================================================
 *
 * pairwise.C	
 *
 * The function identify all pairs of a vector
 *
 * The calling syntax is:
 *
 *		y = pairwise_c(x)
 *
 *  You may also want to look at the corresponding M-code, pairwise.m.
 *
 * This is a MEX-file for MATLAB.  
 * Author Martin du Rocher <martin.durocher@uqtr.ca>
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"

/* Input Arguments */
#define	X_IN	prhs[0]

/* Output Arguments */
#define	Y_OUT	plhs[0]

void pairwise_c(
		   double	x[],
		   double	y[],
           size_t m,
           size_t p
		   )
{
   int ii,jj,kk;
   kk=0;
   for(ii = 0; ii < m-1; ii++){
     for(jj = ii + 1 ; jj < m  ;jj++){
       y[kk] = x[ii]; 
       y[kk+p] = x[jj];
       kk++;
     }
   }
   
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *y; 
    double *x; 
    size_t m,p;
        
    /* Get the dimensions of XY. */ 
    m = mxGetM(X_IN);        
    p = (m-1)*m/2;
    
    /* Create a matrix for the return argument */ 
    Y_OUT = mxCreateDoubleMatrix( (mwSize)p, (mwSize)2, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN);
    y = mxGetPr(Y_OUT);
        
    /* Do the actual computations in a subroutine */
    pairwise_c(x,y,m,p); 
    return;
}