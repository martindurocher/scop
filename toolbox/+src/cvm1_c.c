/*=================================================================
 *
 * cvm1_c.C	
 *
 * SUMAMRY: Calculate the cramer von mises criterion for univariate copula 
 *
 * The calling syntax is:
 *
 *		v = cvm1_c(x,y)
 *
 * This is a MEX-file for MATLAB.  
 * Author Martin du Rocher <martin.durocher@uqtr.ca>
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/* Input Arguments */
#define	X_IN	prhs[0]
#define	Y_IN	prhs[1]

/* Output Arguments */
#define	V_OUT	plhs[0]

void cvm1_c(
		   double	x[],
		   double	y[],
           double	v[],
           size_t mx,
           size_t my
        )
{
    double sum1,sum2,sum3;
    int i,j;
    
    /* Initialize sum*/
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
   
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum1 += 1-MAX(x[i],x[j]);

    for(i=0; i<mx ;i++)
        for (j=0; j<my; j++)
           sum2 += 1-MAX(x[i],y[j]);

    for(i=0; i<my ;i++)
        for (j=0; j<my; j++)
          sum3 += 1-MAX(y[i],y[j]);
    
    v[0] = sum1 / mx / mx - sum2 /mx / my * 2 + sum3 / my / my;

}

void mexFunction( 
        int nlhs, 
        mxArray *plhs[], 
		int nrhs, 
        const mxArray*prhs[] )
{ 
    double *x;
    double *y; 
    double *v;
            
    /* Get the dimensions of XY. */
    size_t mx,my;
    
    mx = mxGetM(X_IN);        
    my = mxGetM(Y_IN);
    
    /* Create a matrix for the return argument */ 
    V_OUT = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    v = mxGetPr(V_OUT);
        
    /* Do the actual computations in a subroutine */
    cvm1_c(x,y,v,mx,my); 
    return;
}