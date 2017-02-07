/*=================================================================
 *
 * cvm1_c.C	
 *
 * SUMAMRY: Calculate the cramer von mises criterion for univariate 
 * copula against the product copula 
 *
 * The calling syntax is:
 *
 *		v = cvm1_ind_c(x)
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

/* Output Arguments */
#define	V_OUT	plhs[0]

void cvm1_ind2_c(
		   double	x[],
           double	v[],
           size_t mx)
{
    double sum1,sum2;
    int i,j;
    
    /* Initialize sum*/
    sum1 = 0;
    sum2 = 0;
   
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum1 += (1-MAX(x[i],x[j]));

    for(i=0; i<mx ;i++)  
        sum2 += x[i]* x[i] * (2*log(x[i])-3);;
    
    v[0] =  sum1 / mx / mx - sum2 / (2*mx) - 47.0/54.0;

}

void mexFunction( 
        int nlhs, 
        mxArray *plhs[], 
		int nrhs, 
        const mxArray*prhs[] )
{ 
    double *x;
    double *v;
            
    /* Get the dimensions of XY. */
    size_t mx;
    
    mx = mxGetM(X_IN);        
    
    /* Create a matrix for the return argument */ 
    V_OUT = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN);
    v = mxGetPr(V_OUT);
        
    /* Do the actual computations in a subroutine */
    cvm1_ind2_c(x,v,mx); 
    return;
}