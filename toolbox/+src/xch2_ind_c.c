/*=================================================================
 *
 * xch2_c.C	
 *
 * Description: calculate the Cramer von mises criterion with the 
 * product copula.
 *
 * The calling syntax is:
 *
 *		v = xch2_ind_c(x1,x2)
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
#define	X1_IN	prhs[0]
#define	X2_IN	prhs[1]

/* Output Arguments */
#define	V_OUT	plhs[0]

void xch2_ind_c(
		   double	x1[],
           double	x2[],
           double	v[],
           size_t mx)
{
    int i,j;
    
    double sum1, sum2, sum3;
    
    /* Initialize sum*/
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum1 += (1-MAX(x1[i],x1[j])) * (1-MAX(x2[i],x2[j]));
    
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum2 += (1-MAX(x1[i],x2[j])) * (1-MAX(x2[i],x1[j]));
          
    for(i=0; i<mx ;i++)
          sum3 += (1-x1[i]*x1[i])*(1-x2[i]*x2[i]);
    
    
    v[0] = (sum1 + sum2) / (2*mx*mx) - (0.5 * sum3 / mx) + (1.0/9.0) ;

    
}

void mexFunction( 
        int nlhs, 
        mxArray *plhs[], 
		int nrhs, 
        const mxArray*prhs[] )
{ 
    double *x1;
    double *x2; 
    double *v;
            
    /* Get the dimensions of XY. */
    size_t mx;
    
    mx = mxGetM(X1_IN);        
    
    /* Create a matrix for the return argument */ 
    V_OUT = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x1 = mxGetPr(X1_IN);
    x2 = mxGetPr(X2_IN);
    v = mxGetPr(V_OUT);
        
    /* Do the actual computations in a subroutine */
    xch2_ind_c(x1,x2,v,mx); 
    return;
}
