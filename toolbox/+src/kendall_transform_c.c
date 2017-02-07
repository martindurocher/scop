/*=================================================================
 *
 * kendall_transform_c.C	
 *
 * summary: calculate the probability integral transformation 
 *    of a sample
 *
 * The calling syntax is:
 *
 *		y = kendall_transform_c(x)
 *
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
#define	V_OUT	plhs[0]

void kendall_transform_c(
		   double	x[],
           double	v[],
           size_t mx
        )
{
    int i,j;
    double w = mx;
    
    w = 1 / w;
       
    for(i = 0; i < mx; i++)
    {
        v[i] = 0;
        for(j = 0; j < mx; j++)
        {
            if((x[j] <= x[i]) & (x[j + mx] <= x[i + mx]))
            {
               v[i] += w;
            }
        }
    }
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
    V_OUT = mxCreateDoubleMatrix( (mwSize)mx, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN);
    v = mxGetPr(V_OUT);
        
    /* Do the actual computations in a subroutine */
    kendall_transform_c(x,v,mx); 
    return;
}