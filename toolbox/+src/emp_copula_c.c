/*=================================================================
 *
 * emp_copula_c.C	
 *
 * description: Calculate the emprirical copula
 *
 * The calling syntax is:
 *
 *		y = emp_copula_c(x)
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
#define	B1_IN	prhs[1]
#define	B2_IN	prhs[2]

/* Output Arguments */
#define	C_OUT	plhs[0]

void emp_copula_c(
		   double	x[],
           double	b1[],
           double	b2[],
           double	c[],
           size_t m1,
           size_t m2,
           size_t n
        )
{
    int i,j,k;
    double sum;
    double w = n;
    
    w = 1.0 / w;
       
    for(i = 0; i < m1; i++)
    {
        for(j = 0; j < m2; j++)
        {
            sum = 0;
            for(k = 0; k < n; k++)
            {
                if((x[k] <= b1[i]) & (x[k+n] <= b2[j]))
                    sum += w;
                
                if(x[k] > b1[i])
                    break;
            }
            c[m1*j+i] = sum;
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
    double *b1;
    double *b2;
    double *c;
            
    /* Get the dimensions of XY. */
    size_t m1,m2,n;
    
    n  = mxGetM(X_IN);
    m1 = mxGetN(B1_IN); 
    m2 = mxGetN(B2_IN);
    
    /* Create a matrix for the return argument */ 
    C_OUT = mxCreateDoubleMatrix( (mwSize)m1, (mwSize)m2, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x = mxGetPr(X_IN);
    b1 = mxGetPr(B1_IN);
    b2 = mxGetPr(B2_IN);
    c = mxGetPr(C_OUT);
        
    /* Do the actual computations in a subroutine */
    emp_copula_c(x,b1,b2,c,m1,m2,n); 
    return;
}