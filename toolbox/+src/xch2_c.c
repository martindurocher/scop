/*=================================================================
 *
 * xch2_c.C	
 *
 * Description: calculate the Cramer von mises criterion between two
 *  sample using exchangeability.
 *
 * The calling syntax is:
 *
 *		v = xch2_c(x1,x2,y1,y2)
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
#define	Y1_IN	prhs[2]
#define	Y2_IN	prhs[3]

/* Output Arguments */
#define	V_OUT	plhs[0]

void xch2_c(
		   double	x1[],
           double	x2[],
		   double	y1[],
           double	y2[],
           double	v[],
           size_t mx,
           size_t my)
{
    int i,j;
    
    double sum1_1, sum1_2, sum1_3, sum1_4, sum2_2;
    double sum2_3, sum2_4, sum3_3, sum3_4, sum4_4;
    
    /* Initialize sum*/
    sum1_1 = 0;
    sum1_2 = 0;
    sum1_3 = 0;
    sum1_4 = 0;
    sum2_2 = 0;
    sum2_3 = 0;
    sum2_4 = 0;
    sum3_3 = 0;
    sum3_4 = 0;
    sum4_4 = 0;
    
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum1_1 += (1-MAX(x1[i],x1[j])) * (1-MAX(x2[i],x2[j]));
    
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum1_2 += (1-MAX(x1[i],x2[j])) * (1-MAX(x2[i],x1[j]));
          
    for(i=0; i<mx ;i++)
        for (j=0; j<my; j++)
          sum1_3 += (1-MAX(x1[i],y1[j])) * (1-MAX(x2[i],y2[j]));
    
    for(i=0; i<mx ;i++)
        for (j=0; j<my; j++)
          sum1_4 += (1-MAX(x1[i],y2[j])) * (1-MAX(x2[i],y1[j]));
          
    for(i=0; i<mx ;i++)
        for (j=0; j<mx; j++)
          sum2_2 += (1-MAX(x2[i],x2[j])) * (1-MAX(x1[i],x1[j]));
    
    for(i=0; i<mx ;i++)
        for (j=0; j<my; j++)
          sum2_3 += (1-MAX(x2[i],y1[j])) * (1-MAX(x1[i],y2[j]));
  
    for(i=0; i<mx ;i++)
        for (j=0; j<my; j++)
          sum2_4 += (1-MAX(x2[i],y2[j])) * (1-MAX(x1[i],y1[j]));
    
    for(i=0; i<my ;i++)
        for (j=0; j<my; j++)
          sum3_3 += (1-MAX(y1[i],y1[j])) * (1-MAX(y2[i],y2[j]));
        
    for(i=0; i<my ;i++)
        for (j=0; j<my; j++)
          sum3_4 += (1-MAX(y1[i],y2[j])) * (1-MAX(y2[i],y1[j]));
       
    for(i=0; i<my ;i++)
        for (j=0; j<my; j++)
          sum4_4 += (1-MAX(y2[i],y2[j])) * (1-MAX(y1[i],y1[j]));
    
    
    v[0] = 0.25 * ((sum1_1 + sum2_2 + 2 * sum1_2) / (mx * mx) -
           2 * (sum1_3 + sum1_4 + sum2_3 + sum2_4) / ( mx * my) +
           (sum3_3 + sum4_4 + 2 * sum3_4) / (my * my));

    
}

void mexFunction( 
        int nlhs, 
        mxArray *plhs[], 
		int nrhs, 
        const mxArray*prhs[] )
{ 
    double *x1;
    double *x2;
    double *y1;
    double *y2; 
    double *v;
            
    /* Get the dimensions of XY. */
    size_t mx,my;
    
    mx = mxGetM(X1_IN);        
    my = mxGetM(Y1_IN);
    
    /* Create a matrix for the return argument */ 
    V_OUT = mxCreateDoubleMatrix( (mwSize)1, (mwSize)1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    x1 = mxGetPr(X1_IN);
    x2 = mxGetPr(X2_IN);
    y1 = mxGetPr(Y1_IN);
    y2 = mxGetPr(Y2_IN);
    v = mxGetPr(V_OUT);
        
    /* Do the actual computations in a subroutine */
    xch2_c(x1,x2,y1,y2,v,mx,my); 
    return;
}
