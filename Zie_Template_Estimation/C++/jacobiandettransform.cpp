/*
 *  jacobiandettransform.cpp
 *  
 *
 *  Created by Akshay on 4/19/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>




//using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	

    
    
    if(nrhs==0)
		mexPrintf("takes 3 arguments ");
    if(nrhs<3)
        mexErrMsgTxt("Number of arguments must be 3");
    // mexPrintf("jac det takes 3 arguments ");
   double   *d1, *d2, *d3;
    
    d1 = static_cast<double*>(mxGetData(prhs[0]));
    d2 = static_cast<double*>(mxGetData(prhs[1]));
	d3 = static_cast<double*>(mxGetData(prhs[2]));
    mwSize dims[2];

    
    dims[0]=mxGetM(prhs[0]);
    dims[1]=1;
    
    double j11,j12,j13,j21,j22,j23,j31,j32,j33;
	
	plhs[0] = mxCreateDoubleMatrix(dims[0], 1, mxREAL);
    double* jacb = static_cast<double*>(mxGetData(plhs[0]));
   
	
	for (int i=0; i<dims[0]; i++)
	{
			 j11=d1[i]; j12=d1[i+dims[0]]; j13=d1[i+2*dims[0]];
             j21=d2[i]; j22=d2[i+dims[0]]; j23=d2[i+2*dims[0]];
             j31=d3[i]; j32=d3[i+dims[0]]; j33=d3[i+2*dims[0]];
            
            j11=j11;j22=j22;j33=j33;
            
            //double a = j11*(j22*j33-j32*j33);
            //double b = j12*(j21*j33-j23*j31);
            
            
			jacb[i]=(j11*(j33*j22-j23*j32))-(j12*(j21*j33-j23*j31))+(j13*(j21*j32-j22*j31));
		
	}
			
};
		