// SplineInterpolation.cpp : Defines the entry point for the console application.
//
//#include "stdafx.h"
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "Kroenecker.h"



#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>




//using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
const mxArray *prhs[])
{

    //bool doDerivative = false;
    //if(nlhs>1)
    //    doDerivative = true;
    if(nrhs==0)
		mexPrintf("Spline interpolation takes 4 arguments ");
    if(nrhs<3)
        mexErrMsgTxt("Number of arguments must be 3");
    // mexPrintf("Spline interpolation takes 4 arguments ");
    double* pts = static_cast<double*>(mxGetData(prhs[0]));
    double* data = static_cast<double*>(mxGetData(prhs[1]));
	double* offset = static_cast<double*>(mxGetData(prhs[2]));
	double* scale;
	if(nrhs>=4){
        scale = static_cast<double*>(mxGetData(prhs[3]));
    }
    else
    {
        scale=new double[3];
        scale[0]=1;
        scale[1]=1;
        scale[2]=1;
    }

	//int* dx = static_cast<double*>(mxGetData(prhs[2]));
 
	int ndim =(int)mxGetNumberOfDimensions(prhs[1]);
	int* dim_pts =(int*)mxGetDimensions(prhs[0]);
	//mexPrintf("%d",dim_pts[0]);
	int* dim_img=(int*)mxGetDimensions(prhs[1]);
	int* dim_offset=(int*)mxGetDimensions(prhs[2]);
	int* dim_scale=(int*)mxGetDimensions(prhs[3]);
	int N=dim_pts[0];
	//size_t N = mxGetNumberOfElements(prhs[0])/3;
    ////Allocate Tc
    mwSize dims[2];
	mwSize dims2[2];

	dims2[0] = 64; 
	dims2[1] = N;
	
	mwSize dims3[3];

	dims3[0] = 64;
	dims3[1] = 3;
	dims3[2] = N;
	dims[0] = N; dims[1] = 1;
	if(ndim>3){
	dims[1] = dim_img[3];
	}
	double* diff1 ;
    double* diff2 ;
    double* diff3;
    double* dfdx;
    
    double* diffp1 ;
    double* diffp2 ;
    double* diffp3;

	bool do_deriv=false;
	if(nlhs>3)do_deriv=true;
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    double* val = static_cast<double*>(mxGetData(plhs[0]));
	plhs[1] = mxCreateNumericArray(2,dims2,mxDOUBLE_CLASS, mxREAL);
    double* dfdp = static_cast<double*>(mxGetData(plhs[1]));
	plhs[2] = mxCreateNumericArray(2,dims2,mxINT32_CLASS, mxREAL);
    int* idx = static_cast<int*>(mxGetData(plhs[2]));
	
	
	plhs[3] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff1 = static_cast<double*>(mxGetData(plhs[3]));
	plhs[4] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff2 = static_cast<double*>(mxGetData(plhs[4]));
	plhs[5] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff3 = static_cast<double*>(mxGetData(plhs[5]));
	plhs[6] = mxCreateNumericArray(3,dims3,mxDOUBLE_CLASS, mxREAL);
    dfdx = static_cast<double*>(mxGetData(plhs[6]));
    
    
    plhs[7] = mxCreateNumericArray(3,dims3,mxDOUBLE_CLASS, mxREAL);
    diffp1 = static_cast<double*>(mxGetData(plhs[7]));
	plhs[8] = mxCreateNumericArray(3,dims3,mxDOUBLE_CLASS, mxREAL);
    diffp2 = static_cast<double*>(mxGetData(plhs[8]));
	plhs[9] = mxCreateNumericArray(3,dims3,mxDOUBLE_CLASS, mxREAL);
    diffp3 = static_cast<double*>(mxGetData(plhs[9]));
    
   
	int s=64;
    int count=0;
	double* p=new double[s];
	int idx_idx=0;
	for(int i=0;i<N;i++){
		val[i]=0;
		if(do_deriv){
		diff1[i]=0;
		diff2[i]=0;
		diff3[i]=0;
		}
		double t=(pts[i]-offset[0])/scale[0]-floor((pts[i]-offset[0])/scale[0]);
		double t2=t*t;
		double t3=t2*t;

		/*double* x={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
		double* dx={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
		*/
		double x[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dx[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i]-offset[0])/scale[0])-1;
		int px[4]={(int)min<double>(max<double>(t-1,0),dim_img[0]-1),(int) max<double>(min<double>(t,dim_img[0]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[0]-1),(int) max<double>(min<double>(t+2,dim_img[0]-1),0)};
		t=(pts[i+N]-offset[1])/scale[1]-floor((pts[i+N]-offset[1])/scale[1]);
		t2=t*t;
		t3=t2*t;

		double y[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dy[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i+N]-offset[1])/scale[1])-1;
		int py[4]={(int)min<double>(max<double>(t-1,0),dim_img[1]-1),(int) max<double>(min<double>(t,dim_img[1]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[1]-1),(int) max<double>(min<double>(t+2,dim_img[1]-1),0)};
		
		t=(pts[i+2*N]-offset[2])/scale[2]-floor((pts[i+2*N]-offset[2])/scale[2]);
		t2=t*t;
		t3=t2*t;
		double z[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
		double dz[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
		t=floor((pts[i+2*N]-offset[2])/scale[2])-1;
		int pz[4]={(int)min<double>(max<double>(t-1,0),dim_img[2]-1),(int) max<double>(min<double>(t,dim_img[2]-1),0),(int)min<double>(max<double>(t+1,0),dim_img[2]-1),(int) max<double>(min<double>(t+2,dim_img[2]-1),0)};
		
        int index=0;
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
					for(int m=0;m<dims[1];m++){
						//mexPrintf("%f %d %d %d \n",data[px[l]+dim_img[0]*py[k]+dim_img[0]*dim_img[1]*pz[j]+dim_img[0]*dim_img[1]*dim_img[2]*m],px[l],py[k],pz[j]);
				
					if(m<1){
						idx[idx_idx]=px[l]+dim_img[0]*py[k]+dim_img[0]*dim_img[1]*pz[j];
						dfdp[idx_idx]=x[l]*y[k]*z[j]/216;
						idx_idx++;

						}
				double dd=data[idx[idx_idx-1]+dim_img[0]*dim_img[1]*dim_img[2]*m];  //purpose of dd?
						
			
				val[i+m*N]+=dfdp[idx_idx-1]*dd;

					if(do_deriv){
									if(m<1){
							index=i*(192)+l+4*k+16*j;
							dfdx[index]=dx[l]*y[k]*z[j]/216;
							dfdx[index+64]=x[l]*dy[k]*z[j]/216;
							dfdx[index+128]=x[l]*y[k]*dz[j]/216;
                            

							
						}
			
					diff1[i+m*N]+=dfdx[index]*dd; // ???? Jacobian
					diff2[i+m*N]+=dfdx[index+64]*dd;  
					diff3[i+m*N]+=dfdx[index+128]*dd;
                    
                    
                    diffp1[index+64*m]=dfdx[index]*dd;
                    diffp2[index+64*m]=dfdx[index+64]*dd;
                    diffp3[index+64*m]=dfdx[index+128]*dd;                        
//                     diffp2[i+index+m*N]=dfdx[index+64]*dd;
//                     diffp3[i+index+m*N]=dfdx[index+128]*dd;   
                    
          
                    
					}
					}
				}
			}
		}
	
	}
	if(nrhs<4) delete scale;
   
    delete p;
	
};