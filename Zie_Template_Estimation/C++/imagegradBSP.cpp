
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>
using namespace std;




void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[])
{
    
    double* pts = static_cast<double*>(mxGetData(prhs[0]));
    double* dataI = static_cast<double*>(mxGetData(prhs[1]));
    double* offset = static_cast<double*>(mxGetData(prhs[2]));
    double* scale = static_cast<double*>(mxGetData(prhs[3]));
    
    
    //double* det = static_cast<double*>(mxGetData(prhs[7]));
    
    int ndim =(int)mxGetNumberOfDimensions(prhs[2]);
    int* dim_pts =(int*)mxGetDimensions(prhs[0]);
    int* dim_img=(int*)mxGetDimensions(prhs[1]);
    int* dim_offset=(int*)mxGetDimensions(prhs[2]);
    int* dim_scale=(int*)mxGetDimensions(prhs[3]);
    int N=dim_pts[0];
    ////Allocate Tc
    mwSize dims[2];
    mwSize dims2[1];
    
    dims2[0] = 1;
    
    
    
    dims[0] = N; dims[1] = 1;
    
    double* diff1 ;
    double* diff2 ;
    double* diff3 ;

    
    
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    double* val = static_cast<double*>(mxGetData(plhs[0]));
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff1 = static_cast<double*>(mxGetData(plhs[1]));
    plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff2 = static_cast<double*>(mxGetData(plhs[2]));
    
    plhs[3] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    diff3 = static_cast<double*>(mxGetData(plhs[3]));
    
    double lval;
   
    
    int s=64;
    int idx_idx=0;
    for(int i=0;i<N;i++){
        
        
        
        diff1[i]=0;
        diff2[i]=0;
        diff3[i]=0;

        double t=(pts[i]-offset[0])/scale[0]-floor((pts[i]-offset[0])/scale[0]);
        double t2=t*t;
        double t3=t2*t;
        
        /*double* x={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
         * double* dx={-pow(t,3)+3*pow(t,2)-3*t+1, 3*pow(t,3)-6*pow(t,2)+4, -3*pow(t,3)+3*pow(t,2)+3*t+1,  pow(t,3)};
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
        

        lval=0;
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                for(int l=0;l<4;l++){
                    
                    int idx=px[l]+dim_img[0]*py[k]+dim_img[0]*dim_img[1]*pz[j];
                    double dfdp=x[l]*y[k]*z[j]/216;
                    lval+=dataI[idx]*dfdp;
                    diff1[i]+=dx[l]*y[k]*z[j]/216*dataI[idx];
                    diff2[i]+=x[l]*dy[k]*z[j]/216*dataI[idx];
                    diff3[i]+=x[l]*y[k]*dz[j]/216*dataI[idx];
                    
                    
                }
            }
        }
        
        val[i]=lval;

    }
};

