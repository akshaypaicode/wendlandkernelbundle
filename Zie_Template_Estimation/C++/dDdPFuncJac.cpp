
//

#include <math.h>
#include <matrix.h>
#include <mex.h>



class interp
{
private:
    //	mwSize* dims;
    double* d ;
    double* dfdp11 ;
    double* dfdp12 ;
    double* dfdp13 ;
    double* dfdp21 ;
    double* dfdp22 ;
    double* dfdp23 ;
    double* dfdp31 ;
    double* dfdp32 ;
    double* dfdp33 ;
    double*	dDdP;
    int* idx ;
    int N;
    int Np;
    int support;
    
    
    
public:
    interp(double* vd,double* vdfdp11,double* vdfdp12,double* vdfdp13,double* vdfdp21,double* vdfdp22,double* vdfdp23,
            double* vdfdp31,double* vdfdp32,double* vdfdp33, double* vdDdP,int* vidx, int vN, int vNp, int vsupport) {
        d=vd;
        dfdp11=vdfdp11;
        dfdp12=vdfdp12;
        dfdp13=vdfdp13;
        dfdp21=vdfdp21;
        dfdp22=vdfdp22;
        dfdp23=vdfdp23;
        dfdp31=vdfdp31;
        dfdp32=vdfdp32;
        dfdp33=vdfdp33;
        dDdP=vdDdP;
        idx=vidx;
        N=vN;
        Np=vNp;
        support=vsupport;
        
        
    }
    
    void operator()(int r1, int r2) const
    {
        
        double* pdDdp=new double[(int)3*Np];
        
        int s=support*support*support;
        int sur=support*support*support*3;
        for(int i=0;i<3*Np;i++)pdDdp[i]=0;
        int j=0;
        int count_a=0;
        int count_b=0;
        //mexPrintf("initialized!!");
        for(int i=r1;i!=r2;i++)
        {
            //j=i>>6;
         if (count_a>(s-1))
         {
             j=j+1;
             count_a=0;
         }
            
         
            
//             pdDdp[idx[i]]+=dfdp11[i]*d[j]+dfdp21[i]*d[j+N]+dfdp31[i]*d[j+2*N];
//             pdDdp[idx[i]+Np]+=dfdp12[i]*d[j]+dfdp22[i]*d[j+N]+dfdp32[i]*d[j+2*N];
//             pdDdp[idx[i]+2*Np]+=dfdp13[i]*d[j]+dfdp23[i]*d[j+N]+dfdp33[i]*d[j+2*N];
            
//             pdDdp[idx[i]]+=dfdp[i]*d[j]+dfdp[i+sur]*d[j+N]+dfdp[i+2*sur]*d[j+2*N];
//             pdDdp[idx[i]+Np]+=dfdp[i+s]*d[j]+dfdp22[i+s+sur]*d[j+N]+dfdp[i]*d[j+2*N];
//             pdDdp[idx[i]+2*Np]+=dfdp[i]*d[j]+dfdp[i]*d[j+N]+dfdp[i]*d[j+2*N];
            
//             
//             
            pdDdp[idx[i]]+=dfdp11[i]*d[j]+dfdp21[i]*d[j+N]+dfdp31[i]*d[j+2*N];
            pdDdp[idx[i]+Np]+=dfdp12[i]*d[j]+dfdp22[i]*d[j+N]+dfdp32[i]*d[j+2*N];
            pdDdp[idx[i]+2*Np]+=dfdp13[i]*d[j]+dfdp23[i]*d[j+N]+dfdp33[i]*d[j+2*N];
            count_a++;
            //mexPrintf("%f\n",dfdp2[i]);
            
            
        }
         for(int j=0;j<3*Np;j++)
      dDdP[j]=dDdP[j]+pdDdp[j];
    
    delete pdDdp;
  }
};




void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[])
{
    
    //bool doDerivative = false;
    //if(nlhs>1)
    //    doDerivative = true;
    //mexPrintf("Spline interpolation takes 4 arguments ");
    double* d = static_cast<double*>(mxGetData(prhs[0]));
    double* dfdp11 = static_cast<double*>(mxGetData(prhs[1]));
    double* dfdp12 = static_cast<double*>(mxGetData(prhs[2]));
    double* dfdp13 = static_cast<double*>(mxGetData(prhs[3]));
    double* dfdp21 = static_cast<double*>(mxGetData(prhs[4]));
    double* dfdp22 = static_cast<double*>(mxGetData(prhs[5]));
    double* dfdp23 = static_cast<double*>(mxGetData(prhs[6]));
    double* dfdp31 = static_cast<double*>(mxGetData(prhs[7]));
    double* dfdp32 = static_cast<double*>(mxGetData(prhs[8]));
    double* dfdp33 = static_cast<double*>(mxGetData(prhs[9]));
    int* idx = static_cast<int*>(mxGetData(prhs[10]));
    double* dim_p = static_cast<double*>(mxGetData(prhs[11]));
    double* sup = static_cast<double*>(mxGetData(prhs[12]));
    
    int support = static_cast<int>(sup[0]);
    int* dim_d =(int*)mxGetDimensions(prhs[0]);
    int N=dim_d[0];
    int Np=(int)dim_p[0];
    ////Allocate Tc
    mwSize dims[2];
    //mexPrintf("%d %f %f",N,Np,dim_p[0]);
    dims[0] = (int)dim_p[0]; dims[1] = 3;
    
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
    double* dDdP = static_cast<double*>(mxGetData(plhs[0]));
    for(int i=0;i<3*Np;i++)dDdP[i]=0;
    
    //parallel_for(blocked_range<int>(0,N*64),interp(d,dfdp1,dfdp2,dfdp3,dDdP,idx,N,Np),auto_partitioner());
    
    interp interpO = interp(d,dfdp11,dfdp12,dfdp13,dfdp21,dfdp22,dfdp23,dfdp31,dfdp32,dfdp33,dDdP,idx,N,Np,support);
    interpO(0,N*(support*support*support));
    
    return;
    
    
    
    
};
