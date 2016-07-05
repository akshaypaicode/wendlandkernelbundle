// wendsqnorm computes ||f||^2 as sum of kernel products.cpp : Defines the entry point for the console application.
//#include "stdafx.h"
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "Kroenecker.h"



//using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[])
{
    
    
    
    double* data = static_cast<double*>(mxGetData(prhs[0]));
    double* spacing = static_cast<double*>(mxGetData(prhs[1]));
    double* sup = static_cast<double*>(mxGetData(prhs[2]));
    double* scaling = static_cast<double*>(mxGetData(prhs[3]));
    
    
    int support = static_cast<int>(sup[0]);
    
    
    
    //Get dimension of data
    int* dim_data=(int*)mxGetDimensions(prhs[0]);
    int N_cts = dim_data[0]*dim_data[1]*dim_data[2];
    
    
    
    mwSize dims_def[1];
    dims_def[0]=1;
//
//
    mwSize dims_der[2];
    dims_der[0]=N_cts;
    dims_der[1]=3;
//
//
    plhs[0] = mxCreateNumericArray(1,dims_def,mxDOUBLE_CLASS, mxREAL);
    double* val = static_cast<double*>(mxGetData(plhs[0]));
    
    plhs[1] = mxCreateNumericArray(2,dims_der,mxDOUBLE_CLASS, mxREAL);
    double* drdp = static_cast<double*>(mxGetData(plhs[1]));
    
    
    mwSize dims_mat[2];
    dims_mat[0]=dim_data[0]*dim_data[1]*dim_data[2];
    dims_mat[1]=dim_data[0]*dim_data[1]*dim_data[2];
    
    plhs[2] = mxCreateNumericArray(2,dims_mat,mxDOUBLE_CLASS, mxREAL);
    double* distMat = static_cast<double*>(mxGetData(plhs[2]));
    
    plhs[3] = mxCreateNumericArray(2,dims_mat,mxDOUBLE_CLASS, mxREAL);
    double* Dval = static_cast<double*>(mxGetData(plhs[3]));
    
    int support_ind[support];
    
    for (int s=0;s<=support;s++)
    {
        support_ind[s]=((-support/2))+s;
    }
    
    val[0]=0;
    double derx,dery,derz;
    
    
    
    for (int i=0;i<dim_data[2];i++)
    {
        for (int j=0; j<dim_data[1]; j++)
        {
            for (int k=0; k<dim_data[0]; k++)
            {
                int i_index=(k)+(dim_data[0]*(j))+(dim_data[1]*dim_data[0]*(i));
                
                //mexPrintf("%d\n",i_index);
                
                drdp[i_index]=0;
                drdp[i_index+N_cts]=0;
                drdp[i_index+2*N_cts]=0;
                
                double ddx=0;
                double ddy=0;
                double ddz=0;
                
                for (int l=0;l<=support;l++)
                {
                    for (int m=0;m<=support;m++)
                    {
                        for (int n=0;n<=support;n++)
                        {
                            
                            int a = k+support_ind[n];//(int)max<double>(min<double>(k+support_ind[n],dim_data[2]-1),0);
                            int b = j+support_ind[m];//(int)max<double>(min<double>(j+support_ind[m],dim_data[1]-1),0);
                            int c = i+support_ind[l];//(int)max<double>(min<double>(i+support_ind[l],dim_data[0]-1),0);
                            
                            if ((a<=dim_data[0]-1 & a>=0) & ((b<=dim_data[1]-1 & b>=0)) & ((c<=dim_data[2]-1 & c>=0)))
                            {
                                
                                int linearindex=(a)+dim_data[0]*(b)+dim_data[0]*dim_data[1]*(c);
                                //mexPrintf("%f\n",data[i_index+N_cts]);
                                
                                
                                double x = support_ind[n];
                                double y = support_ind[m];
                                double z = support_ind[l];
                                
                                double r=sqrt(x*x+y*y+z*z);
                                
                                double su=(support/2);
                                r=r*scaling[0]/(su);
                                
                                
                                double p = max<double>(1-r,0);
                                
                                int cc =  i_index * (dim_data[0]*dim_data[1]*dim_data[2]) + linearindex;
                                
                                
                                double dr = 0;
                                
                                if (su>1)
                                {
                                    
                                    if (r>0 && r<=1)
                                    {
                                        dr = (4)*(p*p*p*-1*(4*r+1)+(p*p*p*p));
                                    }
                                    else
                                    {
                                        dr = 0;
                                    }
                                }
                                
                                else
                                    
                                {
                                    if (r>0 && r<=1)
                                    {
                                        dr = -1;
                                    }
                                    else
                                    {
                                        dr = 0;
                                    }
                                }

                                Dval[cc] = dr*(-r/su);
                                
                                if (su>1)
                                {
                                    val[0] += ((data[i_index]*data[linearindex])+(data[i_index+N_cts]*data[linearindex+N_cts])+(data[i_index+2*N_cts]*data[linearindex+2*N_cts]))*(p*p*p*p*(4*r+1));
                                    distMat[cc] = (p*p*p*p*(4*r+1));
                                }
                                else
                                {
                                    val[0] += ((data[i_index]*data[linearindex])+(data[i_index+N_cts]*data[linearindex+N_cts])+(data[i_index+2*N_cts]*data[linearindex+2*N_cts]))*p;
                                    distMat[cc] = p;
                                }
                                
                                if (support_ind[l]==0 & support_ind[m]==0 & support_ind[n]==0)
                                {
                                    
                                    
                                    if (su>1)
                                    {
                                        derx=2*p*p*p*p*((4*r)+1)*data[linearindex];
                                        dery=2*p*p*p*p*((4*r)+1)*data[linearindex+N_cts];
                                        derz=2*p*p*p*p*((4*r)+1)*data[linearindex+2*N_cts];
                                    }
                                    else
                                    {
                                        derx=2*p*data[linearindex];
                                        dery=2*p*data[linearindex+N_cts];
                                        derz=2*p*data[linearindex+2*N_cts];
                                    }
                                }
                                else
                                {
                                    
                                    if (su>1)
                                    {
                                        derx=p*p*p*p*((4*r)+1)*data[linearindex];
                                        dery=p*p*p*p*((4*r)+1)*data[linearindex+N_cts];
                                        derz=p*p*p*p*((4*r)+1)*data[linearindex+2*N_cts];
                                    }
                                    else
                                    {
                                        derx=p*data[linearindex];
                                        dery=p*data[linearindex+N_cts];
                                        derz=p*data[linearindex+2*N_cts];
                                    }
                                }
                                
                                ddx+=derx;
                                ddy+=dery;
                                ddz+=derz;
                                
                            }
                            
                        }
                        
                        
                        
                    }
                }
                drdp[i_index]+=ddx;
                drdp[i_index+N_cts]+=ddy;
                drdp[i_index+2*N_cts]+=ddz;
                
                
            }
        }
    }
    
};

