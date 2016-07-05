// wendland composition.cpp : Defines the entry point for the console application.
//
// spacing can be anything.
// Input:- points - points to be interpolated
//         data - data at sampling points
//         spacing - spacing of the control points: the boundaries are replicated for support outside the limits. The spacing is assumed to start from 1;
//         dim_data - dimensions of the data - again, important for indexing.
// First form: Akshay Pai, University of Copenhagen - 2015
//#include "stdafx.h"
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Kroenecker.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double* pts = static_cast<double*>(mxGetData(prhs[0]));
    double* spacing = static_cast<double*>(mxGetData(prhs[2]));
    double* sup = static_cast<double*>(mxGetData(prhs[3]));
    double* numComp = static_cast<double*>(mxGetData(prhs[4]));
    double* szParam = static_cast<double*>(mxGetData(prhs[5]));
    double* numelP = static_cast<double*>(mxGetData(prhs[6]));
    double* levels = static_cast<double*>(mxGetData(prhs[7]));
    double* pLevelS = static_cast<double*>(mxGetData(prhs[8]));
    double* voxelSize = static_cast<double*>(mxGetData(prhs[9]));
    
    int pLevels = static_cast<int>(pLevelS[0]);
    
    
    
    const mxArray *data;
    data = prhs[1];
    
//get number of points
    int* dim_pts =(int*)mxGetDimensions(prhs[0]);
    int N_pts=dim_pts[0];
    
//support of the level needed
    
    int supportLast = (int) sup [pLevels-1];
    
    mwSize dims_def[2];
    dims_def[0]=N_pts;
    dims_def[1]=3;
    
    
    
    plhs[0] = mxCreateNumericArray(2,dims_def,mxDOUBLE_CLASS, mxREAL);
    double* newPts = static_cast<double*>(mxGetData(plhs[0]));
    
////define variables
    double px_i,py_i,pz_i,dd,jac3z,x,y,z,r1,dx,dy,dz,dr,v11,v12,v13,pp;
    double* currentData;
    vector< vector <double> > dfdp (supportLast*supportLast*supportLast*3, vector<double>(3));
    vector< vector <double> > temp (supportLast*supportLast*supportLast*3, vector<double>(3));
    int i_idx_idx_1=0;
    int i_idx_idx_2=0;
    int supportCurrent;
    mxArray*cellElement;
    
    for (int i=0;i<N_pts;i++)
    {
        px_i=pts[i];
        py_i=pts[i+N_pts];
        pz_i=pts[i+(2*N_pts)];
        
        for (int c=0; c<(int)numComp[0];c++)
        {
            
            double tx=0;
            double ty=0;
            double tz=0;
            
            
            
            for (int p=0; p<levels[0]; p++)//kernel bundle levels
            {
                
                
                supportCurrent = (int) sup [p];
                int Nsur = supportCurrent*supportCurrent*supportCurrent*3;
                int numelp = static_cast<int>(numelP[p]);
                
                int szP[3];
                szP[0] = static_cast<int>(szParam[(3*p)]);
                szP[1] = static_cast<int>(szParam[1+(3*p)]);
                szP[2] = static_cast<int>(szParam[2+(3*p)]);
                
                double currentSpacing[3];
                currentSpacing[0] = spacing[(3*p)];
                currentSpacing[1] = spacing[1+(3*p)];
                currentSpacing[2] = spacing[2+(3*p)];
                
                cellElement = mxGetCell(data,p);
                currentData = mxGetPr(cellElement);
                
                
                //spacing left and right
                
                double px_s=px_i+currentSpacing[0];
                double py_s=py_i+currentSpacing[1];
                double pz_s=pz_i+currentSpacing[2];
                
                double index_x=floor((px_s)/currentSpacing[0]);
                double index_y=floor((py_s)/currentSpacing[1]);
                double index_z=floor((pz_s)/currentSpacing[2]);
                
                
                int* a_1 = new int[supportCurrent];
                int* b_1 = new int[supportCurrent];
                int* c_1 = new int[supportCurrent];
                
                double* a_2 = new double[supportCurrent];
                double* b_2 = new double[supportCurrent];
                double* c_2 = new double[supportCurrent];
                
                
                
                for (int s=0;s<sup[p];s++)
                {
                    
                    
                    double ss=((-supportCurrent/2)+1)+s;
                    
                    
                    a_1[s] = (int)min<double>(max<double>(index_x-1+ss,0),szP[0]-1);
                    a_2[s] = min<double>(max<double>(index_x+ss,0),szP[0]-1);
                    b_1[s] = (int)min<double>(max<double>(index_y-1+ss,0),szP[1]-1);
                    b_2[s] = min<double>(max<double>(index_y+ss,0),szP[1]-1);
                    c_1[s] = (int)min<double>(max<double>(index_z-1+ss,0),szP[2]-1);
                    c_2[s] = min<double>(max<double>(index_z+ss,0),szP[2]-1);
                }
                
                for (int j=0; j<supportCurrent;j++)
                {
                    for (int k=0; k<supportCurrent;k++)
                    {
                        for (int l=0; l<supportCurrent;l++)
                            
                        {
                            
                            int linearindex=(a_1[l])+szP[0]*(b_1[k])+szP[1]*szP[0]*(c_1[j]);
                            
                            
                            x = ((px_s)/currentSpacing[0])-(a_2[l]);
                            y = ((py_s)/currentSpacing[1])-(b_2[k]);
                            z = ((pz_s)/currentSpacing[2])-(c_2[j]);
                            
                            r1=sqrt(x*x+y*y+z*z);
                            
                            double su=(double)supportCurrent/2;
                            
                            double r=r1/(su);
                            
                            pp = max<double>(1-r,0);
                            
                            if (su>1)
                                
                            {
                                dd=pp*pp*pp*pp*((4*r)+1);
                            }
                            else
                            {
                                dd=pp;
                            }
                            
                            tx+=dd*currentData[linearindex];
                            ty+=dd*currentData[linearindex+numelp];
                            tz+=dd*currentData[linearindex+2*numelp];
                            
                            
                            
                        }
                    }
                }
                
                
                delete a_1;
                delete b_1;
                delete c_1;
                
                delete a_2;
                delete b_2;
                delete c_2;
                
                
                
            }
            
            px_i+=(tx/numComp[0]);
            py_i+=(ty/numComp[0]);
            pz_i+=(tz/numComp[0]);
            
        }
        
        newPts[i]=px_i;
        newPts[i+N_pts]=py_i;
        newPts[i+(2*N_pts)]=pz_i;
    }
    
    
    
    
};

