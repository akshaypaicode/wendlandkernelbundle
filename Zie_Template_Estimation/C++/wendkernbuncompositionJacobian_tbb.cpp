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

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/spin_mutex.h"
using namespace std;
using namespace tbb;
typedef spin_mutex MyMutexType;
MyMutexType MyMutex;
MyMutexType MyMutex2;
MyMutexType MyMutex3;
MyMutexType MyMutex4;
class interp
{
private:
    
    double* pts;
    double* spacing;
    double* sup;
    double* numComp;
    double* szParam;
    double* numelP;
    double* levels;
    double* pLevelS;
    double* voxelSize;
    double* do_derivative;
    double* newPts;
    double* diff1;
    double* diff2;
    double* diff3;
    int N_pts;
    const mxArray *data;
    int supportLast;
    
    
public:
    interp(double *vpts,   double* vspacing, double* vsup, double* vnumComp, double* vszParam, double* vnumelP, double* vlevels, double* vpLevelS, double* vvoxelSize, double* vdo_derivative, double* vnewPts, double* vdiff1, double* vdiff2, double* vdiff3, int vN_pts, const mxArray *vdata, int vsupportLast)
    {
        
    pts=vpts;
    spacing=vspacing;
    sup=vsup;
    numComp=vnumComp;
    szParam=vszParam;
    numelP=vnumelP;
    levels=vlevels;
    pLevelS=vpLevelS;
    voxelSize=vvoxelSize;
    do_derivative=vdo_derivative;
    newPts=vnewPts;;
    diff1=vdiff1;
    diff2=vdiff2;
    diff3=vdiff3;
    N_pts=vN_pts;
    data = vdata;
    supportLast = vsupportLast;
    
    }
    
    void operator()(const blocked_range<int> & r) const
	
    
    {
      
        
        
    ///setup idx and val
    int pLevels = static_cast<int>(pLevelS[0]);
    vector< vector <int> > idx (supportLast*supportLast*supportLast, vector<int>(N_pts));
    
    ////define variables
    double px_i,py_i,pz_i,dd,jac3z,x,y,z,r1,dx,dy,dz,dr,v11,v12,v13,pp;
    double jac[3][3];
    double jacNew[3][3];
    double* currentData;
    vector< vector <double> > dfdp (supportLast*supportLast*supportLast*3, vector<double>(3));
    vector< vector <double> > temp (supportLast*supportLast*supportLast*3, vector<double>(3));
    int i_idx_idx_1=0;
    int i_idx_idx_2=0;
    int supportCurrent;
    mxArray*cellElement;
    
    for(int i=r.begin();i!=r.end();i++)
    {
        px_i=pts[i];
        py_i=pts[i+N_pts];
        pz_i=pts[i+(2*N_pts)];
        
        if (do_derivative[0]==1.0)
        {
            //start with zero jacobian
            diff1[i]=1;diff1[i+N_pts]=0;diff1[i+2*N_pts]=0;
            diff2[i]=0;diff2[i+N_pts]=1;diff2[i+2*N_pts]=0;
            diff3[i]=0;diff3[i+N_pts]=0;diff3[i+2*N_pts]=1;
            
        }
        
        
        for (int c=0; c<(int)numComp[0];c++)
        {
            
            double tx=0;
            double ty=0;
            double tz=0;
            
            for (int row=0; row<3;row++)
            {
                for (int col=0; col<3;col++)
                {
                    jacNew[row][col]=0;
                }
            }

            
            
            for (int p=0; p<levels[0]; p++)//kernel bundle levels
            {
                
                
                
                for (int row=0; row<3;row++)
                {
                    for (int col=0; col<3;col++)
                    {
                        jac[row][col]=0;
                    }
                }
            
                
                
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
                            
                            if (su>1)
                            {
                                
                                if (r>0 && r<=1)
                                {
                                    dr = (4)*(pp*pp*pp*-1*(4*r+1)+(pp*pp*pp*pp));
                                    dx=(dr*x/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[0]*currentSpacing[0]*currentSpacing[0]));
                                    dy=(dr*y/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[1]*currentSpacing[1]*currentSpacing[1]));
                                    dz=(dr*z/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[2]*currentSpacing[2]*currentSpacing[2]));
                                }
                                else
                                {
                                    dr = 0;
                                    dx = dr*x;
                                    dy = dr*y;
                                    dz = dr*z;
                                }
                            }
                            
                            else
                                
                            {
                                if (r>0 && r<=1)
                                {
                                    dr = -1;
                                    dx=(dr*x/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[0]*currentSpacing[0]*currentSpacing[0]));
                                    dy=(dr*y/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[1]*currentSpacing[1]*currentSpacing[1]));
                                    dz=(dr*z/sqrt((x*x)+(y*y)+(z*z)))*(1/(voxelSize[2]*currentSpacing[2]*currentSpacing[2]));
                                }
                                else
                                {
                                    dr = 0;
                                    dx = dr*x;
                                    dy = dr*y;
                                    dz = dr*z;
                                }
                            }
                            
                            //compute spatial derivative
                            
                            if (c==0)
                            {
                                
                                jac[0][0]+=0;
                                jac[0][1]+=0;
                                jac[0][2]+=0;
                                
                                jac[1][0]+=0;
                                jac[1][1]+=0;
                                jac[1][2]+=0;
                                
                                jac[2][0]+=0;
                                jac[2][1]+=0;
                                jac[2][2]+=0;
                                
                            }
                            else
                                
                            {
                                
                                jac[0][0]+=currentData[linearindex]*dx;
                                jac[0][1]+=currentData[linearindex+numelp]*dx;
                                jac[0][2]+=currentData[linearindex+2*numelp]*dx;
                                
                                jac[1][0]+=currentData[linearindex]*dy;
                                jac[1][1]+=currentData[linearindex+numelp]*dy;
                                jac[1][2]+=currentData[linearindex+2*numelp]*dy;
                                
                                jac[2][0]+=currentData[linearindex]*dz;
                                jac[2][1]+=currentData[linearindex+numelp]*dz;
                                jac[2][2]+=currentData[linearindex+2*numelp]*dz;
                                
                            }

                            
                        }
                    }
                }
                
                ///add 1 to jacobian
                
                jacNew[0][0]+=(jac[0][0]/numComp[0]);
                jacNew[0][1]+=(jac[0][1]/numComp[0]);
                jacNew[0][2]+=(jac[0][2]/numComp[0]);
                
                jacNew[1][0]+=(jac[1][0]/numComp[0]);
                jacNew[1][1]+=(jac[1][1]/numComp[0]);
                jacNew[1][2]+=(jac[1][2]/numComp[0]);
                
                jacNew[2][0]+=(jac[2][0]/numComp[0]);
                jacNew[2][1]+=(jac[2][1]/numComp[0]);
                jacNew[2][2]+=(jac[2][2]/numComp[0]);
                
                delete a_1;
                delete b_1;
                delete c_1;
                
                delete a_2;
                delete b_2;
                delete c_2;
                
                //setting relevant variables to zero
                i_idx_idx_1=0;
                i_idx_idx_2=0;
                
            }
            
            //add one to the jacobian
            jacNew[0][0] = jacNew[0][0]+1;
            jacNew[1][1] = jacNew[1][1]+1;
            jacNew[2][2] = jacNew[2][2]+1;
            
            
            
            //temp
            
                
            diff1[i]=(diff1[i]*jacNew[0][0])+(diff2[i]*jacNew[0][1])+(diff3[i]*jacNew[0][2]);
            diff1[i+N_pts]=(diff1[i+N_pts]*jacNew[0][0])+(diff2[i+N_pts]*jacNew[0][1])+(diff3[i+N_pts]*jacNew[0][2]);
            diff1[i+2*N_pts]=(diff1[i+2*N_pts]*jacNew[0][0])+(diff2[i+2*N_pts]*jacNew[0][1])+(diff3[i+2*N_pts]*jacNew[0][2]);
            
            diff2[i]=(diff1[i]*jacNew[1][0])+(diff2[i]*jacNew[1][1])+(diff3[i]*jacNew[1][2]);
            diff2[i+N_pts]=(diff1[i+N_pts]*jacNew[1][0])+(diff2[i+N_pts]*jacNew[1][1])+(diff3[i+N_pts]*jacNew[1][2]);
            diff2[i+2*N_pts]=(diff1[i+2*N_pts]*jacNew[1][0])+(diff2[i+2*N_pts]*jacNew[1][1])+(diff3[i+2*N_pts]*jacNew[1][2]);
            
            diff3[i]=(diff1[i]*jacNew[2][0])+(diff2[i]*jacNew[2][1])+(diff3[i]*jacNew[2][2]);
            diff3[i+N_pts]=(diff1[i+N_pts]*jacNew[2][0])+(diff2[i+N_pts]*jacNew[2][1])+(diff3[i+N_pts]*jacNew[2][2]);
            diff3[i+2*N_pts]=(diff1[i+2*N_pts]*jacNew[2][0])+(diff2[i+2*N_pts]*jacNew[2][1])+(diff3[i+2*N_pts]*jacNew[2][2]);

            
            px_i+=(tx/numComp[0]);
            py_i+=(ty/numComp[0]);
            pz_i+=(tz/numComp[0]);
            
            
            
            
        }
        
        newPts[i]=px_i;
        newPts[i+N_pts]=py_i;
        newPts[i+(2*N_pts)]=pz_i;
    }
    
    }
    
};
    
    
    
    









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
    double* do_derivative = static_cast<double*>(mxGetData(prhs[10]));
    
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
    
    //Output jacobian
    mwSize dims_jac[2];
    dims_jac[0]=N_pts;
    dims_jac[1]=3;
    
    plhs[1] = mxCreateNumericArray(2,dims_jac,mxDOUBLE_CLASS, mxREAL);
    double* diff1 = static_cast<double*>(mxGetData(plhs[1]));
    
    plhs[2] = mxCreateNumericArray(2,dims_jac,mxDOUBLE_CLASS, mxREAL);
    double* diff2 = static_cast<double*>(mxGetData(plhs[2]));
    
    plhs[3] = mxCreateNumericArray(2,dims_jac,mxDOUBLE_CLASS, mxREAL);
    double* diff3 = static_cast<double*>(mxGetData(plhs[3]));
    


    parallel_for(blocked_range<int>(0,N_pts),interp(pts, spacing, sup, numComp, szParam, numelP, levels, pLevelS, voxelSize, do_derivative, newPts, diff1, diff2, diff3, N_pts, data, supportLast),auto_partitioner());
};