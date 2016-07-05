#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <vector>

using namespace std;

class interp {

private:
  double* pts ;
  double* dataI ;
  double* dataR ;
  double* offset ;
  double* range;
  double* scale;
  int ndim ;
  int* dim_pts ;
  int* dim_img;
  int* dim_offset;
  int* dim_scale;
  int N;
  double* no_bins;
  double* pI ;
  double* pR ;
  double* pRI;
  double* det;
  double* sumDet;
  double* dfdx;
  double* diff1 ;
  double* diff2 ;
  double* diff3 ;
  double* val;
  
public:
  interp(double* vpts,double* vdataI,double* vdataR,double* voffset,double* vscale,int vndim,int* vdim_pts,int* vdim_img,int* vdim_offset,int* vdim_scale,
          int vN,mwSize* vdims,double* vrange, double* vno_bins,double* vpI,double* vpR, double* vpRI, double* vdet, double* vsumDet,	double* vdiff1 ,double* vdiff2,	double* vdiff3,double* vval	) {
    pts=vpts;
    dataI=vdataI;
    dataR=vdataR;
    pI=vpI;
    pR=vpR;
    pRI=vpRI;
    no_bins=vno_bins;
    range=vrange;
    offset=voffset;
    scale=vscale;
    ndim=vndim;
    dim_pts=vdim_pts;
    dim_img=vdim_img;
    dim_offset=vdim_offset;
    dim_scale=vdim_scale;
    N=vN;
    det=vdet;
    diff1=vdiff1;
    diff2=vdiff2;
    diff3=vdiff3;
    sumDet=vsumDet;
    val=vval;    
  }
  
  void operator()(int r1, int r2) const {
    // initialize histograms to zero
    double lsumDet=0;
    for(int i=0;i<no_bins[0];i++){
      pI[i]=0;
      for(int j=0;j<no_bins[1];j++){
        pRI[j + i*(int)no_bins[1]]=0;
        pR[j]=0;
      }
    }
    for(int i=r1;i!=r2;i++)
    {
      diff1[i]=0;
      diff2[i]=0;
      diff3[i]=0;
      val[i]=0;
      lsumDet+=1/det[i];
      double t=(pts[i]-offset[0])/scale[0]-floor((pts[i]-offset[0])/scale[0]);
      double t2=t*t;
      double t3=t2*t;
      
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
      
      double valtjeck=0;
      int idxR=(int)floor((dataR[i]-range[2]));
      t=(dataR[i])-floor(dataR[i]);
      t2=t*t;
      t3=t2*t;
      
      double tr_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
      for (int nn=0;nn<4;nn++){
        pR[idxR+nn-1]+=tr_val[nn]/6;
      }
      for(int j=0;j<4;j++){
        for(int k=0;k<4;k++){
          for(int l=0;l<4;l++){
            
            int idx=px[l]+dim_img[0]*py[k]+dim_img[0]*dim_img[1]*pz[j];
            double dfdp=x[l]*y[k]*z[j]/216;
            val[i]+=dataI[idx]*dfdp;
            
            diff1[i]+=dx[l]*y[k]*z[j]/216*dataI[idx];
            diff2[i]+=x[l]*dy[k]*z[j]/216*dataI[idx];
            diff3[i]+=x[l]*y[k]*dz[j]/216*dataI[idx];
            //using Parzen window 3rd order b-spline            
          }
        }
      }
      t=val[i]-floor(val[i]);
      t2=t*t;
      t3=t2*t;
      
      double t_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
      int dd=(int)floor((val[i]-range[0]));

      for(int m=0;m<4;m++)
      {
        pI[dd+m-1]+=1/det[i]*t_val[m]/6;
        for (int nn=0;nn<4;nn++){
          pRI[idxR+nn-1+(int)no_bins[1]*(dd+m-1)]+=1/det[i]*t_val[m]*tr_val[nn]/36;
        }
      }
    } 
    sumDet[0]=lsumDet;    
  }
};


class NMI2 {

private:
  mwSize* dims;
  double* pts ;
  double* dataR ;
  double* dataI ;
  double* offset ;
  double* range;
  double* scale;
  int ndim ;
  int* dim_pts ;
  int* dim_img;
  int* dim_offset;
  int* dim_scale;
  int N;
  double* no_bins;
  double* diff1 ;
  double* diff2 ;
  double* diff3 ;
  double* pI ;
  double* pR;
  double* pRI ;
  double* LogpI ;
  double* LogpRI;
  double* det;
  double* val;
  double HI,HR,HRI;
  
  
public:
  NMI2(double* vpts,double* vdataR,double* vdataI,double* voffset,double* vscale,int vndim,int* vdim_pts,int* vdim_img,int* vdim_offset,int* vdim_scale,
          int vN,	double* vdiff1 ,double* vdiff2,	double* vdiff3 ,mwSize* vdims,double vHI,double vHR,double vHRI,double* vrange, double* vno_bins,double* vLogpI, double* vLogpRI,double* vdet,double* vval) {
    HR=vHR;
    HI=vHI;
    HRI=vHRI;
    LogpI=vLogpI;
    LogpRI=vLogpRI;
    range=vrange;
    no_bins=vno_bins;
    pts=vpts;
    dataR=vdataR;
    dataI=vdataI;
    offset=voffset;
    scale=vscale;
    ndim=vndim;
    dim_pts=vdim_pts;
    dim_img=vdim_img;
    dim_offset=vdim_offset;
    dim_scale=vdim_scale;
    N=vN;
    diff1=vdiff1;
    diff2=vdiff2;
    diff3=vdiff3;
    det=vdet;
    val=vval;
    
    dims=vdims;
  }
  
  void operator()(int r1, int r2) const {
    
    for(int i=r1;i!=r2;i++) {
      int idxR=(int)floor((dataR[i]-range[2]));
      double t=(dataR[i])-floor(dataR[i]);
      double t2=t*t;
      double t3=t2*t;
      double tr_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
      t=val[i]-floor(val[i]);
      t2=t*t;
      t3=t2*t;
      int dd=(int)(floor(val[i])-range[0]);
      double t_val[4]={-t3+3*t2-3*t+1, 3*t3-6*t2+4, -3*t3+3*t2+3*t+1,  t3};
      double dt_val[4]={-3*t2+6*t-3, 9*t2-12*t, -9*t2+6*t+3,  3*t2};
      double dNMIdW=0;
      for(int m=0;m<4;m++) {
        for (int nn=0;nn<4;nn++){
          dNMIdW+=1/det[i]*((-(LogpI[dd+m-1]+1)-(HR+HI)*(-(LogpRI[idxR+nn-1+(int)no_bins[1]*(dd+m-1)]+1))/HRI)/(HRI*N))*dt_val[m]*tr_val[nn]/36;
          //dNMIdW+=1/det[i]*((-(LogpI[dd+m-1]+1)-(HR+HI)*(-(LogpRI[idxR+nn-1+(int)no_bins[1]*(dd+m-1)]+1)*(tr_val[nn]/6))/HRI)/(HRI*N))*dt_val[m]/6;
        }
      }
      diff1[i]=dNMIdW*diff1[i];
      diff2[i]=dNMIdW*diff2[i];
      diff3[i]=dNMIdW*diff3[i];
    }
    
  }
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  if(nrhs==0)
    mexPrintf("NMI3D takes 6 arguments,pts,ref_data, Image, range, nobins, offset, scale (optional) ");
  if(nrhs<6)
    mexErrMsgTxt("Number of arguments must be 6");

  double* pts = mxGetPr(prhs[0]);
  double* dataR = mxGetPr(prhs[1]);
  double* dataI = mxGetPr(prhs[2]);
  double* range = mxGetPr(prhs[3]);
  double* no_bins = mxGetPr(prhs[4]);
  double* offset = mxGetPr(prhs[5]);
  double* scale=new double[2];
  scale[0]=1;
  scale[1]=1;
  
  scale = mxGetPr(prhs[6]);
  double* det = mxGetPr(prhs[7]);
  double detSum=0;
  double* pI=new double[(int)no_bins[0]];
  double* pR=new double[(int)no_bins[1]];
  double* pRI=new double[(int)(no_bins[0]*no_bins[1])];
  double* LogpI=new double[(int)no_bins[1]];
  double* LogpRI=new double[(int)(no_bins[0]*no_bins[1])];
  
  int ndim =(int)mxGetNumberOfDimensions(prhs[2]);
  int* dim_pts =(int*)mxGetDimensions(prhs[0]);
  int* dim_img=(int*)mxGetDimensions(prhs[2]);
  int* dim_offset=(int*)mxGetDimensions(prhs[5]);
  int* dim_scale=(int*)mxGetDimensions(prhs[6]);
  int N=dim_pts[0];
  mwSize dims[2];
  mwSize dims2[1];
  
  dims2[0] = 1;
  dims[0] = N; dims[1] = 1;
  double* vals=new double[N];
  double* diff1 ;
  double* diff2 ;
  double* diff3 ;
  double HRI=0;
  double HI=0;
  double HR=0;
  bool do_deriv=false;
  if (nlhs>3) 
    do_deriv = true;
  plhs[0] = mxCreateNumericArray(1,dims2,mxDOUBLE_CLASS, mxREAL);
  double* val = mxGetPr(plhs[0]);
  plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
  diff1 = mxGetPr(plhs[1]);
  plhs[2] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
  diff2 = mxGetPr(plhs[2]);
  
  plhs[3] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
  diff3 = mxGetPr(plhs[3]);

  interp interpO = interp(pts,dataI,dataR,offset,scale,ndim,dim_pts,dim_img,dim_offset,dim_scale,N,dims,range,no_bins,pI,pR,pRI,det,&detSum,diff1,diff2,diff3,vals);
  interpO(0,N);

  double sumPI;
  
  //take the logarithm
  sumPI=0; // for checking normalizations
  int idx;
  double pIlog1, pInormalize, pRInormalize, pRIlog1, pIregul, pRIregul;

  /* Each histogram bin is added a regularization value so that the sum of these matches regul times the count of 
   * entries put into the histogram. However, the addition to the pI histogram should at most be one vote per bin
   * so the regularization to image votes ratio will be lower with many image votes (= less regularization needed). */    
  double alpha=1; 
  /* Overall factor on histogram regularizations. alpha 1 and 0.01 gives 0.7% and 2.0% cross-platform difference on Single101. 
   * This dictated the choice of 1, that was tested to give same Cohen's D. Possibly higher alpha would work even better */
  pIregul = N/no_bins[0];
  pRIregul = N/(no_bins[0]*no_bins[1]);
  if (pIregul>1) {
    pRIregul /= pIregul;
    pIregul = 1;
  }
  pRIregul *= alpha;
  pIregul *= alpha;
  
  // A few helper variables to speed things up below - avoids many log calls (of same value) 
  pInormalize = 1/(detSum+no_bins[0]*pIregul);
  pIlog1 = log(pIregul * pInormalize);
  pRInormalize = 1/(detSum+no_bins[0]*no_bins[1]*pRIregul);
  pRIlog1 = log(pRIregul * pRInormalize);
  // Loop across histograms and compute NMI components
  for(int i=0;i<no_bins[0];i++){
    if(pI[i]==0){
      pI[i] = pIregul * pInormalize;
      LogpI[i] = pIlog1;
    } else {
      pI[i] = (pI[i]+pIregul) * pInormalize;
      LogpI[i] = log(pI[i]);
    }
    HI+=pI[i]*LogpI[i];
    sumPI+=pI[i];
    for(int j=0;j<no_bins[1];j++){
      idx = j + i*(int)no_bins[1];
      if(pRI[idx]==0){
        pRI[idx] = pRIregul * pRInormalize;
        LogpRI[idx] = pRIlog1;
      } else {
        pRI[idx] = (pRI[idx]+pRIregul) * pRInormalize;
        LogpRI[idx] = log(pRI[idx]);
      }
      HRI += pRI[idx]*LogpRI[idx];
    }
  }

  // printf("Sum pI %f, detSum %f , N %d, nobins %f %f\n",sumPI, detSum, N, no_bins[0], no_bins[1]);
  
  for(int j=0;j<no_bins[1];j++){    
    if(pR[j]==0){
      pR[j] = pIregul * pInormalize;
      HR += pR[j]*pIlog1;
    } else {
      pR[j] = (pR[j]+pIregul) * pInormalize;
      HR += pR[j]*log(pR[j]);
    }
  }
  
  //compute the spatial derivatives of NMI
  NMI2 nmi2O = NMI2(pts,dataR,dataI,offset,scale,ndim,dim_pts,dim_img,dim_offset,dim_scale,N,diff1 ,diff2,diff3 ,dims,HI,HR,HRI,range, no_bins,LogpI,LogpRI,det,vals);
  nmi2O(0,N);
  val[0]=((-HR)+(-HI))/(-HRI);
  
  delete pI;
  delete pR;
  delete pRI;
  delete LogpI;
  delete LogpRI;
  delete vals;
  return;
};
