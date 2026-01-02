#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mathfunc.cuh"

__global__ void findPressureTimeHistory(float* rObs,float* tObs,float* pT,float* pL,float* pA,float* Xham,float* Nham,float* Vham,float* PRham,int nSurf,int* XPtr,int* NPtr,int* VPtr,int* PRPtr,int* E,int nObs,int nTime, float dPsi,float dTau,float a0,float rhoRef,float* vInf, float oM,int *errorFlag);

__device__ float acous_dist(float* src,float* obs,float* vinf,float a0);
__device__ float sumprod(float* vec1,float* vec2);
__device__ float vecmag(float* vec);
__device__ void findpts(float* x,float xn,int* ipts,int nx,int* npi,float lpi);
__device__ int bsearch(float* A,float key, int imin, int imax);

__device__ static int nX = 3, nXham = 4, nNham = 6, nVham = 6, nPRham = 2;
__device__ static int iPR = 0,iPRd = 1,iVX = 0,iVdX = 3,iNX = 0,iNdX = 3,iDS = 3,iX = 0;

__global__ void findPressureTimeHistory(float* rObs,float* tObs,float* pT,float* pL,float* pA,float* Xham,float* Nham,float* Vham,float* PRham,int nSurf,int* XPtr,int* NPtr,int* VPtr,int* PRPtr,int* E,int nObs,int nTime, float dPsi,float dTau,float a0,float rhoRef,float* vInf, float oM,int *errorFlag){

   float tau;
   int iXham;
   int iNham;
   int iVham;
   int iPRham;
   float obs[3],src[3],V[3],U[3],Vdot[3],Udot[3],rvec[3],rhat[3];
   float rho,rhodot,P,Pdot,N[3],Ndot[3],ds;
   float r,Mr,M,Mdot_r,Vn,Vn_dot,Un,Un_dot;
   float Li[3],Lidot[3],Lidot_r,Li_r,Li_m;
   float Mo[3],Modot[3],Modot_r,Mo_r,Mo_m,Ma_n,Man_dot;
   float D1,T1,T2,T3,T4;
   float LN,LF1,LF2,TN,TF2,xx,zz;
   float pThick,pLoad,tObs2;
   // Binning variables
   int npi;
   int mm,mp,nm,np;	 
   float lpi = nTime*dPsi; //2*M_PI;
   int ipts,iptp,iptm;
   float wtm,wtp,ttemp,ttem,ttep;
  
   int itau = blockIdx.x*blockDim.x + threadIdx.x; // time step 

   while(itau<nTime){ 

   int iobs = blockIdx.z*blockDim.z + threadIdx.z; // Observer index  

   //int iobs0 = 1; 
 

   while(iobs<nObs){

       obs[0] = rObs[iobs*nX+0]; obs[1] = rObs[iobs*nX+1]; obs[2] = rObs[iobs*nX+2];
       D1 = 1.0/(4.0*M_PI*a0);
 
       
              
       for (int i=0;i<nSurf;i++)
       {   
   	   int j = blockIdx.y*blockDim.y + threadIdx.y; // elements

 
   //    	   printf("Testing1 %d\n",E[0]);


	   while (j<E[i])
    	   {

 
	//	if(j%1000==0) printf("BlockID: %d %d %d; ThreadID %d %d %d; iobs:%d itau:%d iel:%d E:%d\n",blockIdx.x,blockIdx.y,blockIdx.z,threadIdx.x,threadIdx.y,threadIdx.z,iobs,itau,j,E[0]);
   

		iXham = XPtr[i] + j*nTime*nXham + itau*nXham;
		iNham = NPtr[i] + j*nTime*nNham + itau*nNham;
		iVham = VPtr[i] + j*nTime*nVham + itau*nVham;
		iPRham = PRPtr[i] + j*nTime*nPRham + itau*nPRham;

		tau = itau*dTau;

		ds = Xham[iXham+iDS];
		rho = 1.0;
		rhodot = 1.0;
		P = PRham[iPRham + iPR];
		Pdot = PRham[iPRham + iPRd];

		for (int id=0;id<nX;id++) V[id]    = Vham[iVham + iVX + id];
		for (int id=0;id<nX;id++) Vdot[id] = Vham[iVham + iVdX + id];
        	for (int id=0;id<nX;id++) N[id]    = Nham[iNham+iNX+id];
		for (int id=0;id<nX;id++) Ndot[id] = Nham[iNham+iNdX+id];
		for (int id=0;id<nX;id++) src[id]  = Xham[iXham + iX + id];



		r = acous_dist(src,obs,vInf,a0);

        	for (int id=0;id<nX;id++) rvec[id] = obs[id] - src[id] - vInf[id]*r/a0;
        	for (int id=0;id<nX;id++) rhat[id] = rvec[id]/r;

		// Special modification motion only around x-axis
                //  Ndot[0] = 0.0;
                // End modification

        	Mr     = sumprod(V,rhat)/a0;
        	M      = vecmag(V)/a0;
        	Mdot_r = sumprod(Vdot,rhat)/a0;

        	T1 = D1/(r*pow((1.0 - Mr),2));
        	T2 = T1/(1.0 - Mr);
        	T3 = T1/r;
        	T4 = T3/(1.0 - Mr);

        	Vn = sumprod(V,N);
        	Vn_dot = sumprod(Vdot,N) + sumprod(V,Ndot);

   		int onsurf = 1; // on-surface calculation (impermeable version)
		// set surface velocity equal to flow velocity 
		if (onsurf==1)
		{
		   for (int id=0;id<nX;id++) U[id] = V[id];
		   for (int id=0;id<nX;id++) Udot[id] = Vdot[id];
		}

		Un = sumprod(U,N);
		Un_dot = sumprod(Udot,N) + sumprod(U,Ndot);
		
		// Begin Modification
		 Vn_dot = 0.0; // If Vn is supposed to be a constant 
		 Un_dot = 0.0;
		// End Modification

		Ma_n = (Un - Vn)*rho + rhoRef*Vn;
		Man_dot = (Un_dot - Vn_dot)*rho + (Un - Vn)*rhodot + rhoRef*Vn_dot;
		TF2 = Man_dot*T1*a0*ds + Ma_n*Mdot_r*T2*a0*ds;
		TN = a0*Ma_n*(Mr - M*M)*T4*a0*ds;
	      

		// thickness pressure	
        	pThick = TN + TF2;

		for (int id=0;id<nX;id++) Li[id] = P*N[id];
		for (int id=0;id<nX;id++) Lidot[id] = Pdot*N[id] + P*Ndot[id];
		  
		Lidot_r = sumprod(Lidot,rhat);
		Li_r = sumprod(Li,rhat);
		Li_m = sumprod(Li,V)/a0;

		// load far		
		LF1 = Lidot_r*T1*ds + Li_r*Mdot_r*T2*ds;
	
		// load near
		LN = (Li_r-Li_m)*T3*a0*ds + Li_r*(Mr-M*M)*T4*a0*ds;

		for (int id=0;id<nX;id++) Mo[id] = rhoRef*U[id]*(Un-Vn);
		for (int id=0;id<nX;id++) Modot[id] = rhoRef*Udot[id]*(Un-Vn) + rhoRef*U[id]*(Un_dot-Vn_dot);
		Modot_r = sumprod(Modot,rhat);
		Mo_r = sumprod(Mo,rhat);
		Mo_m = sumprod(Mo,V)/a0;

		LF2 = Modot_r*T1*ds + Mo_r*Mdot_r*T2*ds;

		// load near
		LN = LN + (Mo_r-Mo_m)*T3*a0*ds + Mo_r*(Mr-M*M)*T4*a0*ds;
	   		
		pLoad =  LN + LF1 + LF2; // loading noise

		ttemp = (tau + r/a0)*oM;


		//if(j==10 && iobs ==0 && i==0){
		//
		//		//fprintf(fid,"%d %e %e\n",k+360*iChunk,v5[k],v5d[k]);
		// 	printf("%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",itau,P,Pdot,LN,LF1,LF2,Modot_r,Mo_r,T1,T2,Mo_m,T3,Mr,M,T4,Lidot_r,Li_r);
		//
		//	
		//}

		// Binning Approach //

	
	        findpts(tObs,ttemp,&ipts,nTime,&npi,lpi);
 	
		// compute previous and next time step ttemp
	
		if(itau>0){	
			mm = itau - 1;
			nm = 0;
		} else {
			mm = nTime - 1;
			nm = 1;
		}
		if(itau<nTime-1){	
			mp = itau + 1;
			np = 0;
		} else {
			mp = 0;
			np = 1;
		}
		iXham = XPtr[i] + j*nTime*nXham + mm*nXham;	
		for (int id=0;id<nX;id++) src[id]  = Xham[iXham + iX + id];
		r = acous_dist(src,obs,vInf,a0);
		ttem = ((mm*dTau) + r/a0)*oM - nm*lpi; //2*M_PI; 
		iXham = XPtr[i] + j*nTime*nXham + mp*nXham;	
		for (int id=0;id<nX;id++) src[id]  = Xham[iXham + iX + id];
		r = acous_dist(src,obs,vInf,a0);
		ttep = ((mp*dTau) + r/a0)*oM + np*lpi; //2*M_PI;

		iptm = ipts;
		tObs2 = tObs[iptm];
		
		while(tObs2>ttem+npi*lpi){
			wtm = (ttem+npi*lpi - tObs2)/(ttem - ttemp);
	
			atomicAdd(&pT[iptm*nObs + iobs],pThick*wtm);
			atomicAdd(&pL[iptm*nObs + iobs],pLoad*wtm);
			atomicAdd(&pA[iptm*nObs + iobs],(pThick+pLoad)*wtm);
		
			tObs2 -= dPsi;
			iptm -= 1;
			if(iptm<0) iptm = nTime - 1;

	//		if(pT[iptm*nObs+iobs] > 1){
	//			printf("iptm %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",j,itau,P,Pdot,LN,LF1,LF2,Modot_r,Mo_r,T1,T2,Mo_m,T3,Mr,M,T4,Lidot_r,Li_r,TN,TF2,Ma_n,ds,a0,Man_dot,Mdot_r,pThick,pLoad);
	//		*errorFlag = 1;
	//		return;
       	//		}
		
		}
		
		if(ipts<nTime-1){
                        iptp = ipts+1;
                }else{
                        iptp = 0;
                }

		tObs2 = tObs[iptp];
		if(iptp==0) tObs2 += lpi; //2*M_PI;
		while(tObs2<ttep+npi*lpi){
			wtp = (ttep+npi*lpi - tObs2)/(ttep - ttemp);

   			atomicAdd(&pT[iptp*nObs + iobs],pThick*wtp);
			atomicAdd(&pL[iptp*nObs + iobs],pLoad*wtp);
			atomicAdd(&pA[iptp*nObs + iobs],(pThick+pLoad)*wtp);

			tObs2 += dPsi;
			iptp += 1;
			if(iptp>nTime-1) iptp = 0;
		
	//		if(pT[iptp*nObs+iobs] > 1){
	//			printf("iptp %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",j,itau,P,Pdot,LN,LF1,LF2,Modot_r,Mo_r,T1,T2,Mo_m,T3,Mr,M,T4,Lidot_r,Li_r,TN,TF2,Ma_n,ds,a0,Man_dot,Mdot_r,pThick,pLoad);
	//			*errorFlag = 1;
	//			return;
	//		}
		}


		j += gridDim.y*blockDim.y;

		__syncthreads();

	   } // loop for grid points closed

       } // close surface loop

       iobs += gridDim.z*blockDim.z;

       __syncthreads();
	
   }// iobs 

   itau += gridDim.x*blockDim.x;

   __syncthreads();	

   } // itau

} // close function loop

__device__ float acous_dist(float* src,float* obs,float* vinf,float a0){

   float a,b,c,disc,diff_t;

   a = pow(vinf[0],2) + pow(vinf[1],2) + pow(vinf[2],2) - pow(a0,2);
   b = -2*((obs[0] - src[0])*vinf[0] + (obs[1] - src[1])*vinf[1] + (obs[2] - src[2])*vinf[2]);
   c = pow(obs[0] - src[0],2) + pow(obs[1] - src[1],2) + pow(obs[2] - src[2],2);
   disc = pow(b,2) - 4*a*c;
   if (disc<0.0){
      printf("Error in retarded time!!! %f %f %f\n",a,b,c);
      //exit(1);
   }
   diff_t = (-b - sqrt(disc))/(2.0*a);
   if (diff_t<0.0){
      diff_t = (-b + sqrt(disc))/(2.0*a);
   }
   disc = diff_t*a0;
   return disc;
}

__device__ float sumprod(float* vec1,float* vec2){
   float sp = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]; 
   return sp;
}
__device__ float vecmag(float* vec){
   float vm = sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2)); 
   return vm;
}

__device__ void findpts(float* x,float xn,int* ipts,int nx,int* npi,float lpi){

     int n,inext;
     float newx;


     // Shift xn so that it lies within the t array if t is periodic

     *npi = floor((lpi + x[0] - xn)/ lpi);

     newx = xn + *npi* lpi;

     n = bsearch(x,newx,0,nx);



      // This condition should never be reached
     if (n==-1){
        printf("Something wrong in the bsearch, check\n");
        //exit(1);
     }

     // Assign lower part of array
     for (int i=0;i>=0;i--){
        inext = n  + i;
        if (inext<0){
           //if (periodic==1){
              inext = inext%nx;
        }
        *ipts = inext;
     }

}

// Binary search, used in Lagrange interpolation
__device__ int bsearch(float* A,float key, int imin, int imax){

   int imid;
   int nA = imax;

   while (1==1)
   {
      imid = (imin + imax)/2;
      if (key<A[imid]){
         imax = imid;
      }
      else{
         imin = imid;
      }
      if (imin+1>=imax) break;
   }


   if (imin<=nA && A[imin]<=key)
   {
      imid = imin;
   }
   else
   {
      imid = -1;
   }

   return imid;


}


//                                                      -dj-

