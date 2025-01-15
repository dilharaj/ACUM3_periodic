#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "variables.cuh"
#include "mathfunc.cuh"

void preprocess_per(){
	
	int iXh;
	int iNh;
	int iVh;
	int iUh;
	int iPRh;
	//float n_x0;
	int n_win=13;
   	//int n_win=25;


	float* v1 = (float*)malloc(nTime*sizeof(float));
	float* v2 = (float*)malloc(nTime*sizeof(float));
	float* v3 = (float*)malloc(nTime*sizeof(float));
	float* v4 = (float*)malloc(nTime*sizeof(float));
	float* v5 = (float*)malloc(nTime*sizeof(float));
	float* v1d = (float*)malloc(nTime*sizeof(float));
	float* v2d = (float*)malloc(nTime*sizeof(float));
	float* v3d = (float*)malloc(nTime*sizeof(float));
	float* v4d = (float*)malloc(nTime*sizeof(float));
	float* v5d = (float*)malloc(nTime*sizeof(float));
	float* n_x = (float*)malloc(nTime*sizeof(float));
	float* n_y = (float*)malloc(nTime*sizeof(float));
	float* n_z = (float*)malloc(nTime*sizeof(float));
	float* nd_x = (float*)malloc(nTime*sizeof(float));
	float* nd_y = (float*)malloc(nTime*sizeof(float));
	float* nd_z = (float*)malloc(nTime*sizeof(float));
	float* v_x = (float*)malloc(nTime*sizeof(float));
	float* v_y = (float*)malloc(nTime*sizeof(float));
	float* v_z = (float*)malloc(nTime*sizeof(float));
	float* vd_x = (float*)malloc(nTime*sizeof(float));
	float* vd_y = (float*)malloc(nTime*sizeof(float));
	float* vd_z = (float*)malloc(nTime*sizeof(float));
	float* x = (float*)malloc(nTime*sizeof(float));
	float* y = (float*)malloc(nTime*sizeof(float));
	float* z = (float*)malloc(nTime*sizeof(float));
	float *f_cubic = (float *)malloc((nTime+1)*sizeof(float));
	float *gauss = (float *)malloc(n_win*sizeof(float));
	
	

	for (int i=0;i<nSurf;i++)
	{
		for(int j=0;j<E[i]*nTime;j++)
		{
			iXh = XPtr[i] + j*nXham;
			iNh = NPtr[i] + j*nNham;
			iUh = UPtr[i] + j*nUham;
			iPRh = PRPtr[i] + j*nPRham;

			Xham[iXh+iXham] = Xham[iXh+iXham]*CFDscale;
			Xham[iXh+iYham] = Xham[iXh+iYham]*CFDscale;	
			Xham[iXh+iZham] = Xham[iXh+iZham]*CFDscale;	
			Xham[iXh+iDSham]= Xham[iXh+iDSham]*CFDscale*CFDscale;
			PRham[iPRh+iRHOham] =  PRham[iPRh+iRHOham]*rhoRef;
			Uham[iUh+iUXham] =  Uham[iUh+iUXham]*a0;
			Uham[iUh+iUYham] =  Uham[iUh+iUYham]*a0;
			Uham[iUh+iUZham] =  Uham[iUh+iUZham]*a0;
			PRham[iPRh+iPRham] =  PRham[iPRh+iPRham]*0.5*rhoRef*pow(Mref,2)*pow(a0,2);

		}
	
	   	//
      		// preprocessor of gaussian filter
      		pre_gaussfilt1(gauss, n_win);
	
		for (int j=0;j<E[i];j++)
		{
			for (int k=0;k<nTime;k++)
			{
				iXh = XPtr[i] + j*nTime*nXham + k*nXham;
				iNh = NPtr[i] + j*nTime*nNham + k*nNham;
				iUh = UPtr[i] + j*nTime*nUham + k*nUham;
				iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;

				v1[k] = PRham[iPRh+iRHOham];
				v2[k] = Uham[iUh+iUXham];
				v3[k] = Uham[iUh+iUYham];
				v4[k] = Uham[iUh+iUZham];
				v5[k] = PRham[iPRh+iPRham];
				// modification - if motion is around x-axis 
				//if(k!=0) Nham[iNh+iNXham] = n_x0;   
				//if(k==0) n_x0 = Nham[iNh+iNXham];
				// modification end
				n_x[k] = Nham[iNh+iNXham];
				n_y[k] = Nham[iNh+iNYham];
				n_z[k] = Nham[iNh+iNZham];
				x[k] = Xham[iXh+iXham];
		                y[k] = Xham[iXh+iYham];
                                z[k] = Xham[iXh+iZham];
		   	}
	         
        		// ==========================================
			// Apply filters
			// =========================================

		
//			       make_periodic(v1, nTime, 0, f_cubic);
//			       make_periodic(v2, nTime, 0, f_cubic);
//			       make_periodic(v3, nTime, 0, f_cubic);
//			       make_periodic(v4, nTime, 0, f_cubic);
//			       make_periodic(v5, nTime, 0, f_cubic);
//			


			if(lowpass==1)
			{
			
			     for(int i=0;i<2;i++){	
			       // ==========================================
			       // median filter
			       // =========================================
			       medfilt1(v1, nTime);
			       medfilt1(v2, nTime);
			       medfilt1(v3, nTime);
			       medfilt1(v4, nTime);
			       medfilt1(v5, nTime);


           		       // ========================================
           		       // Gaussian filter
           		       // =========================================
           		       gaussfilt1(v1, gauss, nTime, n_win);
           		       gaussfilt1(v2, gauss, nTime, n_win);
           		       gaussfilt1(v3, gauss, nTime, n_win);
           		       gaussfilt1(v4, gauss, nTime, n_win);
           		       gaussfilt1(v5, gauss, nTime, n_win);
   			
			     }
			}

			     
			for (int k=0;k<nTime;k++)
			{
			  iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;
			  iUh = UPtr[i] + j*nTime*nUham + k*nUham;
			  PRham[iPRh+iRHOham] = v1[k];
			  Uham[iUh+iUXham]    = v2[k];
			  Uham[iUh+iUYham]    = v3[k];
			  Uham[iUh+iUZham]    = v4[k];
		          PRham[iPRh+iPRham]  = v5[k];
			}
         
	         	
			// Density, Flow velocity and Pressure differentials
			if(periodic == 1){
				differentiate(v1,v1d,dTau,nTime); 
				differentiate(v2,v2d,dTau,nTime);
				differentiate(v3,v3d,dTau,nTime);
				differentiate(v4,v4d,dTau,nTime);
				differentiate(v5,v5d,dTau,nTime);
			} else {
				differentiate_biased(v1,v1d,dTau,nTime); 
				differentiate_biased(v2,v2d,dTau,nTime);
				differentiate_biased(v3,v3d,dTau,nTime);
				differentiate_biased(v4,v4d,dTau,nTime);
				differentiate_biased(v5,v5d,dTau,nTime);
			}	
	
			// Surface normal time derivative
		        differentiate(n_x,nd_x,dTau,nTime);
		        differentiate(n_y,nd_y,dTau,nTime);
		        differentiate(n_z,nd_z,dTau,nTime);
		         

			// Surface velocity
			differentiate(x,v_x,dTau,nTime);
			differentiate(y,v_y,dTau,nTime);
			differentiate(z,v_z,dTau,nTime);

			// Surface acceleration
			differentiate(v_x,vd_x,dTau,nTime);
			differentiate(v_y,vd_y,dTau,nTime);
			differentiate(v_z,vd_z,dTau,nTime);


			
			for (int k=0;k<nTime;k++)
			{
			         
				iXh = XPtr[i] + j*nTime*nXham + k*nXham;
				iNh = NPtr[i] + j*nTime*nNham + k*nNham;
				iVh = VPtr[i] + j*nTime*nVham + k*nVham;
				iUh = UPtr[i] + j*nTime*nUham + k*nUham;
				iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;


				  PRham[iPRh+iRHOdham] = v1d[k];
				  Uham[iUh+iUdXham]  = v2d[k];
				  Uham[iUh+iUdYham]  = v3d[k];
				  Uham[iUh+iUdZham]  = v4d[k];
				  PRham[iPRh+iPRdham]  = v5d[k];
				  Nham[iNh+iNdXham]  = nd_x[k];
				  Nham[iNh+iNdYham]  = nd_y[k];
				  Nham[iNh+iNdZham]  = nd_z[k];
				  Vham[iVh+iVXham]  = v_x[k];
				  Vham[iVh+iVYham]  = v_y[k];
				  Vham[iVh+iVZham]  = v_z[k];
				  Vham[iVh+iVdXham]  = vd_x[k];
				  Vham[iVh+iVdYham]  = vd_y[k];
				  Vham[iVh+iVdZham]  = vd_z[k];
				  Xham[iXh+iXham] = x[k];
				  Xham[iXh+iYham] = y[k];
			          Xham[iXh+iZham] = z[k];
    		       }  
               } 
	
	}

	 free(v1); free(v2); free(v3);free(v4); free(v5);
	 free(v1d); free(v2d); free(v3d); free(v4d); free(v5d); 
	 free(n_x); free(n_y); free(n_z); 
	 free(nd_x); free(nd_y); free(nd_z); 
	 free(v_x); free(v_y); free(v_z);
	 free(vd_x); free(vd_y); free(vd_z);	 
	 free(x); free(y); free(z);
	 free(f_cubic); free(gauss);
}

void preprocess_imper(){
	
	int iXh;
	int iNh;
	int iVh;
	int iPRh;
	//float n_x0;
	int n_win=13;

	float* v1 = (float*)malloc(nTime*sizeof(float));
	float* v2 = (float*)malloc(nTime*sizeof(float));
	float* v3 = (float*)malloc(nTime*sizeof(float));
	float* v4 = (float*)malloc(nTime*sizeof(float));
	float* v5 = (float*)malloc(nTime*sizeof(float));
	float* v1d = (float*)malloc(nTime*sizeof(float));
	float* v2d = (float*)malloc(nTime*sizeof(float));
	float* v3d = (float*)malloc(nTime*sizeof(float));
	float* v4d = (float*)malloc(nTime*sizeof(float));
	float* v5d = (float*)malloc(nTime*sizeof(float));
	float* n_x = (float*)malloc(nTime*sizeof(float));
	float* n_y = (float*)malloc(nTime*sizeof(float));
	float* n_z = (float*)malloc(nTime*sizeof(float));
	float* nd_x = (float*)malloc(nTime*sizeof(float));
	float* nd_y = (float*)malloc(nTime*sizeof(float));
	float* nd_z = (float*)malloc(nTime*sizeof(float));
	float* v_x = (float*)malloc(nTime*sizeof(float));
	float* v_y = (float*)malloc(nTime*sizeof(float));
	float* v_z = (float*)malloc(nTime*sizeof(float));
	float* vd_x = (float*)malloc(nTime*sizeof(float));
	float* vd_y = (float*)malloc(nTime*sizeof(float));
	float* vd_z = (float*)malloc(nTime*sizeof(float));
	float* x = (float*)malloc(nTime*sizeof(float));
	float* y = (float*)malloc(nTime*sizeof(float));
	float* z = (float*)malloc(nTime*sizeof(float));
	float *f_cubic = (float *)malloc((nTime+1)*sizeof(float));
	float *gauss = (float *)malloc(n_win*sizeof(float));
	
	

	for (int i=0;i<nSurf;i++)
	{
		for(int j=0;j<E[i]*nTime;j++)
		{
			iXh = XPtr[i] + j*nXham;
			iNh = NPtr[i] + j*nNham;
			iPRh = PRPtr[i] + j*nPRham;

			Xham[iXh+iXham] = Xham[iXh+iXham]*CFDscale;
			Xham[iXh+iYham] = Xham[iXh+iYham]*CFDscale;	
			Xham[iXh+iZham] = Xham[iXh+iZham]*CFDscale;	
			Xham[iXh+iDSham]= Xham[iXh+iDSham]*CFDscale*CFDscale;
			PRham[iPRh+iPRham] =  PRham[iPRh+iPRham]*0.5*rhoRef*pow(Mref,2)*pow(a0,2);

		}
	
	   	//
      		// preprocessor of gaussian filter
      		pre_gaussfilt1(gauss, n_win);
	
		for (int j=0;j<E[i];j++)
		{
			for (int k=0;k<nTime;k++)
			{
				iXh = XPtr[i] + j*nTime*nXham + k*nXham;
				iNh = NPtr[i] + j*nTime*nNham + k*nNham;
				iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;

				v5[k] = PRham[iPRh+iPRham];
				// modification - if motion is around x-axis 
				//if(k!=0) Nham[iNh+iNXham] = n_x0;   
				//if(k==0) n_x0 = Nham[iNh+iNXham];
				// modification end

	
				n_x[k] = Nham[iNh+iNXham];
				n_y[k] = Nham[iNh+iNYham];
				n_z[k] = Nham[iNh+iNZham];
				x[k] = Xham[iXh+iXham];
		                y[k] = Xham[iXh+iYham];
                                z[k] = Xham[iXh+iZham];
		   	}
	         
			// ==========================================
			// Apply filters
			// =========================================
			


		
//     			        make_periodic(v1, nTime, 0, f_cubic);
//				make_periodic(v2, nTime, 0, f_cubic);
//				make_periodic(v3, nTime, 0, f_cubic);
//				make_periodic(v4, nTime, 0, f_cubic);
				make_periodic(v5, nTime, 0, f_cubic);

			if(lowpass==1)
			{
			
			     for(int i=0;i<2;i++){	
			    	// median filter 
			     	medfilt1(v5, nTime);
				// Gaussian filter
           		       	gaussfilt1(v5, gauss, nTime, n_win);
   			
			     }
			}

			     
			for (int k=0;k<nTime;k++)
			{
			  iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;
		          PRham[iPRh+iPRham] =v5[k];
			}
         

			// Pressure time derivative		
			if(periodic == 1){
				differentiate(v5,v5d,dTau,nTime);
			} else {
				differentiate_biased(v5,v5d,dTau,nTime);
			}
 
			// Surface normals time derivative
		        differentiate(n_x,nd_x,dTau,nTime);
		        differentiate(n_y,nd_y,dTau,nTime);
		        differentiate(n_z,nd_z,dTau,nTime);
		         

			// Surface velocity
			differentiate(x,v_x,dTau,nTime);
			differentiate(y,v_y,dTau,nTime);
			differentiate(z,v_z,dTau,nTime);

			// Surface acceleration
			differentiate(v_x,vd_x,dTau,nTime);
			differentiate(v_y,vd_y,dTau,nTime);
			differentiate(v_z,vd_z,dTau,nTime);

//if(j==10){
//	FILE *fid;
//	
//	fid = fopen("test_peri.dat","a");
//	for(int k=0;k<nTime;k++){
//		//fprintf(fid,"%d %e %e\n",k+360*iChunk,v5[k],v5d[k]);
//		fprintf(fid,"%d %e %e %e %e %e %e %e %e %e %e %e %e %e\n",k,x[k],y[k],z[k],v_x[k],v_y[k],v_z[k],vd_x[k],vd_y[k],vd_z[k],nd_x[k],nd_y[k],nd_z[k],v5d[k]);
//
//	}
//	fclose(fid);
//}
			
			for (int k=0;k<nTime;k++)
			{
			         
				iXh = XPtr[i] + j*nTime*nXham + k*nXham;
				iNh = NPtr[i] + j*nTime*nNham + k*nNham;
				iVh = VPtr[i] + j*nTime*nVham + k*nVham;
				iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;


				PRham[iPRh+iPRdham]  = v5d[k];
				Nham[iNh+iNdXham]  = nd_x[k];
				Nham[iNh+iNdYham]  = nd_y[k];
				Nham[iNh+iNdZham]  = nd_z[k];
				Vham[iVh+iVXham]  = v_x[k];
				Vham[iVh+iVYham]  = v_y[k];
				Vham[iVh+iVZham]  = v_z[k];
				Vham[iVh+iVdXham]  = vd_x[k];
				Vham[iVh+iVdYham]  = vd_y[k];
				Vham[iVh+iVdZham]  = vd_z[k];
				Xham[iXh+iXham] = x[k];
				Xham[iXh+iYham] = y[k];
			        Xham[iXh+iZham] = z[k];
    		       }  
               } 
	
	}

	 free(v1); free(v2); free(v3);free(v4); free(v5);
	 free(v1d); free(v2d); free(v3d); free(v4d); free(v5d); 
	 free(n_x); free(n_y); free(n_z); 
	 free(nd_x); free(nd_y); free(nd_z); 
	 free(v_x); free(v_y); free(v_z);
	 free(vd_x); free(vd_y); free(vd_z);	 
	 free(x); free(y); free(z);
	 free(f_cubic); free(gauss);
}
