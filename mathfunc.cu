#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
# define MAX(x,y) ((x) > (y) ? x: y)
# define MIN(x,y) ((x) < (y) ? x: y)



void redistribute(float* t,float* a,float* b,float ttime,int n){

   int ii,imin;
   float tmin;
   float* tmp = (float*)malloc(n*3*sizeof(float));


  
   for (int i=0;i<36;i++)

   for (int i=0;i<n;i++)
   {
      tmp[i*3+0] = fmod(t[i],ttime);
   }
   

   for (int i=0;i<n;i++)
   tmin = tmp[0];
   imin = 0;
   for (int i=1;i<n;i++){
      if (tmp[i*3+0] < tmin){
         tmin = tmp[i*3+0];
         imin = i;
      }
   }
   for (int i=0;i<n;i++){
      ii = (imin+i)%n;
      t[i] = tmp[ii*3+0];
      tmp[i*3+1] = a[ii];
      tmp[i*3+2] = b[ii];
   }

   for (int i=0;i<n;i++){   
      a[i] = tmp[i*3+1];
      b[i] = tmp[i*3+2];
   }
   

   free(tmp);
}

void interp_1d(float* radin,float* datin,float* radout,float* datout,int nrin,int nrout){

     int found,j;

     for (int r=0;r<nrout;r++){
        // Extrapolate if outside limits
        if (radout[r]<radin[0]){
           datout[r] = (radout[r]-radin[1])*(datin[1]-datin[0])/(radin[1]-radin[0])+datin[1];}
        else if (radout[r]>radin[nrin-1]){
           datout[r] = (radout[r]-radin[nrin-1])*(datin[nrin-1]-datin[nrin-2])/(radin[nrin-1]-radin[nrin-2])+datin[nrin-1];}
        else{ 
        // Otherwise interpolate finding the closest points
           found = 0;
           j = 1;
           while (found==0 && j<nrin){
              if ((radin[j]-radout[r])*(radin[j-1]-radout[r])<=0){
                 datout[r] = (datin[j]*(radout[r]-radin[j-1])+datin[j-1]*(radin[j]-radout[r]))/(radin[j]-radin[j-1]);
                 found = 1;
              }
              j++;
           }
        }
     }
}

void phase_shift(float* u, int N,float angle){
		
	float Xreal [N];
	float Ximag [N];
	float areal;
	float aimag;
	float breal; 
        float bimag; 
        float bbreal;
	float creal; 
        float cimag; 

	// Fourier Transform

	for (int j=0; j<N; j++){
		Xreal[j] = 0;
		Ximag[j] = 0;
	}
	for (int k=0; k<N; k++){
		for (int n=0; n<N; n++){
			 areal = u[n]*cos(-2*M_PI*k*n/N);
			 aimag = u[n]*sin(-2*M_PI*k*n/N);	
			Xreal[k] = Xreal[k] + areal;
			Ximag[k] = Ximag[k] + aimag;
		}
	}		

	// Phase shift
	
	for (int j=0; j<N; j++){
		breal = cos(j*angle*M_PI/180);	
		bimag = sin(j*angle*M_PI/180);
		bbreal = Xreal[j]*breal - Ximag[j]*bimag;
		Ximag[j] = Ximag[j]*breal + bimag*Xreal[j];
		Xreal[j] = bbreal;
	}
	
	// Inverse transform

	float x [N];
	for (int j=0; j<N; j++){
		x[j] = 0;
	}

	for (int n=0; n<N; n++){
		for (int k=0; k<N; k++){
			creal = cos(2*M_PI*k*n/N)/N;
			cimag = sin(2*M_PI*k*n/N)/N;
			x[n] = x[n] + Xreal[k]*creal-Ximag[k]*cimag;
		}
	}
	for (int j=0; j<N; j++){
		u[j]=x[j];
	}
}

void differentiate(float* u, float* ud, float dt, int n){
	
	ud[0]   = (u[1] - u[n-1])/(2.0*dt);      // central
	ud[n-1] = (u[0] - u[n-2])/(2.0*dt);    // central
        for(int i=1; i<n-1; i++)
        {
		ud[i] = (u[i+1] - u[i-1])/(2.0*dt);
	}
}

void differentiate_biased(float* u, float* ud, float dt, int n){
	
	ud[0]   = (u[1] - u[0])/(dt);      
	ud[n-1] = (u[n-1] - u[n-2])/(dt);  

	//ud[0]   = (-u[2] +4.0*u[1] -3.0*u[0])/(2.0*dt);      
	//ud[n-1] = (3.0*u[n-1] -4.0*u[n-2] + u[n-3])/(2.0*dt);  

        for(int i=1; i<n-1; i++)
        {
		ud[i] = (u[i+1] - u[i-1])/(2.0*dt);
	}
	ud[0] = 3*ud[1] - 3*ud[2] + ud[3];
	ud[n-1] = 3*ud[n-2] - 3*ud[n-3] + ud[n-4];



}





void make_periodic(float* fun, int n, int is, float *f_cubic)
{
	int i,xmax;	
	float k1,k2,s1,s2,aa,bb,cc,dd;
	
	
	float fun2[n+1];

	for(i=0;i<n;i++)
	{
	  fun2[i+1] = fun[i];
	}
	fun2[0] = fun2[1] - (fun2[2]-fun2[1]);

	// boundary condition
	k2=fun2[n];
	k1=fun2[is];

	s2=fun2[n]-fun2[n-1];
	s1=fun2[1]-fun2[0];
	
	//s1=0;
	//s2=0;

	// coefficient
	xmax = n;
	dd = k1;
	cc = s1-s2;
	bb = -2.0*(s1-s2)/xmax - 3.0*(k1-k2)/xmax/xmax;
	aa = -bb/xmax - (s1-s2)/xmax/xmax - (k1-k2)/xmax/xmax/xmax;

	for (i=0;i<n+1;i++)
	{
		f_cubic[i] = aa*i*i*i + bb*i*i + cc*i + dd;
		fun2[i] = fun2[i] - (f_cubic[i] - k2);
	}
	
	for(i=0;i<n;i++)
	{
		fun[i] = fun2[i+1];
	}
}

void medfilt1(float* func, int n)
{
   int i;      
   float aa,bb,cc,fmax,fmin;
   //float dd,ee;
   float *f_out = (float *) malloc(sizeof(float)*n);  

   for(i=0;i<n;i++)
   {
     
     //
     // 3pts
     //
     if(i==0)
     {
       aa = func[n-1];
       bb = func[i];
       cc = func[i+1];
     }
     else if(i==n-1)
     {
       aa = func[i-1];
       bb = func[i];
       cc = func[0];
     }
     else
     {
       aa = func[i-1];
       bb = func[i];
       cc = func[i+1];
     }
     fmax = MAX(aa,bb);
     fmax = MAX(fmax,cc);
     fmin = MIN(aa,bb);
     fmin = MIN(fmin,cc);
     
     if(aa<fmax && aa>fmin)
     {
       f_out[i] = aa; 
     }
     else if(bb<fmax && bb>fmin)
     {
       f_out[i] = bb; 
     }
     else
     {
       f_out[i] = cc;
     }



 //    // 
 //    //5 stencil
 //    // 
 //        if(i==0)
 //        {
 //           aa = func[n-2];
 //           bb = func[n-1];
 //           cc = func[i];
 //           dd = func[i+1];
 //           ee = func[i+2];
 //        }
 //        else if(i==1)
 //        {
 //           aa = func[n-1];
 //           bb = func[i-1];
 //           cc = func[i];
 //           dd = func[i+1];
 //           ee = func[i+2];
 //        }
 //        else if(i==n-1)
 //        {
 //           aa = func[i-2];
 //           bb = func[i-1];
 //           cc = func[i];
 //           dd = func[0];
 //           ee = func[1];
 //        }
 //        else if(i==n-2)
 //        {
 //           aa = func[i-2];
 //           bb = func[i-1];
 //           cc = func[i];
 //           dd = func[i+1];
 //           ee = func[0];
 //        }
 //        else  
 //        { 
 //           aa = func[i-2];
 //           bb = func[i-1];
 //           cc = func[i];
 //           dd = func[i+1];
 //           ee = func[i+2];
 //        }

 //        fmax = MAX(aa,bb);
 //        fmax = MAX(fmax,cc);
 //        fmax = MAX(fmax,dd);
 //        fmax = MAX(fmax,ee);
 //       		
 //        fmin = MIN(aa,bb);
 //        fmin = MIN(fmin,cc);
 //        fmin = MIN(fmin,dd);
 //        fmin = MIN(fmin,ee);


 //    	if(aa<fmax && aa>fmin)
 //    	{
 //    	  f_out[i] = aa; 
 //    	}
 //    	else if(bb<fmax && bb>fmin)
 //    	{
 //    	  f_out[i] = bb; 
 //    	}
 //    	else if(cc<fmax && cc>fmin)
 //    	{
 //    	  f_out[i] = cc; 
 //    	}
 //    	else if(dd<fmax && dd>fmin)
 //    	{
 //    	  f_out[i] = dd; 
 //    	}
 //    	else
 //    	{
 //    	  f_out[i] = ee;
 //    	}

   }
   
   // copy to original
   for(i=0;i<n;i++) func[i] = f_out[i];

   free(f_out); 
}





//
// pre Gaussisn filter
//
void pre_gaussfilt1(float* gauss, int n_win)
{
   int i,m;      
   float std,sumg,xx,aa;

   std = (n_win-1)/6.0;
   aa  = 1.0/(sqrt(2.0*M_PI)*std);

   m=(n_win-1)/2;
   sumg=0.0;
   for (i=0;i<n_win;i++)
   {
      xx = i - m;
      gauss[i] = aa*exp(-0.5*(xx*xx)/std/std);
      sumg += gauss[i];

   }

   for(i=0;i<n_win;i++) gauss[i] = gauss[i]/sumg;

}

//
// Gaussisn filter
//
void gaussfilt1(float* func, float* gauss, int n, int n_win)
{
   int i,j,k,m;      
   float *f_out = (float *) malloc(sizeof(float)*n);  
   float sample[n_win];
   //float sum;

   for(i=0;i<n;i++)
   {
      // set samples
      m=(n_win-1)/2;
      for (j=0;j<n_win;j++)
      {
         k = i + j - m;
         
	 if(k<0)
         {
            k = n + k;
            sample[j] = func[k];
         }
         else if(k>n-1)
         {
            k = k - n;
            sample[j] = func[k];
         }
	 else
         {
            sample[j] = func[k]; 
         }
      }

      // compute weighted average
      f_out[i]=0.0;
      //sum=0.0;
      for (j=0;j<n_win;j++)
      {
         f_out[i] += gauss[j]*sample[j];
	 //sum += gauss[j];
      }

   }

   // copy to original
   for(i=0;i<n;i++) func[i] = f_out[i];
   
   free(f_out); 
}
