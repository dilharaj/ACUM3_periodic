// Internal functions for BB noise calculation (TBL-LE)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void a_min(float* a, float* amin, int n)
{

for (int i=0; i<n; i++)
{
	if(a[i] < 0.204)
	{
		amin[i] = sqrt(67.552 - 886.788*pow(a[i],2)) - 8.219;
	} else
	{
		if(a[i] > 0.244)
		{
			amin[i] = -142.795*pow(a[i],3) + 103.656*pow(a[i],2) - 57.757*a[i] + 6.006;
		} else
		{
			amin[i] = -32.665*a[i] + 3.981;
		}
	}
}
}

void a_max(float* a, float* amax, int n)
{

for (int i=0; i<n; i++)
{
	if(a[i] < 0.13)
	{
		amax[i] = sqrt(67.552 - 886.788*pow(a[i],2)) - 8.219;
	} else
	{
		if(a[i] > 0.321)
		{
			amax[i] = -4.669*pow(a[i],3) + 3.491*pow(a[i],2) - 16.699*a[i] + 1.149;
		} else
		{
			amax[i] = -15.901*a[i] + 1.098;
		}
	}
}
}

void b_min(float* b, float* bmin, int n)
{

for (int i=0; i<n; i++)
{
	if(b[i] < 0.13)
	{
		bmin[i] = sqrt(16.888 - 886.788*pow(b[i],2)) - 4.109;
	} else
	{
		if(b[i] > 0.145)
		{
			bmin[i] = -817.810*pow(b[i],3) + 355.21*pow(b[i],2) - 135.024*b[i] + 10.619;
		} else
		{
			bmin[i] = -83.607*b[i] + 8.138;
		}
	}
}
}

void b_max(float* b, float* bmax, int n)
{

for (int i=0; i<n; i++)
{
	if(b[i] < 0.1)
	{
		bmax[i] = sqrt(16.888 - 886.788*pow(b[i],2)) - 4.109;
	} else
	{
		if(b[i] > 0.187)
		{
			bmax[i] = -80.541*pow(b[i],3) + 44.174*pow(b[i],2) - 39.381*b[i] + 2.344;
		} else
		{
			bmax[i] = -31.313*b[i] + 1.854;
		}
	}
}
}

float g5comp(float hdstar,float eta)
{

float m,k,mu,eta0,g5;


if(hdstar<=0.25) mu = 0.1211;
else if (hdstar<=0.62) mu = -0.2175*hdstar+0.1755;
else if(hdstar<1.15) mu = -0.0308*hdstar+0.0596;
else mu = 0.0242;

if(hdstar<=0.02) m=0.0;
else if(hdstar<0.5) m=68.724*hdstar-1.35;
else if(hdstar<=0.62) m=308.475*hdstar-121.23;
else if(hdstar<=1.15) m=224.811*hdstar-69.354;
else if(hdstar<1.2) m=1583.28*hdstar-1631.592;
else m=268.344;

if(m<0.0) m=0.0;

eta0 = -sqrt((m*m*(pow(mu,4)))/(6.25+m*m*mu*mu));

k = 2.5*sqrt(1.0-(pow(eta0/mu,2))) -2.5-m*eta0;

if(eta<=eta0) g5=m*eta+k;
else if(eta<=0.0) g5=2.5*sqrt(1.0-pow(eta/mu,2))-2.5;
else if(eta<=0.03616) g5=sqrt(1.5625-1194.99*pow(eta,2))-1.25;
else g5=-155.543*eta+4.375;

return g5;

}


float interp_2d(float* xyzcR,int* xsPtr,int isurf,int isect,float* sect_t,int* sPtr,int is1,int is2,float mach,float* mach_t,int* mPtr,int im1,int im2,float* re,int ptr2,float* re_t,int* rPtr,int ir1,int ir2,int* dPtr,int* nmach_t,int* nre_t,int* naoa_t,float* bbdata_t,int id,int i){

		// i = aoa index
		// id = fata index

	float ws1,wm1,wr1,F1,F2,F3,F4,cl;
	int q1,q2;
//printf("check\n");
	if(is1==is2) ws1=0.5;
	else ws1 = (xyzcR[xsPtr[isurf]+4*isect] - sect_t[sPtr[isurf]+is1])/(sect_t[sPtr[isurf]+is2] - sect_t[sPtr[isurf]+is1]);
	
    	if(im1==im2) wm1=0.5;
	else wm1 = (mach - mach_t[mPtr[isurf]+im1])/(mach_t[mPtr[isurf]+im2] - mach_t[mPtr[isurf]+im1]);	
	if(ir1==ir2) wr1=0.5;
	else wr1 = (re[ptr2] - re_t[rPtr[isurf]+ir1])/(re_t[rPtr[isurf]+ir2] - re_t[rPtr[isurf]+ir1]);	

	q1 = dPtr[isurf]+is1*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im1*nre_t[isurf]*naoa_t[isurf]*7+ir1*naoa_t[isurf]*7+i*7 + id;
	q2 = dPtr[isurf]+is2*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im1*nre_t[isurf]*naoa_t[isurf]*7+ir1*naoa_t[isurf]*7+i*7 + id;

//printf("%d %d %d %d %d %d %d %d %d %d \n",q1,q2,is1,is2,im1,im2,ir1,ir2,i,id);
	F1 = (1-ws1)*bbdata_t[q1] + ws1*bbdata_t[q2];

	q1 = dPtr[isurf]+is1*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im2*nre_t[isurf]*naoa_t[isurf]*7+ir1*naoa_t[isurf]*7+i*7 + id;
	q2 = dPtr[isurf]+is2*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im2*nre_t[isurf]*naoa_t[isurf]*7+ir1*naoa_t[isurf]*7+i*7 + id;
//printf("%d %d\n",q1,q2);
	F2 = (1-ws1)*bbdata_t[q1] + ws1*bbdata_t[q2];

	q1 = dPtr[isurf]+is1*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im1*nre_t[isurf]*naoa_t[isurf]*7+ir2*naoa_t[isurf]*7+i*7 + id;
	q2 = dPtr[isurf]+is2*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im1*nre_t[isurf]*naoa_t[isurf]*7+ir2*naoa_t[isurf]*7+i*7 + id;
//printf("%d %d\n",q1,q2);
	F3 = (1-ws1)*bbdata_t[q1] + ws1*bbdata_t[q2];

	q1 = dPtr[isurf]+is1*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im2*nre_t[isurf]*naoa_t[isurf]*7+ir2*naoa_t[isurf]*7+i*7 + id;
	q2 = dPtr[isurf]+is2*nmach_t[isurf]*nre_t[isurf]*naoa_t[isurf]*7+im2*nre_t[isurf]*naoa_t[isurf]*7+ir2*naoa_t[isurf]*7+i*7 + id;
//printf("%d %d\n",q1,q2);
	F4 = (1-ws1)*bbdata_t[q1] + ws1*bbdata_t[q2];

	cl = (1-wm1)*(1-wr1)*F1 + wm1*(1-wr1)*F2 + (1-wm1)*wr1*F3 + wm1*wr1*F4;


//printf("%f %f %f %f %f %f %f %f %f %f\n",ws1,wm1,wr1,F1,F2,F3,F4,cl,bbdata_t[q1],bbdata_t[q2]);
	return cl;
	

}



