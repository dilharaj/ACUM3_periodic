#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "bb_utils.cuh"
#include "bb_core.cuh"
# define MAX(x,y) ((x) > (y) ? x: y)
# define MIN(x,y) ((x) < (y) ? x: y)
# define FLT_EPSLION 1.19209290e-07F // float

void findBBnoise(int nSurf_BB,int nObs,int nTime,float dTau,float* span,float* omega,int* ccw,float* psi_offset,int* trip,float* rotate,float* translate,int* xsPtr,int* dsPtr,int* nsect,float* xyzcR,float* bbdata,int* nsect_t,int* nmach_t,int* naoa_t,int* nre_t,int* aPtr,int* mPtr,int* rPtr,int* sPtr,int* dPtr,float* sect_t,float* mach_t,float* aoa_t,float* re_t,float* bbdata_t,int t_del,float* vInf,float* rObs, float* tObs,int nf2,float* ff,float a0,float* BBspl,float* BBoaspl)
{



// Preprocessing

float x0[3],x1[3],xt0[3],xt1[3],thx,thy,thz,psi,c,tr[3],v[3],v1[3],err,err_max,cn,ct,clsec;
float trans[3][3],bl_temp[8],mach,lsect,rp[3],aMv,bb,cc,tau,aa,rpnorm,rp2[3],r_e,phi;
int ptr,ptr2,ptr3,ptrm,ptrp,im1,im2,ir1,ir2,is1,is2,it,it_max,ch1,ch2,ch3,a1,a2,a3,iaoa;
float nu = 0.0000181206/1.22500;
float l1,theta,wa1,alpha_eff,cl_alpha0,alpha_cl0,alpha_eff1,delta_0,deltas_0;

//float CT1[nSurf_BB], CT2[nSurf_BB];
//
//for(int i=0; i<nSurf_BB; i++){
//	CT1[i] = 0.0;
//	CT2[i] = 0.0;
//}


//float aoain[] = {-3.7966981914327325, 0.5732962011181550, 3.3347710910654129, 5.5184459670734469, 7.3646080488360974, 7.3937133211004369, 6.0602296177163968, 4.4214833814613881, 2.6454674941799068, 0.9564961989277790, 0.2268483126317935,-0.0193160482059958,-0.2026326418192593,-0.3651270826102695,-0.7966310473752259,-1.4924294151394493};



// bb output saving arrays

float* output_tmp = (float*)malloc(nObs*nTime*nf2*7*sizeof(float));

for(int i=0; i<nObs*nTime*nf2*7; i++) output_tmp[i] = -100.0;

for (int isurf=0; isurf<nSurf_BB; isurf++)
{

    printf("Surf: %d\n",isurf+1);
    printf("   preprocessing...\n");

    float* xyz = (float*)malloc(nTime*nsect[isurf]*3*sizeof(float));
    float* xyzt = (float*)malloc(nTime*nsect[isurf]*3*sizeof(float));
    float* u_ms = (float*)malloc(nTime*nsect[isurf]*sizeof(float));
    float* re = (float*)malloc(nTime*nsect[isurf]*sizeof(float));
    float* alpha = (float*)malloc(nTime*nsect[isurf]*sizeof(float));
    float* bl = (float*)malloc(nTime*nsect[isurf]*4*sizeof(float));
    float* cl_alpha = (float*)malloc(nTime*nsect[isurf]*sizeof(float));
    
    int* km = (int*)malloc(nsect[isurf]*sizeof(int));
    int* kp = (int*)malloc(nsect[isurf]*sizeof(int));


    for (int isect=0; isect<nsect[isurf]; isect++)
    {
    	if (isect==nsect[isurf]-1)
	{
	    km[isect] = isect-1;
	    kp[isect] = isect;
	} else 
	{
	    km[isect] = isect;
	    kp[isect] = isect+1;
	}
    }    
    thx = rotate[3*isurf]*M_PI/180;
    thy = rotate[3*isurf+1]*M_PI/180;
    thz = rotate[3*isurf+2]*M_PI/180;

	
    for (int i=0; i<3; i++) tr[i] = translate[3*isurf+i];

    for (int isect=0; isect<nsect[isurf]; isect++)
    {
	c = xyzcR[xsPtr[isurf] + 4*isect + 3] * span[isurf];
        for (int itime=0; itime<nTime; itime ++)	    
        {
	    // position vectors

		// qtr chrd 
	    for (int i=0; i<3; i++) x0[i] = xyzcR[xsPtr[isurf] + 4*isect + i] * span[isurf];
	    for (int i=1; i<3; i++) x0[i] = xyzcR[xsPtr[isurf] + 4*(int)round(nsect[isurf]*0.75) + i] * span[isurf];



	    if( ccw[isurf]==0) x0[1] = -x0[1]; // y

	    psi = (dTau*itime*omega[isurf] + psi_offset[isurf]*M_PI/180);
	    if(ccw[isurf] == 0) psi = -psi;

	    x1[2] = x0[2];
	    x1[1] = x0[0]*sin(psi) + x0[1]*cos(psi);
	    x1[0] = x0[0]*cos(psi) - x0[1]*sin(psi);
	
		// trailing edge
	    xt0[0] = x0[0];
	    
	    theta = bbdata[dsPtr[isurf] + itime*nsect[isurf]*7 + isect*7 + 6]*M_PI/180;

	    xt0[2] = x0[2] - c*0.75*sin(theta);
	
	    if(ccw[isurf] == 1)
	    {
	    	xt0[1] = x0[1] - c*0.75*cos(theta);
	    } else {
	    	xt0[1] = x0[1] + c*0.75*cos(theta);
	    }

	    xt1[2] = xt0[2];
	    xt1[1] = xt0[0]*sin(psi) + xt0[1]*cos(psi);
	    xt1[0] = xt0[0]*cos(psi) - xt0[1]*sin(psi);
	

	         // axis rotation

	  	// around x axis
	    xt0[0] = xt1[0];
	    xt0[2] = xt1[1]*sin(thx) + xt1[2]*cos(thx);
	    xt0[1] = xt1[1]*cos(thx) - xt1[2]*sin(thx);
	
	    x0[0] = x1[0];
	    x0[2] = x1[1]*sin(thx) + x1[2]*cos(thx);
	    x0[1] = x1[1]*cos(thx) - x1[2]*sin(thx);

		// around y axis

	    xt1[1] =  xt0[1];
	    xt1[2] = -xt0[0]*sin(thy) + xt0[2]*cos(thy);
	    xt1[0] =  xt0[0]*cos(thy) + xt0[2]*sin(thy);
	
	    x1[1] =  x0[1];
	    x1[2] = -x0[0]*sin(thy) + x0[2]*cos(thy);
	    x1[0] =  x0[0]*cos(thy) + x0[2]*sin(thy);

//	if(itime==0) printf("isurf:%d isect:%d x:%f y:%f z:%f xt:%f yt:%f zt:%f\n",isurf,isect,x1[0],x1[1],x1[ptr+2],xyzt[ptr+0],xyzt[ptr+1],xyzt[ptr+2],thx,thy,thz,tr[0],tr[1],tr[2]);
	
		// around z axis

	    xt0[2] = xt1[2];
	    xt0[1] = xt1[0]*sin(thz) + xt1[1]*cos(thz);
	    xt0[0] = xt1[0]*cos(thz) - xt1[1]*sin(thz);
	
	    x0[2] = x1[2];
	    x0[1] = x1[0]*sin(thz) + x1[1]*cos(thz);
	    x0[0] = x1[0]*cos(thz) - x1[1]*sin(thz);
	

		// translation

	    ptr = isect*nTime*3 + itime*3; 
	    for (int i=0; i<3; i++)
	    {
	    	xyz[ptr + i] = x0[i] + tr[i];	
	    	xyzt[ptr + i] = xt0[i] + tr[i];
	    }

	//if(itime==0) printf("isurf:%d isect:%d x:%f y:%f z:%f xt:%f yt:%f zt:%f rx:%f ry:%f rz:%f tx:%f ty:%f tz:%f\n",isurf,isect,xyz[ptr+0],xyz[ptr+1],xyz[ptr+2],xyzt[ptr+0],xyzt[ptr+1],xyzt[ptr+2],thx,thy,thz,tr[0],tr[1],tr[2]);
	}
	    // flow velocity wrt surface + Mach and Re
	for (int itime=0; itime<nTime; itime ++)
	{
	    ptr = isect*nTime*3 + itime*3; 
	    ptr2 = isect*nTime + itime; 
	    if (itime==0)
	    {
   		ptrm = isect*nTime*3 + (nTime-1)*3; 
		ptrp = isect*nTime*3 + (itime+1)*3; 
	    } else if (itime == nTime-1)
	    {
		ptrm = isect*nTime*3 + (itime-1)*3; 
		ptrp = isect*nTime*3; 
	    }else
	    {
	    	ptrm = isect*nTime*3 + (itime-1)*3;
		ptrp = isect*nTime*3 + (itime+1)*3;
	    }

	    
	    for (int i=0; i<3; i++)
	    {
	        v[i] = vInf[i] -(xyz[ptrp+i] - xyz[ptrm+i])/(2*dTau); //flow velocity wrt surface
	    }
	
            // Mach and Re
	    

	    for (int i=0; i<3; i++) v1[i] = xyzt[ptr+i] - xyz[ptr+i];

	    l1 = sqrt(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2));

	    for (int i=0; i<3; i++) v1[i] = v1[i]/l1;

	    u_ms[ptr2] = (v[0]*v1[0] + v[1]*v1[1] + v[2]*v1[2]); // chordwise U

////////////////// HARDCODED////////////
//	   u_ms[ptr2] = 126*1.0*xyzcR[xsPtr[isurf]+4*isect];
//	if(isurf==5) u_ms[ptr2] = 0.197*340;
////////////////////////////////////////////////////
	    mach = fabs(u_ms[ptr2]/a0); 
	    re[ptr2] = fabs(u_ms[ptr2]*c/nu);

	    im2 = 0;
	    for (int j=0; j<nmach_t[isurf]; j++)
	    {
	        if (mach < mach_t[mPtr[isurf] + j]) break;
	        im2 += 1;
	    }
	
	    ir2 = 0;   
	    for (int j=0; j<nre_t[isurf]; j++)
	    {
//if(isurf==0&&itime==0&&isect>6&&isect<9) printf("%f %f %d\n",re[ptr2],re_t[rPtr[isurf]+j],rPtr[isurf]);

	        if (re[ptr2] < re_t[rPtr[isurf] + j]) break;
	        ir2 += 1;
	    }
	    is2 = 0;
	    for (int j=0; j<nsect_t[isurf]; j++)
	    {
	        if (xyzcR[xsPtr[isurf]+4*isect] < sect_t[sPtr[isurf] + j]) break;
	        is2 += 1;
	    }
	    im1 = MAX(0,im2-1);
	    ir1 = MAX(0,ir2-1);
	    is1 = MAX(0,is2-1);
	    im2 = MIN(im2,nmach_t[isurf]-1);
	    ir2 = MIN(ir2,nre_t[isurf]-1);
	    is2 = MIN(is2,nsect_t[isurf]-1);

//if(isurf==0&&itime==0) printf("%d %d %d %d %d %d\n",im1,ir1,is1,im2,ir2,is2);

	    // cl curve 
	    float cl[naoa_t[isurf]]; 
	    for (int i=0; i<naoa_t[isurf]; i++)
	    {
	
//printf("ch\n");
		
//if(isect==5&&itime==10)	
		cl[i] = interp_2d(xyzcR,xsPtr,isurf,isect,sect_t,sPtr,is1,is2,mach,mach_t,mPtr,im1,im2,re,ptr2,re_t,rPtr,ir1,ir2,dPtr,nmach_t,nre_t,naoa_t,bbdata_t,0,i);

//if(isurf==0&&itime==0&&isect>6&&isect<9) printf("%f %f\n",aoa_t[aPtr[isurf]+i],cl[i]);
	
	    } //naoa_t

	    // cl_alpha
	
            ch1 = 0;
	    ch2 = 0;
	    ch3 = 0;
	    for (int i=0; i<naoa_t[isurf]; i++)
	    {
		if (aoa_t[aPtr[isurf]+i] >= 0.0 && ch1==0)
		{
			a1 = i;
			ch1 = 1;
		}
		if (aoa_t[aPtr[isurf]+i] >= 6.0 && ch2==0)
		{
			a2 = i;
			ch2 = 1;
		}
		if (cl[i]>0 && ch3==0)
		{
			a3 = i;
			ch3 = 1;
		}
		if(ch1==1 && ch2==1 && ch3==1) break;
	    }

	    // cl_alpha (per deg)
	    cl_alpha[ptr2] = (cl[a2] - cl[a1])/(aoa_t[aPtr[isurf]+a2] - aoa_t[aPtr[isurf]+a1]);

	    // cl @ alpha=0
	    wa1 = (0-aoa_t[aPtr[isurf]+a1-1])/(aoa_t[aPtr[isurf]+a1] - aoa_t[aPtr[isurf]+a1-1]);
	    cl_alpha0 = cl[a1-1]*(1-wa1) + cl[a1]*wa1;

	    // alpha @ cl=0 (deg)
	    wa1 = (0-cl[a3-1])/(cl[a3]-cl[a3-1]);
	    alpha_cl0 = aoa_t[aPtr[isurf]+a3-1]*(1-wa1) + aoa_t[aPtr[isurf]+a3]*wa1;

//aoa_t[aPtr[isurf]+a2],aoa_t[aPtr[isurf]+a1],wa1,a1);
	    // aoa_eff (deg)

	    it = 0;
	    it_max = 50;
	    err = 1;
	    err_max = 1e-3;
	    alpha_eff = 0.0;

	    cn = bbdata[dsPtr[isurf] + itime*nsect[isurf]*7 + isect*7 + 4];/// pow(xyzcR[xsPtr[isurf] + 4*isect],2);
	    ct = bbdata[dsPtr[isurf] + itime*nsect[isurf]*7 + isect*7 + 5];/// pow(xyzcR[xsPtr[isurf] + 4*isect],2);

//if(itime==0) printf("%d %f %f %f %f %f\n",isect,cl[a2],cl[a1],alpha_cl0,cl_alpha[ptr2],cl_alpha0);

	    while (err>err_max)
	    {
	    		
		clsec = cn*cos(alpha_eff*M_PI/180) - ct*sin(alpha_eff*M_PI/180);
		//cdsec = cn*sin(alpha_eff*M_PI/180) + ct*cos(alpha_eff*M_PI/180);

		alpha_eff1 = (clsec-cl_alpha0)/cl_alpha[ptr2];

		err = abs((alpha_eff1-alpha_eff));

//if(isurf==0&&itime==0&&isect==0) printf("%d %f %f %f %f %f %f %f\n",it,err,alpha_eff,alpha_eff1,clsec,cl_alpha[ptr2],cn,ct);
		alpha_eff = alpha_eff1;

		it += 1;

		if (it>it_max) break;



	    }
//printf("%d %d %f %f %f %f %f %f %f %f\n",isect,it,alpha_eff,alpha_cl0,mach,re[ptr2],cn,ct,clsec,cdsec);
	 
//if(isurf==0&&itime==0) printf("isect:%d alpha_eff:%f alpha_cl0:%f alpha:%f re:%f mach:%f\n",isect,alpha_eff,alpha_cl0,alpha_eff-alpha_cl0,re[ptr2],mach); 
//if(isect==nsect[isurf]-2&&itime==0) printf("%f %f %f %f %f %d %f %f %d\n",alpha_eff,cn,ct,clsec,cl_alpha,it,alpha_cl0,cl_alpha0,dsPtr[isurf]);
	    // total effective alpha for BB calculation (in deg)


//	    theta = bbdata[dsPtr[isurf] + itime*nsect[isurf]*7 + isect*7 + 6]*M_PI/180;
//	    if (isect==0) lsect = (xyzcR[xsPtr[isurf]+4*(isect+1)] - xyzcR[xsPtr[isurf]+4*isect])/2.0;
//	    else if (isect==nsect[isurf]-1) lsect = (xyzcR[xsPtr[isurf]+4*isect] - xyzcR[xsPtr[isurf]+4*(isect-1)])/2.0;
//	    else lsect = (xyzcR[xsPtr[isurf]+4*(isect+1)] - xyzcR[xsPtr[isurf]+4*(isect-1)])/2.0;
//
//if(isurf==0&&itime==0) printf("%d %f %f %f %f %f %f\n",isect,alpha_eff,cn,ct,lsect,clsec,cdsec);
//
//	    CT1[isurf] += 1.0/nTime*(cn*cos(theta)-ct*sin(theta))*lsect*span[isurf]*c*span[isurf]/0.3973;
//	    CT2[isurf] += 1.0/nTime*(clsec*cos(theta-alpha_eff*M_PI/180) - cdsec*sin(theta-alpha_eff*M_PI/180))*lsect*span[isurf]*c*span[isurf]/0.3973;


	    alpha[ptr2] = alpha_eff - alpha_cl0;
	    alpha[ptr2] = fabs(alpha[ptr2]);

//// HARDCODED //////////////////////////////


	   //alpha[ptr2] = fabs((aoain[isect]));
	   //alpha_cl0 = 0.0;

//////////////////////////////////////////


//printf("tdel=%d\n",t_del);
	    // BL thickness
	    if (t_del == 0) // 3D CFD BL data
	    {
	        for (int i=0; i<4; i++) bl[isect*nTime*4+itime*4+i] = bbdata[dsPtr[isurf] + itime*nsect[isurf]*7 + isect*7 + i];
	    } else if (t_del == 1)  // 2D CFD BL data
	    {
		// interpolation for BL thickness
	    
		for (int i=0; i<naoa_t[isurf]; i++)
		{
   		    if (aoa_t[aPtr[isurf]+i] >= alpha_eff)
		    {
			a1 = i;
			break;
		    }
		}
	
		wa1 = (alpha_eff - aoa_t[aPtr[isurf]+a1-1])/(aoa_t[aPtr[isurf]+a1] - aoa_t[aPtr[isurf]+a1-1]);
		
		for (int jj=0; jj<2; jj++)
		{

		    im2 = 0;
		    for (int j=0; j<nmach_t[isurf]; j++)
		    {
		        if (mach < mach_t[mPtr[isurf] + j]) break;
		        im2 += 1;
		    }
		    
		    ir2 = 0;   
		    for (int j=0; j<nre_t[isurf]; j++)
		    {
		        if (re[ptr2] < re_t[rPtr[isurf] + j]) break;
		        ir2 += 1;
		    }
		    is2 = 0;
		    for (int j=0; j<nsect_t[isurf]; j++)
		    {
		        if (xyzcR[xsPtr[isurf]+4*isect] < sect_t[sPtr[isurf] + j]) break;
		        is2 += 1;
		    }
		    im1 = MAX(0,im2-1);
		    ir1 = MAX(0,ir2-1);
		    is1 = MAX(0,is2-1);
		    im2 = MIN(im2,nmach_t[isurf]-1);
		    ir2 = MIN(ir2,nre_t[isurf]-1);
		    is2 = MIN(is2,nsect_t[isurf]-1);
	
		    iaoa = jj + a1-1;

		    for (int i=0; i<4; i++)
		    {
		        bl_temp[jj*4+i] = interp_2d(xyzcR,xsPtr,isurf,isect,sect_t,sPtr,is1,is2,mach,mach_t,mPtr,im1,im2,re,ptr2,re_t,rPtr,ir1,ir2,dPtr,nmach_t,nre_t,naoa_t,bbdata_t,i+3,iaoa);	    
		    }
		}
		for (int i=0; i<4; i++)
		   {
		    bl[isect*nTime*4+itime*4+i] = c*(bl_temp[i]*(1-wa1) + bl_temp[4+i]*wa1);
//if(itime==0&&isect==nsect[isurf]-2) printf("bl=%f", bl[isect*nTime*4+itime*4+i]);
		}

	    } else // Empirical equations
	    { 

		delta_0 = c*pow(10.0,(1.6569-0.9045*log10(re[ptr2])+0.0596*pow(log10(re[ptr2]),2.0)));
		deltas_0 = c*pow(10.0,(3.0187-1.5397*log10(re[ptr2])+0.1059*pow(log10(re[ptr2]),2.0)));

		bl[isect*nTime*4+itime*4+1] = delta_0*pow(10.0,(-0.04175*alpha[ptr2]+0.00106*pow(alpha[ptr2],2.0)));
		bl[isect*nTime*4+itime*4+3]= deltas_0*pow(10.0,(-0.0432*alpha[ptr2]+0.00113*pow(alpha[ptr2],2.0)));


		if((alpha[ptr2])<=7.5) bl[isect*nTime*4+itime*4] = delta_0*pow(10.0,(0.03114*alpha[ptr2]));
		else if((alpha[ptr2])>12.5) bl[isect*nTime*4+itime*4] = delta_0*12.0*pow(10.0,0.0258*alpha[ptr2]);
		else bl[isect*nTime*4+itime*4] = delta_0*0.0303*pow(10.0,0.2336*alpha[ptr2]);



		if((alpha[ptr2])<=7.5) bl[isect*nTime*4+itime*4+2] = deltas_0*pow(10.0,(0.0679*alpha[ptr2]));
		else if((alpha[ptr2])>12.5) bl[isect*nTime*4+itime*4+2] = deltas_0*52.42*pow(10.0,0.0258*alpha[ptr2]);
		else bl[isect*nTime*4+itime*4+2] = deltas_0*0.0162*pow(10.0,0.3066*alpha[ptr2]);

//if(itime==0) printf("BL: %f %f %f %f\n",bl[isect*nTime*4+itime*4],bl[isect*nTime*4+itime*4+1],bl[isect*nTime*4+itime*4+2],bl[isect*nTime*4+itime*4+3]);
//if(itime==0) printf("%d %f %f %f %f\n",isect+1,u_ms[ptr2],re[ptr2],alpha[ptr2],c);


	    }

        }
        
    }


	// core calculation
    
    printf("   calculating noise...\n");
    
    for (int iobs=0; iobs<nObs; iobs++)    
    {
	printf("observer: %d\n",iobs);
	for (int isect=0; isect<nsect[isurf]; isect++)
	{
	//printf("section: %d\n",isect);

	// calculate section length taking data at cell edges
	    if (isect==0) lsect = (xyzcR[xsPtr[isurf]+4*(isect+1)] - xyzcR[xsPtr[isurf]+4*isect])/2.0*span[isurf];
	    else if (isect==nsect[isurf]-1) lsect = (xyzcR[xsPtr[isurf]+4*isect] - xyzcR[xsPtr[isurf]+4*(isect-1)])/2.0*span[isurf];
	    else lsect = (xyzcR[xsPtr[isurf]+4*(isect+1)] - xyzcR[xsPtr[isurf]+4*(isect-1)])/2.0*span[isurf];


	    c = xyzcR[xsPtr[isurf] + 4*isect + 3] * span[isurf];

//if(isect==nsect[isurf]-2)   printf("chord %f %f %f %d %d %d\n",c,xyzcR[xsPtr[isurf] + 4*isect + 3],span[isurf],xsPtr[isurf],isect,isurf);
	    for (int itime=0; itime<nTime; itime ++)	    
	    {
	 	   	
		ptr = isect*nTime*3 + itime*3; 
		ptr2 = isect*nTime + itime; 
		ptr3 = isect*nTime*4 + itime*4; 


		// chordwise vector
		float tempr = 0.0;
		for (int i=0; i<3; i++)
		{
		    trans[0][i] =  xyzt[ptr+i] - xyz[ptr+i];
		    tempr += pow(trans[0][i],2);
		}
		for (int i=0; i<3; i++) trans[0][i] = trans[0][i]/sqrt(tempr);

		// spanwise vector
		tempr = 0.0;
		for (int i=0; i<3; i++) 
		{
		    trans[1][i] = xyz[kp[isect]*nTime*3 + itime*3+i] - xyz[km[isect]*nTime*3 + itime*3+i];

//if(itime==480) printf("p:%f m:%f\n",xyzt[kp[isect]*nTime*3 + itime*3+i],xyzt[km[isect]*nTime*3 + itime*3+i]);
		    tempr += pow(trans[1][i],2);
		}
		for (int i=0; i<3; i++) trans[1][i] = trans[1][i]/sqrt(tempr);
	
		// surface normal vector
		if (ccw[isurf] == 1)
		{
		    trans[2][0] = trans[0][1]*trans[1][2] - trans[0][2]*trans[1][1];
		    trans[2][1] = trans[0][2]*trans[1][0] - trans[0][0]*trans[1][2];
		    trans[2][2] = trans[0][0]*trans[1][1] - trans[0][1]*trans[1][0];

		    tempr = sqrt(pow(trans[2][0],2) + pow(trans[2][1],2) + pow(trans[2][2],2));

	            for (int i=0; i<3; i++) trans[2][i] = trans[2][i]/tempr;

		} else {
		    
		    trans[2][0] = trans[1][1]*trans[0][2] - trans[1][2]*trans[0][1];
		    trans[2][1] = trans[1][2]*trans[0][0] - trans[1][0]*trans[0][2];
		    trans[2][2] = trans[1][0]*trans[0][1] - trans[1][1]*trans[0][0];
		    
		    tempr = sqrt(pow(trans[2][0],2) + pow(trans[2][1],2) + pow(trans[2][2],2));

		    for (int i=0; i<3; i++) trans[2][i] = trans[2][i]/sqrt(tempr);
		}
	
		for (int i=0; i<3; i++) rp[i] = rObs[iobs*3+i] - xyzt[ptr+i];

		aMv = pow(a0,2) - pow(vInf[0],2) - pow(vInf[1],2) - pow(vInf[2],2);
		bb = rp[0]*vInf[0] + rp[1]*vInf[1] + rp[2]*vInf[2];
		cc = rp[0]*rp[0] + rp[1]*rp[1] + rp[2]*rp[2];
		tau = (-bb+sqrt(pow(bb,2)+aa*cc))/aMv;
		for (int i=0; i<3; i++) rp[i] = rObs[iobs*3+i] - xyzt[ptr+i] - vInf[i]*tau;

		rpnorm = sqrt(pow(rp[0],2) + pow(rp[1],2) + pow(rp[2],2));

		for (int i=0; i<3; i++) rp2[i] = rp[0]*trans[i][0]+rp[1]*trans[i][1]+rp[2]*trans[i][2];

		r_e = sqrt(pow(rp2[1],2) + pow(rp2[2],2));

		theta = acos(abs(rp2[0]/rpnorm));
		phi = acos(abs(rp2[1])/r_e);

float test1,test2,test3;
test1=0;
test2=0;
test3=0;
for(int i=0;i<3;i++){
	test1 += trans[0][i]*trans[1][i];
	test2 += trans[0][i]*trans[2][i];
	test3 += trans[1][i]*trans[2][i];
}




		core(rpnorm,theta,phi,c,lsect,fabs(alpha[ptr2]),u_ms[ptr2],nu,a0,trip[isurf],bl[ptr3],bl[ptr3+1],bl[ptr3+2],bl[ptr3+3],isect,itime,iobs,output_tmp,nTime,nf2,nsect[isurf],isurf,ff,cl_alpha[ptr2],span[isurf]);	

//if(isect==5) printf("%d %e\n",itime,output_tmp[iobs*nTime*nf2*7 + itime*nf2*7 + 287]);

	    } // itime

	// Redistribution is avoided assuming periodicity

	

	 } // isect

    } //iobs

    free(xyz);free(xyzt);free(u_ms);free(re);free(alpha);free(km);free(kp);


} // isurf

printf("   aggregating into arrays...\n");
for (int iobs=0; iobs<nObs; iobs++)
{
   for (int i=0; i<nf2; i++)
   {
   	for (int j=0; j<7; j++)
	{
	    for (int itime=0; itime<nTime; itime++)
	    {
	    	BBspl[iobs*nf2*7+i*7+j] = 10.0*log10(pow(10.0,(BBspl[iobs*nf2*7+i*7+j]/10.0)) + pow(10.0,(output_tmp[iobs*nTime*nf2*7 + itime*nf2*7 + i*7 + j]/10.0)));

 
	    }
	    BBspl[iobs*nf2*7+i*7+j] = 10.0*log10(pow(10.0,(BBspl[iobs*nf2*7+i*7+j]/10.0))/(float)nTime);
	}
//printf("%f %f\n",BBoaspl[iobs],BBspl[iobs*nf2*7+i*7]);
	BBoaspl[iobs] = 10.0*log10(pow(10.0,(BBoaspl[iobs]/10.0)) + pow(10.0,(BBspl[iobs*nf2*7+i*7]/10.0)));
    }

}
free(output_tmp); 

}




void core(float rpnorm,float theta,float phi,float c,float lsect,float alpha,float u, float nu, float a0,int trip, float del_s,float del_p,float dis_s,float dis_p,int isect,int itime,int iobs,float* output_tmp,int nTime,int nf2,int nsect,int isurf,float* ff, float cl_alpha, float R)
{

//printf("\n");
float mach,r_c,m_c,d_h,d_l;
int optr;


if(del_s<0.0) del_s = 1.0e-10;
if(del_p<0.0) del_p = 1.0e-10;
if(dis_p<0.0) dis_p = 1.0e-10;
if(dis_s<0.0) dis_s = 1.0e-10;


u = abs(u);

optr = iobs*nTime*nf2*7 + itime*nf2*7;
mach = u/a0;

r_c = u*c/nu;

//if(itime==0&&isurf==0) printf("isect:%d r_c: %e\n",isect,r_c);  

m_c = 0.6*u/a0;

if(theta==0) theta=1e-10;
if(phi==0) phi=1e-10;


d_h = (2.0*pow(sin(theta/2.0),2) * pow(sin(phi),2)) / ((1.0+mach*cos(theta))*pow(1.0+(mach-m_c)*cos(theta),2));



d_l = (pow(sin(theta),2) * pow(sin(phi),2)) / (pow(1.0+mach*cos(theta),4));

//if(isurf==8&&itime==0) printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",isect,rpnorm,theta,phi,c,lsect,alpha,u,del_s,del_p,dis_s,dis_p,r_c,m_c,d_h,d_l);

//if(isect==18&&itime==480) printf("%f %f %f %f %f %f\n",d_h,d_l,theta,phi,mach,m_c);

// TBL-TE noise
float pmin[nf2],pmax[nf2],smin[nf2],smax[nf2],aas[nf2],aap[nf2],b[nf2],bmin[nf2],bmax[nf2],bb[nf2],spl_a[nf2],spl_p[nf2],spl_s[nf2],aminb[nf2],amaxb[nf2],anew2[nf2];

float r_deltas_p,st_1,st_2,st_s[nf2],st_p[nf2],e[nf2],as[nf2],ap[nf2],anew[nf2],a_0[1],b_0[1],min_0[1],max_0[1],a_r,b_r,k_1,deltak_1,gamm,gamm_0,beta,beta_0,k_2,r_c_new,stlprim,stpkprm,r_c0,d,g1[nf2],g2,g3,scal,spl_lam,spl_blunt,spl_tip;

r_deltas_p = u*dis_p/nu;

for (int i=0; i<nf2; i++)
{
    st_p[i] = ff[i]*dis_p/u;
    st_s[i] = ff[i]*dis_s/u;
}

st_1 = 0.02*(pow(mach,-0.6));

if (alpha < 1.33) 
{
	st_2 = st_1;
} else
{
	if (alpha > 12.5)
	{
		st_2 = st_1*4.72;
	} else
	{
		st_2 = st_1* pow(10.0,(0.0054*pow(alpha-1.33,2.0)));
	}
}

for (int i=0; i<nf2; i++)
{
    	as[i] = abs(log10(st_s[i]/st_1));
    	ap[i] = abs(log10(st_p[i]/st_1));
    	anew[i] = abs(log10(st_s[i]/st_2));
}

if( r_c < 9.52*pow(10,4))
{
	a_0[0] = 0.57;
}else
{
	if( r_c > 8.57*pow(10,5)) a_0[0] = 1.13;
	else a_0[0] = (-9.57*pow(10,-13))*(pow((r_c-(8.57*pow(10,5))),2)) + 1.13;
}

// compute function
min_0[0] = 0.0;
max_0[0] = 0.0;
a_min(a_0,min_0,1);
a_max(a_0,max_0,1);
a_min(as,smin,nf2);
a_max(as,smax,nf2);
a_min(ap,pmin,nf2);
a_max(ap,pmax,nf2);

a_r = (-20.0-min_0[0])/(max_0[0]-min_0[0]);

for (int i=0; i<nf2; i++)
{
	aas[i] = smin[i] + a_r*(smax[i]-smin[i]);
	aap[i] = pmin[i] + a_r*(pmax[i]-pmin[i]);

	b[i] = abs(log10(st_s[i]/st_2));
//if(i==10&&itime==0) printf("%f %f %f %f %f %f %f %f\n",aas[i],smin[i],a_r,smax[i],smin[i],st_s[i],st_1,dis_s); 
}


// determining b_0 based on Re

if (r_c < 9.52*pow(10,4))
{
	b_0[0] = 0.3;
} else
{
	if(r_c > 8.57*pow(10,5))
	{
		b_0[0] = 0.56;
	} else
	{
		b_0[0] = (-4.48*pow(10,-13))*(pow((r_c-(8.57*pow(10,5))),2)) + 0.56;
	}
}
// compute function

b_min(b_0,min_0,1);
b_max(b_0,max_0,1);
b_min(b,bmin,nf2);
b_max(b,bmax,nf2);

b_r = (-20.0-min_0[0])/(max_0[0]-min_0[0]);
for (int i=0; i<nf2; i++) bb[i] = bmin[i] + b_r*(bmax[i]-bmin[i]);

if(r_c<2.47*pow(10,5))
{
	k_1 = -4.31*log10(r_c) + 156.3;
} else
{
	if(r_c>8.0*pow(10,5))
	{
		k_1 = 128.5;
	} else
	{
		k_1 = -9.0*log10(r_c) + 181.6;
	}
}

//printf("r_deltas_p = %f u =%f \n",r_deltas_p,u);
if(r_deltas_p <= 5000)
{
	deltak_1 = alpha*(1.43*log10(r_deltas_p) - 5.29);
}else
{
	deltak_1 = 0.0;
}

// gamma and beta based on Mach number

gamm = 27.094*mach + 3.31;
gamm_0 = 23.43*mach + 4.651;
beta = 72.65*mach + 10.74;
beta_0 = -34.19*mach - 13.82;

// determine k_2 based on alpha

if (alpha < gamm_0-gamm)
{
	k_2 = k_1-1000.0;
} else
{
	if(alpha > gamm_0+gamm)
	{
		k_2 = k_1-12.0;
	} else
	{
		k_2 = k_1+sqrt(pow(beta,2) - pow(beta/gamm,2)*pow(alpha-gamm_0,2)) + beta_0;
	}
}
//if(itime==270&&isurf==6) printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",isect,rpnorm,theta,phi,c,lsect,alpha,u,del_s,del_p,dis_s,dis_p,r_c,m_c,d_h,d_l);


//printf("spl_p and spl_s:\n");
for (int i=0; i<nf2; i++)
{	// pressure side
	spl_p[i] = 10.0*log10((dis_p*pow(mach,5)*lsect*d_h)/(pow(rpnorm,2.0))) + aap[i] + (k_1-3.0) + deltak_1;
	// suction side
	spl_s[i] = 10.0*log10((dis_s*pow(mach,5)*lsect*d_h)/(pow(rpnorm,2.0))) + aas[i] + (k_1-3.0);
//if(isurf==6&&itime==269&&i==nf2-1) printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",isect,spl_p[i],spl_s[i],dis_p,dis_s,mach,lsect,d_h,rpnorm,aap[i],aas[i],k_1,deltak_1,r_c);

}


if(alpha >= gamm_0 || alpha > 12.5)
{
	r_c_new = 3.0*r_c;
	if(r_c_new < 9.52*pow(10,4))
	{
		a_0[0] = 0.57;
	} else
	{
		if(r_c_new > 8.57*pow(10,5)) a_0[0] = 1.13;
		else a_0[0] = -9.57*pow(10.0,-13.0)*pow(r_c_new-(8.57*pow(10,5)),2) + 1.13;
	}

	if(r_c_new < 2.47*pow(10,5))
	{
		k_1 = -4.31*log10(r_c_new) + 156.3;
	} else
	{
		if(r_c_new > 8*pow(10,5)) k_1 = 128.5;
		else k_1 = -9.0*log10(r_c_new) + 181.6;
	}

	if(alpha < gamm_0-gamm)
	{
		k_2 = k_1-1000.0;
	} else
	{
		if(alpha > gamm_0+gamm) k_2 = k_1-12.0;
		else k_2 = k_1 + sqrt(pow(beta,2) - (pow(beta/gamm,2))*(pow(alpha-gamm_0,2))) + beta_0;
	}

	a_min(a_0,min_0,1);
	a_max(a_0,max_0,1);
	
	a_r = (-20.0-min_0[0])/(max_0[0]-min_0[0]);

	a_min(anew,aminb,nf2);
	a_max(anew,amaxb,nf2);

	for(int i=0; i<nf2; i++)
	{
		anew2[i] = aminb[i] + a_r*(amaxb[i]-aminb[i]);
		spl_a[i] = 10.0*log10((dis_s*pow(mach,5)*lsect*d_l)/(pow(rpnorm,2))) + anew2[i] +k_2;
	}


} else
{
	for(int i=0; i<nf2; i++) spl_a[i] = 10.0*log10((dis_s*pow(mach,5)*lsect*d_h)/(pow(rpnorm,2))) + bb[i] +k_2;

//printf("spl_a noise-2\n");
//for(int i=0; i<nf2; i++) printf("%f %f %f %f %f %f %f %f\n",spl_a[i],dis_s,mach,lsect,d_h,rpnorm,bb[i],k_2);	
}

//if(isurf==6&&itime==269) printf("%d %f %f %f %f %f %f %f %f %f %f %f %f\n",isect,spl_a[nf2-1],dis_s,mach,lsect,d_h,rpnorm,bb[nf2-1],k_2,anew2[nf2-1],alpha,gamm_0,10.0*log10((dis_s*pow(mach,5)*lsect*d_h)/(pow(rpnorm,2))));


// clipping as the same as original model
// BPM model in NASA report

for (int i=0; i<nf2; i++)
{
	if(spl_a[i] < -100) spl_a[i] = -100.0;
	if(spl_s[i] < -100) spl_s[i] = -100.0;
	if(spl_p[i] < -100) spl_p[i] = -100.0;
}




for (int i=0; i<nf2; i++)
{
	output_tmp[optr + i*7 + 1] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+1]/10.0) + pow(10.0,spl_a[i]/10.0));
	output_tmp[optr + i*7 + 2] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+2]/10.0) + pow(10.0,spl_s[i]/10.0));
	output_tmp[optr + i*7 + 3] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+3]/10.0) + pow(10.0,spl_p[i]/10.0));

	output_tmp[optr + i*7] = 10.0*log10(pow(10.0,output_tmp[optr+i*7]/10.0) + pow(10.0,spl_a[i]/10.0) + pow(10.0,spl_p[i]/10.0) + pow(10.0,spl_s[i]/10.0));

//printf("%f %f %f",spl_a[i],spl_s[i],spl_p[i]);
}



// LBL-VS noise

if(trip==0)
{
	if(r_c < 1.3*pow(10,5))	stlprim = 0.18;
	else if(r_c>1.3*pow(10,5) && r_c <= 4*pow(10,5)) stlprim = 0.001756*pow(r_c,0.3931);
	else stlprim = 0.28;

	stpkprm = pow(10,(-0.04*alpha))*stlprim;

	if(alpha<=3.0) r_c0 = pow(10,(0.215*alpha+4.978));
	else r_c0 = pow(10,(0.120*alpha+5.263));

	// compute peak scaled spectrum level

	d = r_c/r_c0;

	if(d<=0.3237) g2 = 77.852*log10(d)+15.328;
	else if(d<=0.5689) g2 = 65.188*log10(d)+9.125;
	else if(d<=1.7579) g2 = -114.052*pow(log10(d),2);
	else if(d<=3.0889) g2 = -65.188*log10(d)+9.125;
	else g2 = -77.852*log10(d)+15.328;

	g3 = 171.04 - 3.03*alpha;
	scal = 10.0*log10(del_p*pow(mach,5)*d_h*lsect/(pow(rpnorm,2)));

	for (int i=0; i<nf2; i++)
	{
		st_p[i] = ff[i]*del_p/u;
		e[i] = st_p[i]/stpkprm;
	}

	for (int i=0; i<nf2; i++)
	{
		if(e[i]<=0.5974) g1[i] = 39.8*log10(e[i])-11.12;
		else if(e[i]<=0.8545) g1[i] = 98.409*log10(e[i])+2.0;
		else if(e[i]<=1.17) g1[i] = -5.076+sqrt(2.484-506.25*(pow(log10(e[i]),2)));
		else if(e[i]<=1.674) g1[i] = -98.409*log10(e[i])+2.0;
		else g1[i] = -39.80*log10(e[i])-11.12;

		spl_lam = g1[i] +g2+g3+scal;
	
		output_tmp[optr + i*7 + 4] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+4]/10.0) + pow(10.0,spl_lam/10.0));
		output_tmp[optr + i*7] = 10.0*log10(pow(10.0,output_tmp[optr+i*7]/10.0) + pow(10.0,spl_lam/10.0));

//if(isect==4&&isurf==0&itime==300) printf("%d %d %f %f %f %f %f %f %f\n",isect,itime,ff[i],spl_lam,g1[i],e[i],st_p[i],u,del_p);
//if(isurf==0&&itime==0&&i==34) printf("%d %f %f %f %f %f %f %f\n",isect,spl_lam,r_c,alpha,del_p,theta,phi,mach);




	}	
}

// TEB-VS

float hdstar,hdstarl,hdstarp,aterm,dstarh,dstravg,g4,g5,g50,g514,f4temp,stpeak,stppp[nf2],eta[nf2];
float h = 0.0005; //[m]
float psi = 14.0;

dstravg = (dis_s+dis_p)/2.0; // average displacement thickness

hdstar = h/dstravg;
dstarh = 1.0/hdstar;

// compute peak strouhal number

aterm = 0.212 - 0.0045*psi;
if(hdstar >= 0.2) stpeak = aterm/(1.0+0.235*dstarh-0.0132*pow(dstarh,2));
else stpeak = 0.1*hdstar + 0.095 - 0.00243*psi;

// compute scaled spectrum level

if(hdstar <= 5.0) g4 = 17.5*log10(hdstar)+157.5-1.114*psi;
else g4 = 169.7 - 1.114*psi;

// strouhal number

for (int i=0; i<nf2; i++) 
{
	stppp[i] = ff[i]*h/u;
	eta[i] = log10(stppp[i]/stpeak);
}

for (int i=0; i<nf2; i++) 
{
	hdstarl = hdstar;
	g514 = g5comp(hdstarl,eta[i]);

	hdstarp = 6.724*(pow(hdstar,2)) - 4.019*hdstar + 1.107;

	g50 = g5comp(hdstarp,eta[i]);

	g5 = g50 + 0.0714*psi*(g514-g50);

	if(g5>0) g5=0.0;

	f4temp = g5comp(0.25,eta[i]);

	if(g5>f4temp) g5=f4temp;

	scal = 10.0*log10(pow(mach,5.5)*h*d_h*lsect/(pow(rpnorm,2)));

	spl_blunt = g4+g5+scal;
	output_tmp[optr + i*7 + 5] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+5]/10.0) + pow(10.0,spl_blunt/10.0));
	output_tmp[optr + i*7] = 10.0*log10(pow(10.0,output_tmp[optr+i*7]/10.0) + pow(10.0,spl_blunt/10.0));


}



// TVF (tip noise)
float mm,um,term,tip_l,stpp[nf2],AR,alprat;

AR = R/c;
float AR0[6] = {2.0,2.67,4.0,6.0,12.0,24.0};
float alprat0[6] = {0.54,0.62,0.71,0.79,0.89,0.95};



if(AR < AR0[0]){
	alprat = alprat0[0];
}else if(AR > AR0[5]){
	alprat = alprat0[5];
}else{
	int kk = 0;
	for(int i=0; i<6; i++){
		if(AR < AR0[i]){
			break;
		}
		kk += 1;
	}
	alprat = alprat0[kk-1]*(AR0[kk]-AR)/(AR0[kk]-AR0[kk-1]) + alprat0[kk]*(AR-AR0[kk-1])/(AR0[kk]-AR0[kk-1]);
}


		
		
if(isect == nsect-1)
{
	int iround = 0;
	

	float alptipp = alpha*alprat; // deg

	if(iround==1)
	{
		tip_l = 0.008*alptipp*c;
	}else
	{
		if(abs(alptipp)<=2.0) tip_l=(0.023+0.0169*alptipp)*c;
		else tip_l=(0.0378+0.0095*alptipp)*c;
	}

	mm = (1.0+0.036*alptipp)*mach;
	um = mm*a0;

	term = mach*mach*pow(mm,3)*pow(tip_l,2)*d_h/(pow(rpnorm,2));

	term = MAX(term,2.2204460492503131e-16);

	scal = 10.0*log10(term);

//printf("tip noise\n");
//printf("%f %f %f %f %f\n\n",scal,term,um,alpha,mach);
	for (int i=0; i<nf2; i++)
	{
		stpp[i] = ff[i]*tip_l/um;
		spl_tip = 126.0-30.5*pow((log10(stpp[i])+0.3),2) + scal;
		output_tmp[optr + i*7 + 6] = 10.0*log10(pow(10.0,output_tmp[optr+i*7+6]/10.0) + pow(10.0,spl_tip/10.0));
		output_tmp[optr + i*7] = 10.0*log10(pow(10.0,output_tmp[optr+i*7]/10.0) + pow(10.0,spl_tip/10.0));


	}

}

//free(pmin);free(pmax);free(smin);free(smax);free(aas);free(aap);free(b);
//free(bmin);free(bmax);free(bb);free(spl_a);free(spl_p);free(spl_s);
}









