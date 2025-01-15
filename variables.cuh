// size parameters
extern int nXtot,nNtot,nVtot,nUtot,nPRtot,tX,tN,tV,tU,tPR;
extern int* XPtr,* NPtr,* VPtr,* UPtr,* PRPtr,* E;

extern int nX;
extern int nXham;
extern int nNham;
extern int nVham;
extern int nUham;
extern int nPRham;

extern int iXham;
extern int iYham;     
extern int iZham;     
extern int iDSham;    
extern int iNXham;    
extern int iNYham;
extern int iNZham;
extern int iNdXham;
extern int iNdYham;
extern int iNdZham;
extern int iVXham;    
extern int iVYham;    
extern int iVZham;    
extern int iVdXham;
extern int iVdYham;
extern int iVdZham;
extern int iUXham;    
extern int iUYham;    
extern int iUZham;    
extern int iUdXham;
extern int iUdYham;
extern int iUdZham;
extern int iPRham;    
extern int iPRdham;   
extern int iRHOham;   
extern int iRHOdham;  

// input parameters
extern float a0,rhoRef,pRef,oM,dPsi,Minf,Mref,xy_angle,xz_angle,CFDscale,OBSscale;
extern int nSurf,periodic,lowpass,impermeable;

// code variables
extern float vInf[3];  
extern float a0I,oMI,dTau;                       
extern int nObs,nTime;                            
extern float* rObs,* tObs;
extern float* Xham,* Nham,* Vham,* Uham,* PRham;


// Broadband noise
extern int* xsPtr,* dsPtr,* aPtr,* mPtr,* rPtr,* sPtr,* dPtr,* ccw,* trip,* nsect,* nsect_t,* nmach_t,* naoa_t,* nre_t;
extern float* xyzcR,* bbdata,* span,* omega,* psi_offset,* rotate,* translate,* sect_t,* aoa_t,* mach_t,* re_t,* bbdata_t,* ff;
extern int BB_noise,nSurf_BB,nf2,t_del,totsect0,totmach,totsect,totaoa,totre,totdata;






