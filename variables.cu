#include "variables.cuh"

#ifndef VARS_H
#define VARS_H 

// size parameters
int nXtot,nNtot,nVtot,nUtot,nPRtot,tX,tN,tV,tU,tPR,nEtot;
int* qPtr,* XPtr,* NPtr,* VPtr,* UPtr,* PRPtr,* E;

int nX = 3;
int nXham = 4;
int nNham = 6;
int nVham = 6;
int nUham = 6;
int nPRham;

int iXham      = 0;
int iYham      = 1;
int iZham      = 2;
int iDSham     = 3;
int iNXham     = 0;
int iNYham     = 1;
int iNZham     = 2;
int iNdXham    = 3;
int iNdYham    = 4;
int iNdZham    = 5;
int iVXham     = 0;
int iVYham     = 1;
int iVZham     = 2;
int iVdXham    = 3;
int iVdYham    = 4;
int iVdZham    = 5;
int iUXham     = 0;
int iUYham     = 1;
int iUZham     = 2;
int iUdXham    = 3;
int iUdYham    = 4;
int iUdZham    = 5;
int iPRham     = 0;
int iPRdham    = 1;
int iRHOham    = 2;
int iRHOdham   = 3;

// input parameters
float a0,rhoRef,pRef,oM,dPsi,Minf,Mref,xy_angle,xz_angle,CFDscale,OBSscale;
int nSurf,nSurf_dum,periodic,lowpass,impermeable;

// code variables

float vInf[3];
float a0I,oMI,dTau;
int nObs,nTime;
float* rObs,* tObs;
float* Xham,* Nham,* Vham,* Uham,* PRham;

// Broadband noise
//int* xsPtr,* dsPtr*, aPtr,* mPtr,* rPtr,* sPtr,* dPtr,* ccw,* trip,* nsect,* nsect_t,* nmach_t;
int* nre_t,* naoa_t,* nsect,* nsect_t,* nmach_t,* xsPtr,* dsPtr,* aPtr,* mPtr,* rPtr,* sPtr,* dPtr;
int* ccw,* trip;
float* xyzcR,* bbdata,* span,* omega,* psi_offset,* rotate,* translate,* sect_t,* aoa_t,* mach_t;
float* re_t,* bbdata_t,* ff;
int BB_noise,nSurf_BB,nf2,t_del,totsect0,totmach,totsect,totaoa,totre,totdata;

#endif
