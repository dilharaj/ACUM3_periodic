#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
using std::ifstream;
#include <string>
using std::string;
#include "variables.cuh"
#include "mathfunc.cuh"

void read_inputs(const char* inp_file){

     string dat;

     ifstream inpfile(inp_file);

     if (!inpfile){
        printf("No input file found\n");
        exit(1);
     }
 
     printf("*** Reading input file\n\n");
     
     // in case, input file does not have the variable
     
     // Default values
     a0 = 340.0;
     rhoRef = 1.2;
     pRef = 0.0; // Compute from other variables
     nTime = 360;
     dPsi = 1.0;
     xy_angle = 0.0;
     xz_angle = 0.0;
     periodic = 1;
     lowpass = 1;
     BB_noise = 0;
     Mref = 0.08; 
     int n = 0;
     while (!inpfile.eof()){
        n++;
        inpfile >> dat;
	
       	if (dat=="a0") inpfile >> a0;
	else if (dat=="impermeable") inpfile >> impermeable;
        else if (dat=="rhoRef") inpfile >> rhoRef;
        else if (dat=="pRef") inpfile >> pRef;
        else if (dat=="oM") inpfile >> oM;
        else if (dat=="nSurf") inpfile >> nSurf; // Number of surfaces
	else if (dat=="nTime") inpfile >> nTime; // Number of Time steps
	else if (dat=="dPsi") inpfile >> dPsi;
	else if (dat=="Minf") inpfile >> Minf; // Free stream Mach number
	else if (dat=="Mref") inpfile >> Mref; // Free stream Mach number
	else if (dat=="xy_angle") inpfile >> xy_angle;//Minf with x-axis in xy plane
	else if (dat=="xz_angle") inpfile >> xz_angle;//Project Minf on xy-plane
	else if (dat=="CFDscale") inpfile >> CFDscale;//length dimensionalization ratio
	else if (dat=="OBSscale") inpfile >> OBSscale; // Reference length to dimentionalize observer position vectors
	else if (dat=="periodic") inpfile >> periodic;
	else if (dat=="lowpass") inpfile >> lowpass;
	else if (dat=="BB_noise") inpfile >> BB_noise;
	else if (dat=="") continue;
	else if (dat=="ACUM-3_INPUTS") continue;
	else if (dat=="***Environmental") continue;
	else if (dat=="***Geometric") continue;
	else if (dat=="***Data") continue;
        else {printf("Input parameter not defined (%d)\n",n); exit(1);}
     }
    

    if(impermeable==0) nPRham = 4;
    else nPRham = 2;  

    oM = abs(oM); 
    dTau = 2*M_PI*dPsi/oM/360;
    dPsi = dPsi*M_PI/180.0;
    a0I = 1.0/a0;
    //vInf[0] = Minf*a0*cos(xz_angle*M_PI/180.0)*cos(xy_angle*M_PI/180.0);
    //vInf[1] = Minf*a0*cos(xz_angle*M_PI/180.0)*sin(xy_angle*M_PI/180.0);
    //vInf[2] = Minf*a0*sin(xz_angle*M_PI/180.0);
   
    vInf[0] = Minf*a0*cos(xy_angle*M_PI/180.0)*cos(xz_angle*M_PI/180.0);
    vInf[1] = Minf*a0*sin(xy_angle*M_PI/180.0);
    vInf[2] = Minf*a0*cos(xy_angle*M_PI/180.0)*sin(xz_angle*M_PI/180.0);
    
 
    printf("	Number of Surfaces = %d\n",nSurf);
    printf("   	Time Step Size = %f\n",dTau);
    printf("  	Number of tipe steps = %d\n",nTime);
    printf("	Freestream Mach Number = %f\n	Angle of Airflow with\n\t\tx-y plane           = %f degrees\n\t\tx-axis in x-y plane = %f degrees\n",Minf,xz_angle,xy_angle);
    printf("        CFD Length Scale = %f\n",CFDscale);
    printf("        Reference Length for Observers = %f\n\n",OBSscale);

    printf(" 	Periodic Inputs  = %d\n",periodic);
    printf("	Low Pass Filter = %d\n",lowpass);


    
     

      printf("\n\tDone reading input file\n\n");
}



// Read file with observers location
void read_observers(const char* filname){

      FILE *fid = fopen(filname,"r");
      if (!fid)
      {
      	printf("Unable to open observer file: %s\n",filname);
      	exit(1);
      }

      fscanf(fid,"%d",&nObs);
      printf("*** Reading observers (%d)\n\n",nObs);
      rObs = (float*)malloc(nObs*nX*sizeof(float));
     
      float rx,ry,rz;
      float cs = 1.0;
      float ss = 0.0;
    
      for (int i=0;i<nObs;i++){
		 fscanf(fid,"%f %f %f",&rx,&ry,&rz);
		 rObs[i*nX+0] = rx*OBSscale;
		 rObs[i*nX+1] = (ry*cs - rz*ss)*OBSscale;
		 rObs[i*nX+2] = (ry*ss + rz*cs)*OBSscale;
		 
      }
	
      fclose(fid);
      
      printf("\tDone reading observers\n\n");

}

void writeTimeHistory(float* pT,float* pL,float* pA)
{
      
     FILE *fid;
     
// thickness
      fid = fopen("pT.out","w");
      fprintf(fid,"# ntime nobs\n");
      fprintf(fid,"%d %d\n",nTime,nObs);
      
      for (int i=0;i<nTime;i++)
      {
         fprintf(fid," %14.8E",i*dPsi);
         for (int j=0;j<nObs;j++)
	 {
            fprintf(fid," %14.8E",pT[i*nObs+j]);
         }
         fprintf(fid,"\n");
      }
 
      fclose(fid);

//Loading
         fid = fopen("pL.out","w");
	 fprintf(fid,"# ntime nobs\n");
	 fprintf(fid,"%d %d\n",nTime,nObs);


         for (int i=0;i<nTime;i++)
	 {
            fprintf(fid," %14.8E",i*dPsi);
            for (int j=0;j<nObs;j++)
	    {
               fprintf(fid," %14.8E",pL[i*nObs+j]);
            }
            fprintf(fid,"\n");
         }
 
         fclose(fid);


//Total

         fid = fopen("pA.out","w");
	 fprintf(fid,"# ntime nobs\n");
	 fprintf(fid,"%d %d\n",nTime,nObs);


         for (int i=0;i<nTime;i++)
	 {
            fprintf(fid," %14.8E",i*dPsi);
            for (int j=0;j<nObs;j++){
               fprintf(fid," %14.8E",pA[i*nObs+j]);
            }
            fprintf(fid,"\n");
         }
 
         fclose(fid);


//	fid = fopen("test.dat","w");
//	for(int i=0;i<nTime;i++)
//	{
//		fprintf(fid,"%d 1.0 %f %f 1.0\n",i+1,pL[i*nObs+45],pL[i*nObs+90]);
//	}

}


void read_surf(int nSurf)
{
	#define filenametemplate "surfaces/surface_%d.dat"
	XPtr = (int*)malloc(nSurf*sizeof(int));
	XPtr[0] = 0;
	NPtr = (int*)malloc(nSurf*sizeof(int));
	NPtr[0] = 0;
	VPtr = (int*)malloc(nSurf*sizeof(int));
	VPtr[0] = 0;
	UPtr = (int*)malloc(nSurf*sizeof(int));
	UPtr[0] = 0;
	PRPtr = (int*)malloc(nSurf*sizeof(int));
	PRPtr[0] = 0;


	E = (int*)malloc(nSurf*sizeof(int));
	int etot,iXh,iNh,iUh,iPRh,nt;
	//int iVh;
	float* Area = (float*)malloc(nSurf*sizeof(float));
	etot = 0;
	//static const char* const filename[] = {"surfaces/surface_1.dat","surfaces/surface_2.dat","surfaces/surface_3.dat","surfaces/surface_4.dat","surfaces/surface_5.dat","surfaces/surface_6.dat","surfaces/surface_7.dat","surfaces/surface_8.dat","surfaces/surface_9.dat","surfaces/surface_10.dat"};
	
	for (int i=0; i<nSurf; i++)
	{
		char fname[100];
		sprintf(fname,filenametemplate,i+1);

		printf("%s\n",fname);
		FILE* fid = fopen(fname,"r");
		fscanf(fid, "%d %d\n", &E[i],&nt);
		fclose(fid);
		if (nt != nTime)
		{
			printf("	Inconsistent number of steps in %s\n",fname);
		}
		
		etot += E[i];
		if(i<nSurf-1){
			XPtr[i+1] = XPtr[i] + nXham*E[i]*nTime;
			NPtr[i+1] = NPtr[i] + nNham*E[i]*nTime;
			VPtr[i+1] = VPtr[i] + nVham*E[i]*nTime;
			UPtr[i+1] = UPtr[i] + nUham*E[i]*nTime;
			PRPtr[i+1] = PRPtr[i] + nPRham*E[i]*nTime;
		}
	}

 	nXtot = etot*nXham*nTime;
	Xham = (float*)malloc(nXtot*sizeof(float));
	nNtot = etot*nNham*nTime;
	Nham = (float*)malloc(nNtot*sizeof(float));
	nVtot = etot*nVham*nTime;
	Vham = (float*)malloc(nVtot*sizeof(float));
	nUtot = etot*nUham*nTime;
	Uham = (float*)malloc(nUtot*sizeof(float));
	nPRtot = etot*nPRham*nTime;
	PRham = (float*)malloc(nPRtot*sizeof(float));

	if(impermeable==0){	
	  
	  for (int i=0; i<nSurf;i++)
	  {
		
       		char fname[100];
		sprintf(fname,filenametemplate,i+1);


		FILE* fid = fopen(fname,"r");

	    fscanf(fid, "%d %d\n", &E[i],&nt);

	    for (int k=0; k<nTime; k++)
	    {  
		for (int j=0; j<E[i]; j++)
		{
		  iXh = XPtr[i] + j*nTime*nXham + k*nXham;
		  iNh = NPtr[i] + j*nTime*nNham + k*nNham;
  		  //iVh = VPtr[i] + j*nTime*nVham + k*nVham;
		  iUh = UPtr[i] + j*nTime*nUham + k*nUham;
		  iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;

     	          fscanf(fid,"%f %f %f %f %f %f %f %f %f %f %f %f\n",&Xham[iXh+iXham],&Xham[iXh+iYham],&Xham[iXh+iZham],&Uham[iUh+iUXham],&Uham[iUh+iUYham],&Uham[iUh+iUZham],&PRham[iPRh+iPRham],&PRham[iPRh+iRHOham],&Xham[iXh+iDSham],&Nham[iNh+iNXham],&Nham[iNh+iNYham],&Nham[iNh+iNZham]);

		}
	    }
	    fclose(fid);
	    printf("	Done reading %s with %d elements\n",fname,E[i]);
	
          }

	}else if(impermeable==1){

	  for (int i=0; i<nSurf;i++)
	  {
		
            	char fname[100];
		sprintf(fname,filenametemplate,i+1);


		FILE* fid = fopen(fname,"r");

	    fscanf(fid, "%d %d\n", &E[i],&nt);

	    for (int k=0; k<nTime; k++)
	    {  
		for (int j=0; j<E[i]; j++)
		{
		  iXh = XPtr[i] + j*nTime*nXham + k*nXham;
		  iNh = NPtr[i] + j*nTime*nNham + k*nNham;
  		  //iVh = VPtr[i] + j*nTime*nVham + k*nVham;
	//	  iUh = UPtr[i] + j*nTime*nUham + k*nUham;
		  iPRh = PRPtr[i] + j*nTime*nPRham + k*nPRham;

     	          fscanf(fid,"%f %f %f %f %f %f %f %f\n",&Xham[iXh+iXham],&Xham[iXh+iYham],&Xham[iXh+iZham],&PRham[iPRh+iPRham],&Xham[iXh+iDSham],&Nham[iNh+iNXham],&Nham[iNh+iNYham],&Nham[iNh+iNZham]);

		}


	    }
	    fclose(fid);
	    printf("	Done reading %s with %d elements\n",fname,E[i]);

	
	
	  }
	}
	free(Area);	
	printf("\n");
}



void read_BB_inputs(const char* inp_file)
{
	#define fnametemplate "BB_inputs/BB_surface_%d.dat"
	#define surf_inp "***Surface-%d"
	#define table_inp "BB_inputs/BB_table_%d.dat"

	//// READ INPUT FILE ////

	
	string dat;
	
	ifstream inpfile(inp_file);

     	if (!inpfile){
        	printf("No broadband input file found\n");
        	exit(1);
    	}
 
     	printf("*** Reading broadband input file\n\n");

	nSurf_BB = nSurf; // initialize nSurf_BB to be equal to nSurf
     	t_del    = 0; // 0 to read BL thicknesses from surface file; 1 to read from table	
	// Allocate memory
	span = (float*)malloc(nSurf_BB*sizeof(float));
	omega = (float*)malloc(nSurf_BB*sizeof(float));
	ccw = (int*)malloc(nSurf_BB*sizeof(int));
	psi_offset = (float*)malloc(nSurf_BB*sizeof(float));
	trip = (int*)malloc(nSurf_BB*sizeof(int));
	rotate = (float*)malloc(3*nSurf_BB*sizeof(float));
	translate= (float*)malloc(3*nSurf_BB*sizeof(float));

//	int* check = (int*)malloc(nSurf_BB*sizeof(int));

	// Initialize
	for (int i=0; i<nSurf_BB; i++){
		span[i] = 1.0;
		omega[i] = 1.0;
		ccw[i] = 1;
		psi_offset[i] = 0.0;
		trip[i] = 0;
		for (int j=0; j<3; j++){
			rotate[3*i+j] = 0.0;
			translate[3*i+j] = 0.0;
		}
	}

	// Read data

	int n = 0;
	int is = 0;
	char surf[15];

	while (!inpfile.eof()){
	   n++;
	   inpfile >> dat;
	   sprintf(surf,surf_inp,is+1);
	   if (dat=="ACUM-3_Broadband_Inputs") continue;
	   else if (dat=="") continue;
	   else if (dat=="nSurf_BB"){ inpfile >> nSurf_BB; if (nSurf_BB>nSurf){ printf("Number of surfaces considered for broadband noise cannot be higher than the available surfaces\n"); exit(1);}}
	   else if (dat=="t_del") inpfile >> t_del;
	   else if (dat==surf) is += 1;
	   else if (dat=="span(m)") inpfile >> span[is-1];
	   else if (dat=="omega(rad/s)") inpfile >> omega[is-1];   
	   else if (dat=="ccw") inpfile >> ccw[is-1];
	   else if (dat=="psi_offset") inpfile >> psi_offset[is-1];
	   else if (dat=="trip") inpfile >> trip[is-1];
	   else if (dat=="Initial_Position") continue;
	   else if (dat=="rotate"){
	   	for (int i=0; i<3; i++){
			inpfile >> rotate[3*(is-1)+i];
		}
	   }
  	   else if (dat=="translate"){
	   	for (int i=0; i<3; i++){
			inpfile >> translate[3*(is-1)+i];
		}
		if (is==nSurf_BB) break;
	   }

  	   else {printf("Input parameter not defined (%d)\n",n); exit(1);}
	}

	
	
	// print out inmput parameters
	printf("Surfaces considered for Broadband noise calculation: %d\n",nSurf_BB);
	if (t_del==1) printf("\nBL thicknesses taken from tables.\n\n");

	for (int i=0; i<nSurf_BB; i++){
		printf("\tSurface-%d:\n",i+1);
		printf("\t\tspan (m)         = %f\n",span[i]);
		printf("\t\tomega (rad/s)    = %f\n",omega[i]);
		printf("\t\tccw              = %d\n",ccw[i]);
		printf("\t\tpsi offset (deg) = %f\n",psi_offset[i]);
		printf("\t\tBL trip          = %d\n",trip[i]);
		printf("\t\tInitial Position:\n");
		printf("\t\t\tRotate    = [%f %f %f]\n",rotate[3*i],rotate[3*i+1],rotate[3*i+2]);
		printf("\t\t\tTranslate = [%f %f %f]\n",translate[3*i],translate[3*i+1],translate[3*i+2]);	


		if(omega[i]!=0.0) printf("cn and ct normalization will be adjusted based on span\n");	
	
	}

		
	
	// Read BB_surface data

	xsPtr = (int*)malloc(nSurf_BB*sizeof(int));
	xsPtr[0] = 0;
	dsPtr = (int*)malloc(nSurf_BB*sizeof(int));
	dsPtr[0] = 0;

	nsect = (int*)malloc(nSurf_BB*sizeof(int));
	int nt,lptr;

	totsect0 = 0;
	for (int i=0; i<nSurf_BB; i++)
	{
		char fname[100];
		sprintf(fname,fnametemplate,i+1);
		printf("%s\n",fname);
		FILE* fid = fopen(fname,"r");
		fscanf(fid, "%d %d\n", &nt,&nsect[i]);
		fclose(fid);
		totsect0 = totsect0 + nsect[i];

		if(i<nSurf_BB-1){
			xsPtr[i+1] = xsPtr[i] + nsect[i]*4;
			dsPtr[i+1] = dsPtr[i] + nsect[i]*7*nTime;
		}

		if (nt != nTime)
		{
			printf("	Inconsistent number of steps in %s\n",fname);
		}
	}
				
	xyzcR = (float*)malloc(totsect0*4*sizeof(float)); // y and z locs are also included
	bbdata = (float*)malloc(totsect0*nTime*7*sizeof(float));
        
	float* cn = (float*)malloc(20*sizeof(float));
	float* cc = (float*)malloc(20*sizeof(float));

	for(int ij=0; ij<20; ij++){
		cn[ij] = 0.0;
		cc[ij] = 0.0;
	}


	 
	for (int i=0; i<nSurf_BB; i++)
	{
		char fname[100];
		sprintf(fname,fnametemplate,i+1);
		FILE* fid = fopen(fname,"r");
		fscanf(fid, "%d %d\n", &nt,&nsect[i]);
		
		for (int j=0; j<nsect[i]; j++) fscanf(fid, "%e %e %e %e\n",&xyzcR[xsPtr[i]+4*j],&xyzcR[xsPtr[i]+4*j+1],&xyzcR[xsPtr[i]+4*j+2],&xyzcR[xsPtr[i]+4*j+3]);
		
		for (int j=0; j<nTime; j++)
		{
		    for (int k=0; k<nsect[i]; k++)
		    {
			lptr = dsPtr[i] + j*nsect[i]*7 + k*7;
			fscanf(fid,"%e %e %e %e %e %e %e\n",&bbdata[lptr],&bbdata[lptr+1],&bbdata[lptr+2],&bbdata[lptr+3],&bbdata[lptr+4],&bbdata[lptr+5],&bbdata[lptr+6]);
		        if(i==0){
			//	printf("%d %f %f\n",j,bbdata[lptr+4],bbdata[lptr+5]);
			cn[k] += bbdata[lptr+4]/nTime;
			cc[k] += bbdata[lptr+5]/nTime;
			}
	// multipy by r/R**2 to adjust for Vref change assuming CFD CP uses reference tip Mach number
			if(omega[i]!=0.0){	
	   		  bbdata[lptr+4] = bbdata[lptr+4]/ pow(xyzcR[xsPtr[i] + 4*k],2);
	    		  bbdata[lptr+5] = bbdata[lptr+5]/ pow(xyzcR[xsPtr[i] + 4*k],2);
			}
		    }
		}
		fclose(fid);
	//	for(int ij=0; ij<20; i++){
	//		printf("%d %f %f\n",ij,cn[ij],cc[ij]);
	//	}
	

	}

	// Read aerodynamic tables
	nsect_t = (int*)malloc(nSurf_BB*sizeof(int));
	nmach_t = (int*)malloc(nSurf_BB*sizeof(int));
	naoa_t = (int*)malloc(nSurf_BB*sizeof(int));
	nre_t = (int*)malloc(nSurf_BB*sizeof(int));

	aPtr = (int*)malloc(nSurf_BB*sizeof(int));
	aPtr[0] = 0;
	mPtr = (int*)malloc(nSurf_BB*sizeof(int));
	mPtr[0] = 0;
	rPtr = (int*)malloc(nSurf_BB*sizeof(int));
	rPtr[0] = 0;
	sPtr = (int*)malloc(nSurf_BB*sizeof(int));
	sPtr[0] = 0;
	dPtr = (int*)malloc(nSurf_BB*sizeof(int));
	dPtr[0] = 0;


	totmach = 0;
	totaoa = 0;
	totre = 0;
	totsect = 0;
	totdata = 0;

	
	for (int i=0; i<nSurf_BB; i++)
	{
		char fname[100];
		sprintf(fname,table_inp,i+1);
		printf("%s\n",fname);

		FILE* fid = fopen(fname,"r");


		fscanf(fid, "%*[^\n]\n");
		fscanf(fid, "%d %d %d %d\n",&nsect_t[i],&naoa_t[i],&nmach_t[i],&nre_t[i]);   // correct
		//fscanf(fid, "%d %d %d %d\n",&nsect_t[i],&naoa_t[i],&nre_t[i],&nmach_t[i]);   // wrong

		fclose(fid);	

		totmach += nmach_t[i];
		totaoa += naoa_t[i];
		totre += nre_t[i];
		totsect += nsect_t[i];
		totdata += nsect_t[i]*naoa_t[i]*nmach_t[i]*nre_t[i]*7;

		if(i<nSurf_BB-1){
			aPtr[i+1] = aPtr[i] + naoa_t[i];
			sPtr[i+1] = sPtr[i] + nsect_t[i];
			mPtr[i+1] = mPtr[i] + nmach_t[i];
			rPtr[i+1] = rPtr[i] + nre_t[i];
			dPtr[i+1] = dPtr[i] + nsect_t[i]*naoa_t[i]*nmach_t[i]*nre_t[i]*7;
		} 

	}

	printf("\tReading Broadband Table Inputs\n\n");
	for (int i=0; i<nSurf_BB; i++){
		printf("\tSurface-%d:\n",i+1);
		printf("\t\t# of spanwise sections = %d\n",nsect_t[i]);
		printf("\t\t# of Mach numbers      = %d\n",nmach_t[i]);
		printf("\t\t# of Re numbers        = %d\n",nre_t[i]);
		printf("\t\t# of angles of attack  = %d\n",naoa_t[i]);
	}
		
	mach_t = (float*)malloc(totmach*sizeof(float));
	sect_t = (float*)malloc(totsect*sizeof(float));
	re_t = (float*)malloc(totre*sizeof(float));
	aoa_t = (float*)malloc(totaoa*sizeof(float));
	bbdata_t = (float*)malloc(totdata*sizeof(float));

        int lPtr;
	for (int i=0; i<nSurf_BB; i++)
	{
		char fname[100];
		sprintf(fname,table_inp,i+1);
		FILE* fid = fopen(fname,"r");
		fscanf(fid, "%*[^\n]\n");
		fscanf(fid, "%d %d %d %d\n",&nsect_t[i],&naoa_t[i],&nmach_t[i],&nre_t[i]); // correct
		//fscanf(fid, "%d %d %d %d\n",&nsect_t[i],&naoa_t[i],&nre_t[i],&nmach_t[i]);  // wrong

		for (int j=0; j<nsect_t[i]; j++)
		{
		    fscanf(fid, "%f\n",&sect_t[sPtr[i]+j]);
		    
		    for (int k=0; k<nmach_t[i]; k++)
		    {
		        for (int l=0; l<nre_t[i]; l++)
		        {
		            for (int m=0; m<naoa_t[i]; m++)
		            {
		    		lPtr = dPtr[i] + j*nmach_t[i]*nre_t[i]*naoa_t[i]*7 + k*nre_t[i]*naoa_t[i]*7 + l*naoa_t[i]*7 + m*7;
				
				fscanf(fid, "%e %e %e %e %e %e %e %e %e %e\n",&aoa_t[aPtr[i] + m],&mach_t[mPtr[i] + k],&bbdata_t[lPtr+0],&bbdata_t[lPtr+1],&bbdata_t[lPtr+2],&bbdata_t[lPtr+3],&bbdata_t[lPtr+4],&bbdata_t[lPtr+5],&bbdata_t[lPtr+6],&re_t[rPtr[i] + l]);

//if(i==0) printf("%d %f\n",l,re_t[rPtr[i]+l]);
//printf("%f %f %f %f %f %f %f %f %f %f\n",aoa_t[aPtr[i] + m],mach_t[mPtr[i] + k],bbdata_t[lPtr+0],bbdata_t[lPtr+1],bbdata_t[lPtr+2],bbdata_t[lPtr+3],bbdata_t[lPtr+4],bbdata_t[lPtr+5],bbdata_t[lPtr+6],re_t[rPtr[i] + l]);

//if(lPtr==1428||lPtr==2499||lPtr==1071||lPtr==2142) printf("Hi %d\n",lPtr);	
			    }
			}
		    }
		}		
		fclose(fid);	

		// span section check
		for(int j=0; j<nsect_t[i]; j++){
                  if(sect_t[sPtr[i]+j]<0.0||sect_t[sPtr[i]+j]>1.0){
			printf("Check span section %d of surface %d :%f\n",j,i,sect_t[sPtr[i]+j]);
			exit(1);
		  }

		  for(int k=0; k<nsect_t[i]; k++){
		    if(j==k) continue;
		    if(sect_t[sPtr[i]+j]==sect_t[sPtr[i]+k]){
			printf("Check span sections %d and %d of surface %d\n",j,k,i);
			exit(1);
		    }
		  }	
 		}

		// Mach number check
		for(int j=0; j<nmach_t[i]; j++){
                  if(mach_t[mPtr[i]+j]<0.0){
			printf("Check Mach number %d of surface %d\n",j,i);
			exit(1);
		  }

		  for(int k=0; k<nmach_t[i]; k++){
		    if(j==k) continue;
		    if(mach_t[mPtr[i]+j]==mach_t[mPtr[i]+k]){
			printf("Check Mach numbers %d and %d of surface %d\n",j,k,i);
			exit(1);
		    }
		  }	
 		}
	
		// Re check
		for(int j=0; j<nre_t[i]; j++){
                  if(re_t[rPtr[i]+j]<0.0){
			printf("Check Mach number %d of surface %d\n",j,i);
			exit(1);
		  }

		  for(int k=0; k<nre_t[i]; k++){
		    if(j==k) continue;
		    if(re_t[rPtr[i]+j]==re_t[rPtr[i]+k]){
			printf("Check Re numbers %d and %d of surface %d\n",j,k,i);
			exit(1);
		    }
		  }	
 		}
		
		// AoA check
		for(int j=0; j<naoa_t[i]; j++){
               
		  for(int k=0; k<naoa_t[i]; k++){
		    if(j==k) continue;
		    if(aoa_t[aPtr[i]+j]==aoa_t[aPtr[i]+k]){
			printf("Check AoAs %d and %d of surface %d\n",j,k,i);
			exit(1);
		    }
		  }	
 		}
	}

//for(int i=0; i<naoa_t[0]; i++) printf("%f ",aoa_t[aPtr[0]+i]);	
//printf("%f %f %f %f\n",bbdata_t[1428],bbdata_t[2499],bbdata_t[1071],bbdata_t[2142]);	
	printf("Broadband Noise inputs reading completed.\n");

	nf2 = 42;
	ff = (float*)malloc(nf2*sizeof(float));
	float temp[42] = {1.0,1.25,1.6,2.0,2.5,3.15,4.0,5.0,6.3,8.0,10.0,12.5,16.0,20.0,25.0,31.5,40.0,50.0,63.0,80.0,100.0,125.0,160.0,200.0,250.0,315.0,400.0,500.0,630.0,800.0,1000.0,1250.0,1600.0,2000.0,2500.0,3150.0,4000.0,5000.0,6300.0,8000.0,10000.0,20000.0};
	for (int i=0; i<nf2; i++)
	{
		ff[i] = temp[i];
	}
}

void writeBBoutput(float* BBoaspl,float* BBspl)
{

     FILE *fid;
     
// OASPL

      fid = fopen("bb_oaspl.out","w");
      fprintf(fid,"nobs: %d\n",nObs);
      
      for (int i=0;i<nObs;i++) fprintf(fid," %14.8E\n",BBoaspl[i]);
       
      fclose(fid);

// SPL

      fid = fopen("bb_spl.dat","w");
      //fprintf(fid,"nobs: %d  nf: %d\n",nObs,nf2);
      //fprintf(fid,"Total\n");

//      fprintf(fid,"TITLE    = \"BB_SPL\"\nVARIABLES = \"Frequency (Hz)\"\n\"Total (dB)\"\n\"SPL_A (dB)\"\n\"SPL_S (dB)\"\n\"SPL_P (dB)\"\n\"SPL_LBL-VS (dB)\"\n\"SPL_TEB-VS (dB)\"\n\"SPL_TVF (dB)\"\n");
//      
//      for (int i=0;i<nObs;i++)
//      {
//	fprintf(fid,"ZONE T=\"Observer-%d\"\n STRANDID=0, SOLUTIONTIME=0\n I=%d, J=1, K=1, ZONETYPE=Ordered\n DATAPACKING=POINT\n DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )\n",i+1,nf2);
//	
//	for (int j=0; j<nf2; j++){
//	    fprintf(fid," %14.8E",ff[j]);	
//	    for(int k=0; k<7; k++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+k]);
//            fprintf(fid,"\n");
//        }
//      }


	fprintf(fid," # nfrequency nobs\n    %d      %d\n",nf2,nObs);
	fprintf(fid," ZONE T=\"TOTAL\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7]);
		}
		fprintf(fid,"\n");
	}
	fprintf(fid," ZONE T=\"TBL_TE_P\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+3]);
		}
		fprintf(fid,"\n");
	}
	fprintf(fid," ZONE T=\"TBL_TE_S\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+2]);
		}
		fprintf(fid,"\n");
	}
	fprintf(fid," ZONE T=\"TBL_SS\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+1]);
		}
		fprintf(fid,"\n");
	}
 	fprintf(fid," ZONE T=\"TEB_VS\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+5]);
		}
		fprintf(fid,"\n");
	}
	fprintf(fid," ZONE T=\"LBL-VS\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+4]);
		}
		fprintf(fid,"\n");
	}
	fprintf(fid," ZONE T=\"TVF\"\n");
	for(int j=0; j<nf2; j++){
		fprintf(fid,"   %e",ff[j]);
		for(int i=0; i<nObs; i++){
			fprintf(fid,"  %e",BBspl[i*nf2*7+j*7+6]);
		}
		fprintf(fid,"\n");
	}
        fclose(fid);

	free(BBspl);free(BBoaspl);free(span);free(omega);free(ccw);free(psi_offset);free(trip);free(rotate);free(translate);free(xsPtr);free(dsPtr);free(nsect);free(xyzcR);free(bbdata);free(nsect_t);free(nmach_t);free(naoa_t);free(nre_t);free(aPtr);free(mPtr);free(rPtr);free(sPtr);free(dPtr);free(sect_t);free(mach_t);free(aoa_t);free(re_t);free(bbdata_t);free(ff);


//      fprintf(fid,"SPL_A\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+1]);
//        fprintf(fid,"\n");
//      }
//      
//      fprintf(fid,"SPL_S\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+2]);
//        fprintf(fid,"\n");
//      }
// 
//      fprintf(fid,"SPL_P\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+3]);
//        fprintf(fid,"\n");
//      }
//
//      fprintf(fid,"SPL_LBL-VS\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+4]);
//        fprintf(fid,"\n");
//      }
//
//      fprintf(fid,"SPL_TEB-VS\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E",BBspl[i*nf2*7+j*7+5]);
//        fprintf(fid,"\n");
//      }
//
//      fprintf(fid,"SPL_TVF\n");
//      for (int i=0;i<nObs;i++)
//      {
//	for (int j=0; j<nf2; j++) fprintf(fid," %14.8E\n",BBspl[i*nf2*7+j*7+6]);
//        fprintf(fid,"\n");
//      }


}	
