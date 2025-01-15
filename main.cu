#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "variables.cuh"
#include "io.cuh"
#include "preprocess.cuh"
#include "core.cuh"
#include "bb_core.cuh"
#include "bb_core_gpu.cuh"
int main(void){

    printf("\n\n\n");
    printf("    |--------------------------------------------------|\n");
    printf("    |----- STARTING GPU ACOUSTIC CODE - VERSION 3 -----|\n");
    printf("    |--------------------------------------------------|\n");

    printf("\n\n\n");
    printf("         @           @@@@   @      @  @@       @@        @@@@@\n");
    printf("        @ @         @    @  @      @  @ @     @ @       @     @\n");
    printf("       @   @       @        @      @  @  @   @  @              @\n");
    printf("      @     @     @         @      @  @   @ @   @             @\n");
    printf("     @       @    @         @      @  @    @    @  @@@      @@\n");
    printf("    @ @ @ @ @ @   @         @      @  @         @             @\n");
    printf("   @           @   @        @      @  @         @              @\n");
    printf("  @             @   @    @   @    @   @         @       @     @\n");
    printf(" @               @   @@@@     @@@@    @         @        @@@@@\n\n\n");





    cudaDeviceProp prop;
    int dev = 0; // Choose device here
    cudaDeviceReset();
    cudaSetDevice(dev);
    if (cudaGetDeviceProperties(&prop,dev) == cudaSuccess)
    {
        printf("\n*** Using device %d: %s\n\n",dev,prop.name);
    }
    // Prepare timing
    cudaEvent_t start,stop, stop1, stop2, stop3, stop4,stop5;
    cudaEventCreate(&start);cudaEventCreate(&stop);cudaEventCreate(&stop1);cudaEventCreate(&stop2);    cudaEventCreate(&stop3);cudaEventCreate(&stop4);cudaEventCreate(&stop5);
    float elapsedTime;

    // Start time
    cudaEventRecord(start, 0);

    const char* inp_file = "input";
    const char* observer_file = "obs.dat";
    const char* inp_file_BB = "input_BB";


   
    // reading 
    read_inputs(inp_file);
    if (BB_noise==1) read_BB_inputs(inp_file_BB);

    cudaEventRecord(stop1, 0); cudaEventSynchronize(stop1);
    cudaEventElapsedTime(&elapsedTime, start, stop1);
    printf("        time: %lf ms\n\n",elapsedTime);

    
    read_observers(observer_file);


  
   
    cudaEventRecord(stop2, 0);
    cudaEventSynchronize(stop2);
    cudaEventElapsedTime(&elapsedTime, stop1, stop2);
    printf("        time: %lf ms\n\n",elapsedTime);


    printf("*** Reading Surface(s)\n\n");

    // Read surface_*.dat files
    read_surf(nSurf);

    cudaEventRecord(stop3, 0);
    cudaEventSynchronize(stop3);
    cudaEventElapsedTime(&elapsedTime, stop2, stop3);
    printf("        time: %lf ms\n\n",elapsedTime);


    printf("*** Preprocessing variables\n\n");
    if(impermeable==0) preprocess_per();
    else if(impermeable==1) preprocess_imper();

    printf("	Done.\n\n");
    
    cudaEventRecord(stop4, 0);
    cudaEventSynchronize(stop4);
    cudaEventElapsedTime(&elapsedTime, stop3, stop4);
    printf("        time: %lf ms\n\n",elapsedTime);


 
    // Variable Declaration    
    size_t size_tobs  = nTime*nObs*sizeof(float);
    size_t size_time  = nTime*sizeof(float);
    size_t size_obs  = nObs*3*sizeof(float);
    size_t size_surf = nSurf*sizeof(int);
    size_t size_bb   = nf2*nObs*7*sizeof(float);
    size_t size_bb0   = nObs*sizeof(float);

    float* pT = (float*)malloc(size_tobs);
    float* pL = (float*)malloc(size_tobs);
    float* pA = (float*)malloc(size_tobs);
    float* tObs = (float*)malloc(size_time);
    float* BBspl = (float*)malloc(size_bb);
    float* BBoaspl = (float*)malloc(size_bb0);

    float* pT_gpu; cudaMalloc(&pT_gpu,size_tobs);
    float* pL_gpu; cudaMalloc(&pL_gpu,size_tobs);
    float* pA_gpu; cudaMalloc(&pA_gpu,size_tobs);
    
    float* rObs_gpu; cudaMalloc(&rObs_gpu,size_obs);
    float* tObs_gpu; cudaMalloc(&tObs_gpu,size_time);
    float* vInf_gpu; cudaMalloc(&vInf_gpu,3*sizeof(float));
 
    // Initialization
    for(int t=0;t<nTime;t++){ 
        tObs[t] = t*dPsi;
        for(int iobs=0;iobs<nObs;iobs++){
            pT[t*nObs+iobs] = 0.0;
            pL[t*nObs+iobs] = 0.0;
            pA[t*nObs+iobs] = 0.0;
        }
    }

    for(int iobs=0; iobs<nObs; iobs++)
    {
	BBoaspl[iobs] = -100.0;
	for (int i=0; i<nf2; i++)
	{
	    for (int j=0; j<7; j++)
	    {
	    	BBspl[iobs*nf2*7+i*7+j] = -100.0;
	    }
	}
    }

	// CPU Broadband noise calculation

    if(BB_noise == 1) 
    {
	printf("\n*** Broadband Noise calculation ...\n");

	findBBnoise(nSurf_BB,nObs,nTime,dTau,span,omega,ccw,psi_offset,trip,rotate,translate,xsPtr,dsPtr,nsect,xyzcR,bbdata,nsect_t,nmach_t,naoa_t,nre_t,aPtr,mPtr,rPtr,sPtr,dPtr,sect_t,mach_t,aoa_t,re_t,bbdata_t,t_del,vInf,rObs,tObs,nf2,ff,a0,BBspl,BBoaspl);
    }



 
    cudaMemcpy(vInf_gpu,vInf,3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(pT_gpu,pT,size_tobs,cudaMemcpyHostToDevice);
    cudaMemcpy(pL_gpu,pL,size_tobs,cudaMemcpyHostToDevice);
    cudaMemcpy(pA_gpu,pA,size_tobs,cudaMemcpyHostToDevice);
    cudaMemcpy(tObs_gpu,tObs,size_time,cudaMemcpyHostToDevice);
    cudaMemcpy(rObs_gpu,rObs,size_obs,cudaMemcpyHostToDevice);
    
    
    int* XPtr_g; cudaMalloc(&XPtr_g,size_surf);
    int* NPtr_g; cudaMalloc(&NPtr_g,size_surf);
    int* VPtr_g; cudaMalloc(&VPtr_g,size_surf);
    int* PRPtr_g; cudaMalloc(&PRPtr_g,size_surf);
    int* E_g; cudaMalloc(&E_g,size_surf);
    
    cudaMemcpy(XPtr_g,XPtr,size_surf,cudaMemcpyHostToDevice);
    cudaMemcpy(NPtr_g,NPtr,size_surf,cudaMemcpyHostToDevice);
    cudaMemcpy(VPtr_g,VPtr,size_surf,cudaMemcpyHostToDevice);
    cudaMemcpy(PRPtr_g,PRPtr,size_surf,cudaMemcpyHostToDevice);
    cudaMemcpy(E_g,E,size_surf,cudaMemcpyHostToDevice);

    float* Xham_gpu; cudaMalloc(&Xham_gpu,nXtot*sizeof(float));
    float* Nham_gpu; cudaMalloc(&Nham_gpu,nNtot*sizeof(float));
    float* Vham_gpu; cudaMalloc(&Vham_gpu,nVtot*sizeof(float));
    float* PRham_gpu; cudaMalloc(&PRham_gpu,nPRtot*sizeof(float));
    
    cudaMemcpy(Xham_gpu,Xham,nXtot*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Nham_gpu,Nham,nNtot*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Vham_gpu,Vham,nVtot*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(PRham_gpu,PRham,nPRtot*sizeof(float),cudaMemcpyHostToDevice);

    if(impermeable==0){

    	int* UPtr_g; cudaMalloc(&UPtr_g,size_surf);
	cudaMemcpy(UPtr_g,UPtr,size_surf,cudaMemcpyHostToDevice);
	float* Uham_gpu; cudaMalloc(&Uham_gpu,nUtot*sizeof(float));
	cudaMemcpy(Uham_gpu,Uham,nUtot*sizeof(float),cudaMemcpyHostToDevice);
    }

    // Parallelization
    int ntpb_x = 8; // time steps
    int ntpb_y = 4; // observers
    int ntpb_z = 8;  // elements
    int nbx_max = 512;
    int nby_max = 512;
    int nbz_max = 512;


    // Time steps
    int nb_t = round(float(nTime/ntpb_x) + 0.5f);
    if (nb_t>nbx_max) nb_t = nbx_max;
    int ntpb_t = ntpb_x;
    if (nTime<ntpb_t){
        ntpb_t = nTime;
        nb_t = 1;
    }
    // Observers
    int nb_ob = round(float(nObs/ntpb_y) + 0.5f);
    if (nb_ob>nby_max) nb_ob = nby_max;
    int ntpb_ob = ntpb_y;
    if (nObs<ntpb_ob){
        ntpb_ob = nObs;
        nb_ob = 1;
    }
    // elements
    int nb_el = nbz_max;
    int ntpb_el = ntpb_z;


    printf("\n*** GPU calculation ...\n");
   
    printf("\n      BLOCK/THREAD splitting \n\n");
     
    printf("          Time steps: %d / %d\n",nb_t,ntpb_t);
    printf("          Observers:  %d / %d\n",nb_ob,ntpb_ob);
    printf("          Elements:   %d / %d\n\n",nb_el,ntpb_el);

    dim3 blocksPerGrid(nb_t,nb_ob,nb_el);
    dim3 threadsPerBlock(ntpb_t,ntpb_ob,ntpb_el);

    
    if(impermeable==1) findPressureTimeHistory<<<blocksPerGrid,threadsPerBlock>>>(rObs_gpu,tObs_gpu,pT_gpu,pL_gpu,pA_gpu,Xham_gpu,Nham_gpu,Vham_gpu,PRham_gpu,nSurf,XPtr_g,NPtr_g,VPtr_g,PRPtr_g,E_g,nObs,nTime,dPsi,dTau,a0,rhoRef,vInf_gpu,oM);
   
    cudaEventRecord(stop5, 0);
    cudaEventSynchronize(stop5);
    cudaEventElapsedTime(&elapsedTime, stop4, stop5);
    printf("      GPU time: %lf ms\n\n",elapsedTime);


    // copy solutions from device to host
    cudaMemcpy(pT,pT_gpu,size_tobs,cudaMemcpyDeviceToHost);
    cudaMemcpy(pL,pL_gpu,size_tobs,cudaMemcpyDeviceToHost);
    cudaMemcpy(pA,pA_gpu,size_tobs,cudaMemcpyDeviceToHost);
   



     
    printf("*** Writing output files\n\n");
    writeTimeHistory(pT,pL,pA);
    
    if(BB_noise==1) writeBBoutput(BBoaspl,BBspl);  
 
   // Deallocate memory
    cudaFree(Xham_gpu);cudaFree(Nham_gpu);cudaFree(Vham_gpu);cudaFree(PRham_gpu);
    cudaFree(pT_gpu);cudaFree(pL_gpu);cudaFree(pA_gpu);cudaFree(rObs_gpu);cudaFree(tObs_gpu);cudaFree(vInf_gpu);
    cudaFree(XPtr_g);cudaFree(NPtr_g);cudaFree(VPtr_g);cudaFree(PRPtr_g);cudaFree(E_g);
    //if(impermeable==0) cudaFree(Uham_gpu);cudaFree(UPtr_g);
    
    free(Xham);free(Nham);free(Vham);free(Uham);free(PRham);
    free(pT);free(pL);free(pA);free(tObs);
    free(XPtr);free(NPtr);free(VPtr);free(UPtr);free(PRPtr);free(E);

//    if(BB_noise==1) 
//    {
//        free(BBspl);free(BBoaspl);free(span);free(omega);free(ccw);free(psi_offset);free(trip);free(rotate);free(translate);free(xsPtr);free(dsPtr);free(nsect);free(xyzcR);free(bbdata);free(nsect_t);free(nmach_t);free(naoa_t);free(nre_t);free(aPtr);free(mPtr);free(rPtr);free(sPtr);free(dPtr);free(sect_t);free(mach_t);free(aoa_t);free(re_t);free(bbdata_t);free(ff);
//    }


    // Stop time
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    printf("\n*** Total time: %lf ms\n\n",elapsedTime);




    printf("    DONE.....!\n");


}// close main


