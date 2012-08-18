/*============================*
 * Monte Carlo for one medium *
 *============================*/

#define PI 3.14159265
#include "mex.h"
#include "blas.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>
#define SIGN(x) ((x)>=0 ? 1:-1)
#include "ran1.c"                /*RANDOM NUMBER GENERATOR*/
#include "specular.c"            /*SPECULAR REFLECTION*/
#include "move_s.c"                /*MOVE PHOTON*/
#include "absorption.c"          /*ABSORPTION*/
#include "spin_thph.c"                /*SCATTERING*/
#include "roulette.c"            /*ROULETTE*/
#include "reflection.c"          /*INTERNAL REFLECTION*/
/*MONTE CARLO*/
void mc1m(int num_photons,double ma, double ms, double g, double *x0,double *x1, double *x2, double *u0, double *u1, double *u2,
        double *na_fiber, double r_max, double x2_max,double n_tissue, double *n_fiber, double *fiber_radii, int t_grid,
        double *th, double *ph,double tissue[],double  fiber_tissue[], double r[],double depth[],double weight[],double path[],
        double z[],double num_scatt[]) {
      
    int     i;
    double  s = 0;
    /*Fiber NA*/
    float   na_clad = na_fiber[1];    
    float   na_core = na_fiber[0];    
    /*Fiber refractive indices*/
    float   n_clad = n_fiber[1];    
    float   n_core = n_fiber[0]; 
    /*Maximum angle for detection*/
    double  theta_clad = asin(na_clad/n_tissue);
    double  Ralphi;
    /*ran1 seed*/
    srand ( time(NULL) );       
    long idum = -(rand() % 1000 +1);
    
    /*Image array*/
    int     t_grid0 = t_grid*r_max;
    int     t_grid2 = t_grid*x2_max;
    double  dim_arr = t_grid0*t_grid2;
    double  *aux_tissue = malloc(dim_arr*sizeof(double));
    double  r_tissue;
    /*Fibre geometry*/
    double  r_dclad = fiber_radii[2];  
    double  r_clad = fiber_radii[1];
    double  r_core = fiber_radii[0];  
    
    /*Tissue and fiber_tissue initialization for each monte carlo run*/
    int it,jt, count=0;                 
    for (it = 0; it< t_grid2; it++){            
       for (jt = 0; jt< t_grid0; jt++){
            *(tissue+count) = 0;
            *(fiber_tissue+count) = 0;
            count++;}}
    /*Auxiliar arrays initialization*/
    int ji,count2=0;
    for (ji = 0; ji< num_photons; ji++){    
    	*(r+count2) = 0;
     	*(depth+count2) = 0;
    	*(weight+count2) = 0;
     	*(path+count2) = 0;
    	*(z+count2) = 0;
        *(num_scatt+count2)=0;
     	count2++;}      
    
    /*LAUNCHING*/                 
    for (i=0;i<num_photons;i++) {
        int status=1;                   /*Photon alive*/
        double x[]={x0[i],x1[i],x2[i]}; /*Position*/        
        double u[]={u0[i],u1[i],u2[i]}; /*Direction*/
        double path_pp=0;               /*Pathlength*/
        double weight_pp=1;             /*Weigth*/
        double num_pp = 0;
        double depth_pp=0;              /*Maximum depth*/
        double dw=0;                    /*Weigth lost*/
        double fiber_boundary;
        double z_na, r_na, cone;
        
        /*Specular reflection*/
        specular(n_core, n_tissue,&weight_pp); 
       
        /*aux_ tissue initialization for each photon*/
        int count1 = 0,j;
        for (it = 0; it< t_grid2; it++){       
            for (jt = 0; jt< t_grid0; jt++){
            *(aux_tissue+count1)=0;
            count1++;}}
      
        while (status==1) {             
            /*Propagation distance*/
            if (s == 0){s=log(ran1(&idum))/(-(ma+ms));}
            int tz0=0;
            int tz2=0;
                        
            /*HIT FIBER BOUNDARY*/            
            fiber_boundary = -x[2]*(ma+ms)/u[2];
            /*Position before updating*/
            z_na = x[2];                   
            r_na = r_tissue;
            cone = theta_clad*z_na + r_clad;
            
            if (u[2] < 0  && fiber_boundary <= s && z_na >0) {     
               move_s(x, u, fiber_boundary, &depth_pp, &path_pp, &s, &r_tissue,&num_pp);  
               /*whether the photon is internally reflected*/ 
               reflection(u, n_core, n_tissue, &Ralphi);               
               if (ran1(&idum)<= Ralphi){u[2]=-u[2];}   
            }      
            else{                                              
                move_s(x, u, s, &depth_pp, &path_pp, &s, &r_tissue, &num_pp);
                absorption(ma, ms,&weight_pp,&dw);
                spin_thph(u, g, &idum,th,ph);            
               
                /*Absortion array updating*/
                tz0=floor(r_tissue*t_grid);
                tz2=floor(x[2]*t_grid);
                	if(tz0>t_grid0-1){tz0=t_grid0-1;}
                	if(tz0<0){tz0=0;}
                	if(tz2>t_grid2-1){tz2=t_grid2-1;}
                	if(tz2<0){tz2=0;}
                    *(tissue+(tz0*t_grid2+tz2))+=dw;
                	*(aux_tissue+(tz0*t_grid2+tz2))+=dw;                   
             } /*else (the photon didn't reach the boundary)*/
            
            /*PHOTON TERMINATION*/
            if(x[2] > x2_max ){status=0;}
            if(status==1 &&  weight_pp<0.00001 /*threshold*/ ){roulette(&weight_pp,&idum,&status);}                
       
       /*The photon reach the sensing area?*/
       if(x[2] < 0 && z_na >0){
            for(j=0;j<dim_arr;j++){fiber_tissue[j]+=aux_tissue[j];}
            /*Final status of detected photons*/
            r[i]      = r_tissue;
            z[i]      = x[2];
            path[i]   = path_pp;
            weight[i] = weight_pp;
            depth[i]  = depth_pp;    
            num_scatt[i] = num_pp;}
                        
         } /*while status is alive*/   
     } /*for i num_photons launched*/
    return;
} /*mcml2m*/


/*MEX FUNCTION*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const int MAX_PHOTONS=10000000;
    int num_photons,t_grid;
    double ma, ms, g, r_max,x2_max,n_tissue;      
    double *nphot,*mua,*grid,*mus,*gs, *nt, *x0,*x1,*x2,*u0,*u1,*u2,*na_fiber,*rmax,*x2max,*r,*depth,*weight,
            *path,*z,*num_scatt,*n_fiber,*fiber_radii,*th,*ph; /* *th_tofiber pointers to input & output matrices*/
    double* tissue;
    double* fiber_tissue;
    
    /*Inputs*/
    nphot = mxGetPr(prhs[0]);    
    mua = mxGetPr(prhs[1]);
    mus = mxGetPr(prhs[2]);
    gs = mxGetPr(prhs[3]);
    x0 = mxGetPr(prhs[4]);
    x1 = mxGetPr(prhs[5]);
    x2 = mxGetPr(prhs[6]);
    u0 = mxGetPr(prhs[7]);
    u1 = mxGetPr(prhs[8]);
    u2 = mxGetPr(prhs[9]);
    na_fiber = mxGetPr(prhs[10]);
    rmax = mxGetPr(prhs[11]);
    x2max = mxGetPr(prhs[12]);
    nt = mxGetPr(prhs[13]);    
    n_fiber = mxGetPr(prhs[14]);
    fiber_radii = mxGetPr(prhs[15]);    
    grid =mxGetPr(prhs[16]);
    th = mxGetPr(prhs[17]);
    ph =mxGetPr(prhs[18]);
    /*********/
    num_photons = (int)nphot[0];     
    ma = (double)mua[0]; 
    ms = (double)mus[0]; 
    g = (double)gs[0]; 
    r_max = (double)rmax[0];
    x2_max = (double)x2max[0];
    n_tissue = (double)nt[0];
    t_grid = (int)grid[0];
            
    if (num_photons > MAX_PHOTONS ) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                "Number of photons exceeds limit.");
    }    
    if (num_photons < 1 ) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                "Too few photons.");
    }
     
    mwSignedIndex  dims[2];
    dims[0]=t_grid*x2_max;
    dims[1]=t_grid*r_max;
    
    /*Outouts*/
    plhs[0] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    tissue = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS,mxREAL);
    fiber_tissue = mxGetPr(plhs[1]);        
    plhs[2] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    r = mxGetPr(plhs[2]);    
    plhs[3] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    depth = mxGetPr(plhs[3]);    
    plhs[4] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    weight = mxGetPr(plhs[4]);    
    plhs[5] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    path = mxGetPr(plhs[5]);    
    plhs[6] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    z = mxGetPr(plhs[6]);        
    plhs[7] = mxCreateDoubleMatrix(num_photons, 1, mxREAL);
    num_scatt = mxGetPr(plhs[7]);
    /*********/
    
    /*run C code*/
    mc1m(num_photons,ma,ms,g,x0,x1,x2,u0,u1,u2,na_fiber,r_max,x2_max,n_tissue,n_fiber,fiber_radii,
            t_grid, th, ph,tissue, fiber_tissue,r,depth,weight, path, z,num_scatt);
  return;
}
