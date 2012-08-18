/*SCATTERING*/
void spin_thph(double u[], double g,long *idum, double th[], double ph[]) {
    double cos_theta,sin_theta, psi,cos_psi,sin_psi,theta;
    double up0,up1,up2, eme;
    int idx;
    /*Henyey-Greenstein
    if(g==0){
        cos_theta=2*ran1(idum)-1;
    }
    else {
        cos_theta=(1/(2*g))*(1+pow(g,2)-pow(((1-pow(g,2))/(1-g+(2*g*ran1(idum)))),2)); 
    }
        sin_theta=sin(acos(cos_theta));
    */    
    
    
    /*Toublanc D. Henyey-greenstein and mie phase functions in monte carlo radiative transfer computations. Appl Opt. 1996;35(18):3270-4.*/
    eme = ran1(idum); 
    idx = -1;
    do{idx = idx++;}
    while(ph[idx] < eme);
    theta = th[idx];
    /***********************************************************/
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    psi=2*PI*ran1(idum);    
    cos_psi=cos(psi);
    sin_psi=sin(psi);
  
      if (abs(u[2])>0.99999){
      	up0=sin_theta*cos_psi;
      	up1=sin_theta*sin_psi;
      	up2=SIGN(u[2])*cos_theta;
      }
      else{
      	up0=(sin_theta/sqrt(1-pow(u[2],2)))*(u[0]*u[2]*cos_psi-u[1]*sin_psi)+u[0]*cos_theta;        
      	up1=(sin_theta/sqrt(1-pow(u[2],2)))*(u[1]*u[2]*cos_psi+u[0]*sin_psi)+u[1]*cos_theta;
        up2=-sin_theta*cos_psi*sqrt(1-pow(u[2],2))+u[2]*cos_theta;
      }      
    u[0]=up0;
    u[1]=up1;
    u[2]=up2;
    return;
}
