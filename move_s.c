/*MOVE PHOTON*/
void move_s(double x[], double u[], double d, double *depth_pp, double *path_pp, double *s, double *r_tissue, double *num_pp) {
    double xp0,xp1,xp2;
    double v2;
    xp0=x[0]+u[0]*d;
    xp1=x[1]+u[1]*d;
    xp2=x[2]+u[2]*d;
    x[0]=xp0;
    x[1]=xp1;   
    x[2]=xp2;
    
    v2 = d*sqrt(pow(u[0],2)+pow(u[1],2)+pow(u[2],2));
    
    *s-=d;
    *path_pp += v2;
    *num_pp += 1;
    *r_tissue = sqrt(pow((x[0]),2)+pow(x[1],2));          /*Position updated after move*/ 
    if (xp2>*depth_pp){*depth_pp=xp2;}                    /*maximum depth*/
    return;
}