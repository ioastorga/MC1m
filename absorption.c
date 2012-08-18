/*ABSORPTION*/
void absorption(double ma, double ms, double *weight_pp, double *dw)
{
    *dw =  (ma/(ma+ms))* (*weight_pp);
    *weight_pp-=*dw; 
    return;
}