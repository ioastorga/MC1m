 /*SPECULAR REFLECTION*/
 void specular(float na, float n1, double *weight_pp){
 double Rsp;
 Rsp=pow((na-n1),2)/pow((na+n1),2); 
 *weight_pp-=Rsp;
 /*decrease of the photon to enter the medium*/
 }