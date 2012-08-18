 /*SURFACE REFLECTION*/
 void reflection(double u[],double n_trans, double n_inc, double *Ralphi){
     double alphi,alpht;
        alphi=acos(abs(u[2]));
        alpht=asin((n_inc/n_trans)*sin(alphi));
     if (alphi > asin(n_trans/n_inc)){*Ralphi=1;}
     else {*Ralphi=0.5*(pow(sin(alphi-alpht),2)/pow(sin(alphi+alpht),2)+pow(tan(alphi-alpht),2)/(pow(tan(alphi+alpht),2)));}
 }
 
 
