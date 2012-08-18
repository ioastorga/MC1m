/*ROULETTE*/
void roulette(double *weight_pp,long *idum,int *status){
            if(ran1(idum)<=0.1/*chance*/){ *weight_pp/=0.1; }
            else{*status=0;}
    return;
}