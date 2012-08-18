
/**random number generator**/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2E-07
#define RNMX 1.0-EPS

float ran1(long *idum) {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idum <= 0 || !iy) { /*Initialize*/
        if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0*/
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) { /*Load the shue table(after 8 warm-ups)*/
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

#undef IA
#undef IM
#undef IQ
#undef IR
#undef NTAB
#undef EPS
#undef MAX
#undef MIN
/****************************/

