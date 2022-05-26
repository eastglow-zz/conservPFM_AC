#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

int main()
{   
    //////////////////
    double complex a[10]={1.,2.,3.,4.,5.,6.,7.,8.,9+1.*I,10*I};
    int size;
    //////////////////
    size = sizeof(a)/sizeof(double complex);
    double complex *b;
    b = new double complex[size];
    printf("array size = %d\n",size);
    for (int i=0;i<size;i++) {
        b[i] = a[i]/((double complex) size);
        printf("a[%d] = %lf+%lfi, b[%d]= %lf+%lfi\n",i,creal(a[i]),cimag(a[i]),i,creal(b[i]),cimag(b[i]));
    }

    return 0;
}
