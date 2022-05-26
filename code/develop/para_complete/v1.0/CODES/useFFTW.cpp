#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <fftw3-mpi.h>


void fourier_ready(fftw_plan *planRtoK, fftw_plan *planKtoR, fftw_complex *data, int nx, int ny, int nz)
{   
    fftw_mpi_init();
    *planRtoK = fftw_mpi_plan_dft_3d(nx,ny,nz,data,data,MPI_COMM_WORLD,FFTW_FORWARD, FFTW_ESTIMATE);
    *planKtoR = fftw_mpi_plan_dft_3d(nx,ny,nz,data,data,MPI_COMM_WORLD,FFTW_BACKWARD, FFTW_ESTIMATE);
    
}

void fourier_end(fftw_plan *planRtoK, fftw_plan *planKtoR)
{
   fftw_destroy_plan(*planRtoK);
   fftw_destroy_plan(*planKtoR);
}

void fourierRtoK(fftw_complex *data, fftw_plan *planRtoK)
{
    fftw_execute(*planRtoK);
}

void fourierKtoR(fftw_complex *data, fftw_plan *planKtoR)
{
    //////////////////////
    int datsize;
    double complex fnorm;
    //////////////////////
    datsize = sizeof(data)/sizeof(fftw_complex);
    fnorm = 1./((double complex) datsize);
    
    fftw_execute(*planKtoR);
    for (int i=0;i<datsize;i++) {
        data[i]*=fnorm;
    }
}

void del_buff_fcx(fftw_complex **p)
{
    delete[] *p;
    *p = NULL;
}
