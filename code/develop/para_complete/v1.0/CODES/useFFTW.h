#define FFTW_REAL 0
#define FFTW_IMAG 1

void fourier_ready(fftw_plan *planRtoK, fftw_plan *planKtoR, fftw_complex *data, int nx, int ny, int nz);
void fourier_end(fftw_plan *planRtoK, fftw_plan *planKtoR);
void fourierRtoK(fftw_complex *data, fftw_plan *planRtoK);
void fourierKtoR(fftw_complex *data, fftw_plan *planKtoR);
void del_buff_fcx(fftw_complex **p);
