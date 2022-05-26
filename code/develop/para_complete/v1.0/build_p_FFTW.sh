#bin/sh

mpicc -c ./CODES/useFFTW.cpp -O2 -lm -lfftw3_mpi -lfftw3
ar rcv ./CODES/libFFTW3.a ./useFFTW.o
rm -f ./useFFTW.o

mpicc -c ./CODES/main_p_FFTW.cpp -O2 -lm -lfftw3_mpi -lfftw3
mpicc -c ./CODES/common_functions.cpp -O2 -lm

mpicc -o set65pFFTW.exe ./*.o ./CODES/libFFTW3.a -O2 -lm -lfftw3_mpi -lfftw3
rm *.o
