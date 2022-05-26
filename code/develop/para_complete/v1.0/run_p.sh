#mpirun -machinefile ./mf -np 8 ./set65p.exe
nohup mpirun -machinefile ./mf -np 64 ./set65pFFTW.exe &
#mpirun -machinefile ./mf -np 64 ./set65pFFTW.exe
