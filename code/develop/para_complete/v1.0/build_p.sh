mpicc -c ./*.cpp -lm
mpicc -o set65p.exe ./*.o -lm
rm *.o
