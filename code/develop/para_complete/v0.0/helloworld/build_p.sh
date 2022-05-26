mpicc -c ./*.cpp -lm
mpicc -o first.exe ./*.o -lm
rm *.o
