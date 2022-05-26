#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define nxtot 52
#define nytot 52
#define nztot 52
#define nx 48
#define ny 48
#define nz 48

int shrink_buffer(int *ista, int*iend,int n1, int n2, int nprocs, int myrank);
void getRANK_NPROCS(int *rank, int *size);

const int nORI = nz;

int nzshr = 1;
int RANK, NPROCS;
int kmin,kmax;

int main(int argc, char *argv[])
{
    ///////////////////
    int *p;
    int i;
    int ical;
    int gosignal=1;
    MPI_Status status1;
    ///////////////////
    MPI_Init(&argc,&argv);
    getRANK_NPROCS(&RANK,&NPROCS);
    nzshr=shrink_buffer(&kmin,&kmax,1,nORI,NPROCS,RANK);
    nzshr+=4; //for boundary condition
    ///////////////////
    p = new int[nxtot*nytot*nzshr]; 
    for (i=0;i<nzshr;i++) {
        ical = RANK*nzshr + i;
        p[i]=ical;
    }

    if (RANK!=0) {
        printf("before Recv, RANK=%d\n",RANK);
        MPI_Recv(&gosignal,1,MPI_INT,RANK-1,gosignal,MPI_COMM_WORLD,&status1);
        printf("after Recv, RANK=%d\n",RANK);
    }
    for (i=0;i<nzshr;i++) {
        printf("p[%d]=%d\n",kmin+i,p[i]);
    }
    if (RANK!=NPROCS-1) {
        printf("before Ssend, RANK=%d\n",RANK);
        MPI_Ssend(&gosignal,1,MPI_INT,RANK+1,gosignal,MPI_COMM_WORLD);
        printf("after Ssend, RANK=%d\n",RANK);
    }

    MPI_Finalize();

    return 0;
}

int shrink_buffer(int *ista, int*iend,int n1, int n2, int nprocs, int myrank)
{
    int iwork1, iwork2;
    iwork1 = (n2-n1+1)/nprocs;
    iwork2 = (n2-n1+1)%nprocs;
    *ista = myrank*iwork1 + n1 + (myrank<iwork2?myrank:iwork2);
    *iend = *ista + iwork1 - 1;
    if (iwork2>myrank) *iend = *iend + 1;
    return *iend - *ista + 1;
}
void getRANK_NPROCS(int *rank, int *size)
{
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD,size);
}
