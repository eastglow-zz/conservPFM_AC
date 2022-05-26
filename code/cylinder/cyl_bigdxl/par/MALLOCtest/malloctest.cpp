#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int* makeINT_1D(int size);
double* makeDOUBLE_1D(int size);
double** makeDOUBLE_2D(int size1,int size2);
double*** makeDOUBLE_3D(int size1,int size2,int size3);
void erase1D(void *p1);
void eraseDOUBLE_2D(double **p2,int size1,int size2);
void eraseDOUBLE_3D(double ***p3,int size1,int size2,int size3);
int main()
{
    ///////////////
    double *p1;
    double **p2;
    double ***p3;
    double ****p4;
    int np,nx,ny,nz;
    int nn,i,j,k;
    ///////////////
    np = 3;
    nx = 2;
    ny = 2;
    nz = 2;

    p1 = makeDOUBLE_1D(np);
    for (nn=0;nn<np;nn++) {
        p1[nn]=nn+1;
        printf("%d %lf %ld\n",nn,p1[nn],&p1[nn]);
    }
    erase1D(p1);

    p2 = makeDOUBLE_2D(np,nx);
    for (nn=0;nn<np;nn++) {
        for (i=0;i<nx;i++) {
            p2[nn][i] = i;
            printf("p2[%d][%d]=%lf,%ld\n",nn,i,p2[nn][i],&p2[nn][i]);
        }
    }
    eraseDOUBLE_2D(p2,np,nx);

    p3 = makeDOUBLE_3D(np,nx,ny);
    for (nn=0;nn<np;nn++) {
        for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            p3[nn][i][j]= j;
            printf("p3[%d][%d][%d]=%lf,%ld\n",nn,i,j,p3[nn][i][j],&p3[nn][i][j]);
        }
        }
    }
    eraseDOUBLE_3D(p3,np,nx,ny);
    return 0;
}

int* makeINT_1D(int size)
{
    int *p1;
    p1 = (int *)calloc(size,sizeof(int));
    return p1;
}

double* makeDOUBLE_1D(int size)
{
    double *p1;
    p1 = (double *)calloc(size,sizeof(double));
    return p1;
}

void erase1D(void *p1)
{
    free(p1);
    p1=NULL;
}

double** makeDOUBLE_2D(int size1,int size2)
{
    double **p2;
    int i;
    p2 = (double **)calloc(size1,sizeof(double));
    p2[0] = (double *)calloc(size1*size2,sizeof(double));
    for (i=1;i<size1;i++) {
        p2[i] = p2[0] + i*size2;
    }
    return p2;
}

double*** makeDOUBLE_3D(int size1,int size2,int size3)
{
    double ***p3;
    int i,j;
    p3 = (double ***)calloc(size1,sizeof(double));
    p3[0] = (double **)calloc(size1*size2,sizeof(double));
    p3[0][0] = (double *)calloc(size1*size2*size3,sizeof(double));
    for (i=1;i<size1;i++) {
        p3[i]=p3[0] + i*size2;
    }
    for (i=0;i<size1;i++) {
        for (j=1;j<size2;j++) { 
            p3[i][j] = p3[i][0] + j*size3;
        }
    }
 /*
    p3[0][0] = (double *)calloc(size1*size2*size3,sizeof(double));
    for (i=0;i<size1;i++) {
        for (j=0;j<size2;j++) {
            p3[i][j] = (p3[0]+i*size1)+j*size2;
        }
    }
 */
    return p3;
}

void eraseDOUBLE_2D(double **p2,int size1,int size2)
{
    free(p2[0]);
    free(p2);
    p2=NULL;
}


void eraseDOUBLE_3D(double ***p3,int size1,int size2,int size3)
{
    free(p3[0][0]);
    free(p3[0]);
    free(p3);
    p3=NULL;
}
