#include <stdio.h>
#include <stdlib.h>

#define nx 5
#define ny 6
#define nz 7

int main()
{
    int i,j,k;
    int index;
    int is,js,ks;
    int nxn=0,nyn=0,nzn=0;
    index = 0;
    for (i=0;i<nx;i++) {
        if (i==nx-1)nxn++;
        for (j=0;j<ny;j++) {
            if (j==ny-1)nyn++;
            for (k=0;k<nz;k++) {
                if (k==nz-1)nzn++;
                printf("%d %d %d %d, %d %d %d\n",i,j,k,index,nxn,nyn,nzn);
                index++;
            }
        }
    }

    return 0;
}

