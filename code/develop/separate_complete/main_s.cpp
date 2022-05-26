#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defined_values.h"
#include "common_functions.h"

int nx = nxcal;
int ny = nycal;
int nz = nzcal;
int ntot = nxtot*nytot*nztot;

int nmax = max3(nx,ny,nz);

extern double time, dt, dxl, dyl, dzl;
extern int itime;

int main()
{
    /////////////////////
    void fileout_g(void);
    void fileout_dg(void);
    void initialize(double *p[nphase]);
    void BC_periodic3D(double *p);
    void BC_adiabatic3D(double *p);
    void BC_userdefine3D(double *p);
    void calc_lamb(double *pnew[nphase], double *pold[nphase]);
    void fileout_p(double *p[nphase]);
    double calc_frac_p(double *p);
    void frac_out(double *p[nphase]);
    void clean_RR_HH(void);
    /////////////////////
    double *p[2][nphase]; //phase field old/new switchable
    double *pinter[nphase];
    for (int np=0;np<nphase;np++) {
        p[OLD(itime)][np]=new double[ntot];
        p[NEW(itime)][np]=new double[ntot];
        pinter[np]=new double[ntot];
    }
    /////////////////////
    //fileout_g();
    //fileout_dg();
        
    initialize(p[OLD(0)]);
    //temporal iteration
    for (itime=0, time=0.; 1 ;itime++,time+=dt) {
        //Applying boundary condition on p
        {
            int np;
            for (np=0;np<nphase;np++) {
                //BC_periodic3D(p[OLD(itime)][np]);
                BC_adiabatic3D(p[OLD(itime)][np]);
                //BC_userdefine3D(p[OLD(itime)][np]);
            }
        }

        //Solving for PF eqn.; FTCS scheme
        //solvPF_CH(p[NEW(itime)],p[OLD(itime)]);
        //solvPF_CA(p[NEW(itime)],p[OLD(itime)]);

        if(1){
            //clean_RR_HH();
            solvPF_CA(pinter,p[OLD(itime)]);
            int np;
            for (np=0;np<nphase;np++) {
                //BC_periodic3D(p[np]);
                BC_adiabatic3D(pinter[np]);
               //BC_userdefine3D(p[np]);
            }
            calc_lamb(pinter,p[OLD(itime)]);
            drivPF_CA(p[NEW(itime)],pinter);
        } 


		
        if (itime%100  ==0) {
            fileout_p(p[OLD(itime)]);
            frac_out(p[OLD(itime)]);
            //system("pause");
        }

    }

    for (int np=0;np<nphase;np++) {
        del_buff_d(&p[OLD(itime)][np]);
        del_buff_d(&p[NEW(itime)][np]);
        del_buff_d(&pinter[np]);
    }
    //system("pause");
    return 0;
}

void initialize(double *p[nphase])
{
        ///////////
        void clean_all(double *p[nphase]);
        void clean_one(double *p);
        void relax_p(double *p[nphase], int nrelax);
        void copy_p(double *presult[nphase], double *psource[nphase]);
        void draw_cube(double *p,int ix, int iy, int iz, int size);
        void draw_sphere(double *p,int ix, int iy, int iz, int diameter);
        void draw_cylinder(double *p,int ix, int iy, int iz, char axis, int rad, int height);
        void draw_rectangular(double *p,int ix, int iy, int iz, int sizex, int sizey,int sizez);
        void fill_rectangular(int phaseindex,double *p[nphase],int ix, int iy, int iz, int sizex, int sizey, int sizez);
        void fill_cylinder(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int rad, int height);
        void fill_environ(int phaseindex,double *p[nphase]);
	///////////
	clean_all(p); 
        //draw_cylinder(p[2],nx/2,ny/2,nz/2+6,'z',16,10);
        //fill_cylinder(1,p,nx/2,ny/2,nz/2-6,'z',16,10);
        draw_rectangular(p[1],nx/2-24,ny/2,1 ,48,64,1);
        draw_rectangular(p[2],nx/2+24,ny/2,1 ,48,64,1);
        fill_environ(0,p);
        relax_p(p, 50);
        //clean_one(p[0]);
        //fill_cylinder(2,p,nx/2,ny/2,10+10,'z',20,20);
        //fill_environ(0,p);
}


void BC_adiabatic3D(double *p)
{
	//////////////////
	int l,m,n;
	int lim;
	//////////////////
	for (l=1;l<=nmax;l++) {
	    //layer
	    for (m=1;m<=nmax;m++) {
	        //i-layer (j,k varying)
                if (l<=ny && m <= nz) {
                    p[ic(ii(1)-2,ii(l),ii(m))]=p[ic(ii(2),ii(l),ii(m))];
                    p[ic(ii(1)-1,ii(l),ii(m))]=p[ic(ii(1),ii(l),ii(m))];
                    p[ic(ii(nx)+1,ii(l),ii(m))]=p[ic(ii(nx),ii(l),ii(m))];
                    p[ic(ii(nx)+2,ii(l),ii(m))]=p[ic(ii(nx-1),ii(l),ii(m))];
                }
                //j-layer (k,i varying)
                if (l<=nz && m <= nx) {
                    p[ic(ii(m),ii(1)-2,ii(l))]=p[ic(ii(m),ii(2),ii(l))];
                    p[ic(ii(m),ii(1)-1,ii(l))]=p[ic(ii(m),ii(1),ii(l))];
                    p[ic(ii(m),ii(ny)+1,ii(l))]=p[ic(ii(m),ii(ny),ii(l))];
                    p[ic(ii(m),ii(ny)+2,ii(l))]=p[ic(ii(m),ii(ny-1),ii(l))];
                }
                //k-layer (i,j varying)
                if (l<=nx && m <= ny) {
                    p[ic(ii(l),ii(m),ii(1)-2)]=p[ic(ii(l),ii(m),ii(2))];
                    p[ic(ii(l),ii(m),ii(1)-1)]=p[ic(ii(l),ii(m),ii(1))];
                    p[ic(ii(l),ii(m),ii(nz)+1)]=p[ic(ii(l),ii(m),ii(nz))];
                    p[ic(ii(l),ii(m),ii(nz)+2)]=p[ic(ii(l),ii(m),ii(nz-1))];
                }
            }
	    //edge
	    //i-direction (i varying)
	    if (l<=nx) {
	    	p[ic(ii(l),ii(1)-1,ii(1)-1)]=p[ic(ii(l),ii(1),ii(1))];
		p[ic(ii(l),ii(1)-1,ii(nz)+1)]=p[ic(ii(l),ii(1),ii(nz))];
		p[ic(ii(l),ii(ny)+1,ii(1)-1)]=p[ic(ii(l),ii(ny),ii(1))];
		p[ic(ii(l),ii(ny)+1,ii(nz)+1)]=p[ic(ii(l),ii(ny),ii(nz))];
	    }
	    //j-direction (j varying)
	    if (l<=ny) {
		p[ic(ii(1)-1,ii(l),ii(1)-1)]=p[ic(ii(1),ii(l),ii(1))];
		p[ic(ii(1)-1,ii(l),ii(nz)+1)]=p[ic(ii(1),ii(l),ii(nz))];
		p[ic(ii(nx)+1,ii(l),ii(1)-1)]=p[ic(ii(nx),ii(l),ii(1))];
		p[ic(ii(nx)+1,ii(l),ii(nz)+1)]=p[ic(ii(nx),ii(l),ii(nz))];
	    }
	    //k-direction (k varying)
	    if (l<=nz) {
	 	p[ic(ii(1)-1,ii(1)-1,ii(l))]=p[ic(ii(1),ii(1),ii(l))];
		p[ic(ii(1)-1,ii(ny)+1,ii(l))]=p[ic(ii(1),ii(ny),ii(l))];
		p[ic(ii(nx)+1,ii(1)-1,ii(l))]=p[ic(ii(nx),ii(1),ii(l))];
		p[ic(ii(nx)+1,ii(ny)+1,ii(l))]=p[ic(ii(nx),ii(ny),ii(l))];
	    }
        }
	p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1))];
	p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1))];
	p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1))];
	p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1))];

	p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz))];
	p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz))];
	p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz))];
	p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz))];
}

void BC_periodic3D(double *p)
{
	//////////////////
	int l,m,n;
	int lim;
	//////////////////
	
	for (l=1;l<=nmax;l++) {
		//layer
		for (m=1;m<=nmax;m++) {
			//i-layer (j,k varying)
			if (l<=ny && m <= nz) {
				p[ic(ii(1)-2,ii(l),ii(m))]=p[ic(ii(nx-1),ii(l),ii(m))];
				p[ic(ii(1)-1,ii(l),ii(m))]=p[ic(ii(nx),ii(l),ii(m))];
				p[ic(ii(nx)+1,ii(l),ii(m))]=p[ic(ii(1),ii(l),ii(m))];
				p[ic(ii(nx)+2,ii(l),ii(m))]=p[ic(ii(2),ii(l),ii(m))];
			}
			//j-layer (k,i varying)
			if (l<=nz && m <= nx) {
				p[ic(ii(m),ii(1)-2,ii(l))]=p[ic(ii(m),ii(ny-1),ii(l))];
				p[ic(ii(m),ii(1)-1,ii(l))]=p[ic(ii(m),ii(ny),ii(l))];
				p[ic(ii(m),ii(ny)+1,ii(l))]=p[ic(ii(m),ii(1),ii(l))];
				p[ic(ii(m),ii(ny)+2,ii(l))]=p[ic(ii(m),ii(2),ii(l))];
			}
			//k-layer (i,j varying)
			if (l<=nx && m <= ny) {
				p[ic(ii(l),ii(m),ii(1)-2)]=p[ic(ii(l),ii(m),ii(nz-1))];
				p[ic(ii(l),ii(m),ii(1)-1)]=p[ic(ii(l),ii(m),ii(nz))];
				p[ic(ii(l),ii(m),ii(nz)+1)]=p[ic(ii(l),ii(m),ii(1))];
				p[ic(ii(l),ii(m),ii(nz)+2)]=p[ic(ii(l),ii(m),ii(2))];
			}
		}
		//edge
		//i-direction (i varying)
		if (l<=nx) {
			p[ic(ii(l),ii(1)-1,ii(1)-1)]=p[ic(ii(l),ii(ny),ii(nz))];
			p[ic(ii(l),ii(1)-1,ii(nz)+1)]=p[ic(ii(l),ii(ny),ii(1))];
			p[ic(ii(l),ii(ny)+1,ii(1)-1)]=p[ic(ii(l),ii(1),ii(nz))];
			p[ic(ii(l),ii(ny)+1,ii(nz)+1)]=p[ic(ii(l),ii(1),ii(1))];
		}
		//j-direction (j varying)
		if (l<=ny) {
			p[ic(ii(1)-1,ii(l),ii(1)-1)]=p[ic(ii(nx),ii(l),ii(nz))];
			p[ic(ii(1)-1,ii(l),ii(nz)+1)]=p[ic(ii(nx),ii(l),ii(1))];
			p[ic(ii(nx)+1,ii(l),ii(1)-1)]=p[ic(ii(1),ii(l),ii(nz))];
			p[ic(ii(nx)+1,ii(l),ii(nz)+1)]=p[ic(ii(1),ii(l),ii(1))];
		}
		//k-direction (k varying)
		if (l<=nz) {
			p[ic(ii(1)-1,ii(1)-1,ii(l))]=p[ic(ii(nx),ii(ny),ii(l))];
			p[ic(ii(1)-1,ii(ny)+1,ii(l))]=p[ic(ii(nx),ii(1),ii(l))];
			p[ic(ii(nx)+1,ii(1)-1,ii(l))]=p[ic(ii(1),ii(ny),ii(l))];
			p[ic(ii(nx)+1,ii(ny)+1,ii(l))]=p[ic(ii(1),ii(1),ii(l))];
		}
	}
	//coners
	p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1)-1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1)-1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1)-1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1)-1)];

	p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz)+1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz)+1)];
}

void BC_userdefine3D(double *p)
{
        //i,j: periodic
        //k: adiabatic
	//////////////////
	int l,m,n;
	int lim;
	//////////////////
	
	for (l=1;l<=nmax;l++) {
		//layer
		for (m=1;m<=nmax;m++) {
			//i-layer (j,k varying)
			if (l<=ny && m <= nz) {
				p[ic(ii(1)-2,ii(l),ii(m))]=p[ic(ii(nx-1),ii(l),ii(m))];
				p[ic(ii(1)-1,ii(l),ii(m))]=p[ic(ii(nx),ii(l),ii(m))];
				p[ic(ii(nx)+1,ii(l),ii(m))]=p[ic(ii(1),ii(l),ii(m))];
				p[ic(ii(nx)+2,ii(l),ii(m))]=p[ic(ii(2),ii(l),ii(m))];
			}
			//j-layer (k,i varying)
			if (l<=nz && m <= nx) {
				p[ic(ii(m),ii(1)-2,ii(l))]=p[ic(ii(m),ii(ny-1),ii(l))];
				p[ic(ii(m),ii(1)-1,ii(l))]=p[ic(ii(m),ii(ny),ii(l))];
				p[ic(ii(m),ii(ny)+1,ii(l))]=p[ic(ii(m),ii(1),ii(l))];
				p[ic(ii(m),ii(ny)+2,ii(l))]=p[ic(ii(m),ii(2),ii(l))];
			}
			//k-layer (i,j varying)
			if (l<=nx && m <= ny) {
				p[ic(ii(l),ii(m),ii(1)-2)]=p[ic(ii(l),ii(m),ii(2))];
				p[ic(ii(l),ii(m),ii(1)-1)]=p[ic(ii(l),ii(m),ii(1))];
				p[ic(ii(l),ii(m),ii(nz)+1)]=p[ic(ii(l),ii(m),ii(nz))];
				p[ic(ii(l),ii(m),ii(nz)+2)]=p[ic(ii(l),ii(m),ii(nz-1))];
			}
		}
		//edge
		//i-direction (i varying)
		if (l<=nx) {
			p[ic(ii(l),ii(1)-1,ii(1)-1)]=p[ic(ii(l),ii(ny),ii(nz))];
			p[ic(ii(l),ii(1)-1,ii(nz)+1)]=p[ic(ii(l),ii(ny),ii(1))];
			p[ic(ii(l),ii(ny)+1,ii(1)-1)]=p[ic(ii(l),ii(1),ii(nz))];
			p[ic(ii(l),ii(ny)+1,ii(nz)+1)]=p[ic(ii(l),ii(1),ii(1))];
		}
		//j-direction (j varying)
		if (l<=ny) {
			p[ic(ii(1)-1,ii(l),ii(1)-1)]=p[ic(ii(nx),ii(l),ii(nz))];
			p[ic(ii(1)-1,ii(l),ii(nz)+1)]=p[ic(ii(nx),ii(l),ii(1))];
			p[ic(ii(nx)+1,ii(l),ii(1)-1)]=p[ic(ii(1),ii(l),ii(nz))];
			p[ic(ii(nx)+1,ii(l),ii(nz)+1)]=p[ic(ii(1),ii(l),ii(1))];
		}
		//k-direction (k varying)
		if (l<=nz) {
			p[ic(ii(1)-1,ii(1)-1,ii(l))]=p[ic(ii(nx),ii(ny),ii(l))];
			p[ic(ii(1)-1,ii(ny)+1,ii(l))]=p[ic(ii(nx),ii(1),ii(l))];
			p[ic(ii(nx)+1,ii(1)-1,ii(l))]=p[ic(ii(1),ii(ny),ii(l))];
			p[ic(ii(nx)+1,ii(ny)+1,ii(l))]=p[ic(ii(1),ii(1),ii(l))];
		}
	}
	//coners
	p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1)-1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1)-1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1)-1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1)-1)];

	p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz)+1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz)+1)];
}

void fileout_p(double *p[nphase])
{
	///////////
	char fname[100];
	FILE *out;
	int i,j,k;
	///////////
	sprintf(fname,"p%010d.plt",itime);
	out = fopen(fname,"w");
	fprintf(out,"zone, i=%d, j=%d, k=%d\n",nx,ny,nz);
	for (k=1;k<=nz;k++) {
	for (j=1;j<=ny;j++) {
	for (i=1;i<=nx;i++) {
		fprintf(out,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",i,j,k,p[0][ic(ii(i),ii(j),ii(k))],p[1][ic(ii(i),ii(j),ii(k))],p[2][ic(ii(i),ii(j),ii(k))]);
		//fprintf(out,"%d\t%d\t%d\t%lf\t%lf\n",i,j,k,p[0][ii(i)][ii(j)][ii(k)],p[1][ii(i)][ii(j)][ii(k)]);
	}
	}
	}
	fclose(out);
}

double calc_frac_p(double *p)
{
	/////////////
	double sum=0.;
	int i,j,k;
	/////////////
	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
		sum+=p[ic(ii(i),ii(j),ii(k))];
	}
	}
	}

	return sum/(double)(nxcal*nycal*nzcal);
}



void fileout_g(void)
{
    ////////////////
    double Pi,Pj;
    double dp  = 0.01;
    int npoint = (int)(1./dp) + 1;
    FILE *out;
    ////////////////
    out = fopen("doublewell_g.plt","w");
    fprintf(out,"zone,i=%d, j=%d\n",npoint,npoint);
    for (Pj=0.;Pj<1.+dp;Pj+=dp ) {
    for (Pi=0.;Pi<1.+dp;Pi+=dp ) {
        fprintf(out,"%lf\t%lf\t%lf\n",Pi,Pj,g(Pi,Pj));
    }
    }
    fclose(out);
}
void fileout_dg(void)
{
    ////////////////
    double Pi,Pj;
    double dp  = 0.01;
    int npoint = (int)(1./dp) + 1;
    FILE *out;
    ////////////////
    out = fopen("dg.plt","w");
    fprintf(out,"zone,i=%d, j=%d\n",npoint,npoint);
    for (Pj=0.;Pj<1.+dp;Pj+=dp ) {
    for (Pi=0.;Pi<1.+dp;Pi+=dp ) {
        fprintf(out,"%lf\t%lf\t%lf\n",Pi,Pj,dg(Pi,Pj));
    }
    }
    fclose(out);
}

void frac_out(double *p[nphase])
{
    /////////////////
    double calc_frac_p(double *p);
    /////////////////
    FILE *out;
    int j;
    /////////////////
    out = fopen("fraction.txt","a+");
    fprintf(out,"%d\t",itime);
    for (j=0;j<nphase;j++) {
        fprintf(out,"%lf%s",calc_frac_p(p[j]),(j==nphase-1)?"\n":"\t");
    }
    fclose(out);
}

void draw_cube(double *p,int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (abs(i-ix)<=size/2 && abs(j-iy)<=size/2 && abs(k-iz)<=size/2) {
            p[ic(ii(i),ii(j),ii(k))]=1.;
        }
    }
    }
    }

}

void draw_sphere(double *p,int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if ((i-ix)*(i-ix)+(j-iy)*(j-iy)+(k-iz)*(k-iz) <= (size/2)*(size/2)) {
            p[ic(ii(i),ii(j),ii(k))]=1.;
        }
    }
    }
    }

}
void draw_cylinder(double *p,int ix, int iy, int iz, char axis, int rad, int height)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    int nr;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (axis=='z') {
            nr = (i-ix)*(i-ix)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(k-iz)<=height/2) {
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='x'){
            nr = (j-iy)*(j-iy)+(k-iz)*(k-iz);
            if (nr<=rad*rad && abs(i-ix)<=height/2) {
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='y') {
            nr = (k-iz)*(k-iz)+(i-ix)*(i-ix);
            if (nr<=rad*rad && abs(j-iy)<=height/2) {
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else{
            printf("draw_cylinder(): Wrong axis.\n");
        }
    }
    }
    }

}

void draw_rectangular(double *p,int ix, int iy, int iz, int sizex, int sizey,int sizez)
{
    // make a cube at position (i,j,k)
    // dimension: sizex*sizey*sizez
    //////////////////////////////////
    int i,j,k;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (abs(i-ix)<=sizex/2 && abs(j-iy)<=sizey/2 && abs(k-iz)<=sizez/2) {
            p[ic(ii(i),ii(j),ii(k))]=1.;
        }
    }
    }
    }

}

void fill_environ(int phaseindex,double *p[nphase])
{
    ///////////////////////
    int i,j,k,np;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        sum=0.;
        for (np=0;np<nphase;np++) {
            if (np!=phaseindex) {
                 sum += p[np][ic(ii(i),ii(j),ii(k))];
            }
        }
        p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
    }
    }
    }
}
void fill_rectangular(int phaseindex,double *p[nphase],int ix, int iy, int iz, int sizex, int sizey, int sizez)
{
    ///////////////////////
    int i,j,k,np;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (abs(i-ix)<=sizex/2 && abs(j-iy)<=sizey/2 && abs(k-iz)<=sizez/2) {
            sum=0.;
            for (np=0;np<nphase;np++) {
                if (np!=phaseindex) {
                   sum += p[np][ic(ii(i),ii(j),ii(k))];
                }
            }
            p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
        }
    }
    }
    }
}
void fill_cylinder(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int rad, int height)
{
    ///////////////////////
    int i,j,k,np;
    int nr;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (axis=='z') {
            nr = (i-ix)*(i-ix)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(k-iz)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        } else if (axis=='x') {
            nr = (k-iz)*(k-iz)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(i-ix)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        } else if (axis=='y') {
            nr = (i-ix)*(i-ix)+(k-iz)*(k-iz);
            if (nr<=rad*rad && abs(j-iy)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        }else{
            printf("fill_cylinder(): Wrong axis.\n");
        }
    }
    }
    }
}
void clean_all(double *p[nphase])
{
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        for (np=0;np<nphase;np++) {
            p[np][ic(ii(i),ii(j),ii(k))]=0.;
        }
    }
    }
    }
    
}

void clean_one(double *p)
{
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        p[ic(ii(i),ii(j),ii(k))]=0.;
    }
    }
    }
    
}






void calc_lamb(double *pnew[nphase], double *pold[nphase])
{
    //////////////
    double df(int,double *p[nphase],int,int,int);
    //////////////
    extern double lamb[nphase];
    extern double V0[nphase];
    extern double RR[nphase];
    extern double HH[nphase];
    extern double dxl,dyl,dzl;
    extern double EPSMAX, EMM;
    extern int pfix;

    double R[nphase]={0.};
    double H[nphase]={0.};
    int np;
    int i,j,k;
    int iw,jw,kw;
    int plist[nphase];
    //////////////
    dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

    for (np=0;np<nphase;np++) {
        for (i=1;i<=nx;i++) {
        for (j=1;j<=ny;j++) {
        for (k=1;k<=nz;k++) {
            iw=ii(i);
            jw=ii(j);
            kw=ii(k);
            //makephaselist(plist,pold,iw,jw,kw);
            if (itime==0)V0[np]+=pold[np][ic(iw,jw,kw)];
            R[np]+=(pnew[np][ic(iw,jw,kw)])/dt;
            H[np]+=df(np,pnew,iw,jw,kw);
            if(pfix>=0 && pfix<nphase && fabs(pnew[pfix][ic(iw,jw,kw)])>0.) {
                H[np]-=df(np,pnew,iw,jw,kw);
            }
        }
        } 
        }
        R[np]-=V0[np]/dt;
    }
    
    //printf("itime = %d, R[0]=%le, R[1]=%le, R[2]=%le\n",itime,R[0],R[1],R[2]);
    //printf("itime = %d, H[0]=%le, H[1]=%le, H[2]=%le\n",itime,H[0],H[1],H[2]);
    //printf("itime = %d, HH[0]=%le, HH[1]=%le, HH[2]=%le\n",itime,HH[0],HH[1],HH[2]);
    for (np=0;np<nphase;np++) {
        if (fabs(H[np])>= 1.0e-15) {
            //lamb[np] =  R[np]/H[np];
            lamb[np] =   R[np]/H[np];
        }else{
            lamb[np] = 0.;
        }
    }
}

void clean_RR_HH()
{
    ////////////
    extern double RR[nphase];
    extern double HH[nphase];
    int np;
    ////////////
    for (np=0;np<nphase;np++) {
        RR[np]=0.;
        HH[np]=0.;
    }
}
void relax_p(double *p[nphase], int nrelax)
{
    //////////////////
    extern int nx,ny,ny,ntot;
    ///////////////////
    void relaxPF_CH(double *pnew[nphase], double *pold[nphase]);
    void relaxPF_CA(double *pnew[nphase], double *pold[nphase]);
    void copy_p(double *presult[nphase], double *psource[nphase]);
    void update_p(double *pold[nphase], double *pnew[nphase]);
    void BC_periodic3D(double *p);
    void BC_adiabatic3D(double *p);
    ///////////////////
    double *pdum[nphase];
    for (int np=0;np<nphase;np++) {
        pdum[np]=new double[ntot];
    }
    int ir;
    ///////////////////
    copy_p(pdum,p);
    for (ir=0;ir<nrelax;ir++) {
        //update pnew --> p
        update_p(p,pdum);
        //Applying boundary condition on p
        {
                int np;
                for (np=0;np<nphase;np++) {
                        //BC_periodic3D(p[np]);
                        BC_adiabatic3D(p[np]);
                }
        }

        //relaxPF_CH(pdum,p);
        relaxPF_CA(pdum,p);
    }
    copy_p(p,pdum);
}
