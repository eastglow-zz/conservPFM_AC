#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defined_values.h"
#include "common_functions.h"

#include <mpi.h>
int shrink_buffer(int *ista, int*iend,int n1, int n2, int nprocs, int myrank);
void getRANK_NPROCS(int *rank, int *size);
MPI_Request SEND_layer_d( double *p, int index ,int dest, int tag);
MPI_Request RECV_layer_d( double *p, int index ,int from, int tag);
MPI_Request SEND_layer_i( int *p, int index ,int dest, int tag);
MPI_Request RECV_layer_i( int *p, int index ,int from, int tag);
int RANK, NPROCS;
int kmin,kmax;
int RANK_L;
int RANK_U;

int nx = nxcal;
int ny = nycal;
int nz = nzcal;
int ntot = nxtot*nytot*nztot;

int nmax = max3(nx,ny,nz);

extern double time, dt, dxl, dyl, dzl;
extern int itime;

int main(int argc, char *argv[])
{
    /////////////////////
    extern int RANK,NPROCS,RANK_L,RANK_U,kmin,kmax;
    /////////////////////
    void fileout_g(void);
    void fileout_dg(void);
    void initialize_para(double *p[nphase]);
    void BC_periodic3D_para(double *p);
    void BC_adiabatic3D_para(double *p);
    void BC_userdefine3D(double *p);
    void calc_lamb_para(double *pnew[nphase], double *pold[nphase]);
    void fileout_p_para(double *p[nphase]);
    double calc_frac_p(double *p);
    void frac_out_para(double *p[nphase]);
    void clean_RR_HH(void);
    /////////////////////
    double *p[2][nphase]; //phase field old/new switchable
    double *pinter[nphase];
    MPI_Init(&argc,&argv);
    getRANK_NPROCS(&RANK,&NPROCS);
    RANK_L = RANK==0?NPROCS-1:RANK-1;
    RANK_U = RANK==(NPROCS-1)?0:RANK+1;
    nz=shrink_buffer(&kmin,&kmax,1,nzcal,NPROCS,RANK);
    ntot=(nx+4)*(ny+4)*(nz+4);
    /////////////////////
    for (int np=0;np<nphase;np++) {
        p[OLD(itime)][np]=new double[ntot];
        p[NEW(itime)][np]=new double[ntot];
        pinter[np]=new double[ntot];
    }
    /////////////////////
    //fileout_g();
    //fileout_dg();
     
    initialize_para(p[OLD(0)]);
    //temporal iteration
    for (itime=0, time=0.; 1 ;itime++,time+=dt) {
        //Applying boundary condition on p
        {
            int np;
            for (np=0;np<nphase;np++) {
                //BC_periodic3D(p[OLD(itime)][np]);
                BC_adiabatic3D_para(p[OLD(itime)][np]);
                //BC_userdefine3D(p[OLD(itime)][np]);
            }
        }

        //Solving for PF eqn.; FTCS scheme
        //solvPF_CH(p[NEW(itime)],p[NEW(itime)]);
        //solvPF_CA(p[NEW(itime)],p[OLD(itime)]);

        if(1){
            //clean_RR_HH();
            solvPF_CA(pinter,p[OLD(itime)]);
            int np;
            for (np=0;np<nphase;np++) {
                //BC_periodic3D(p[np]);
                BC_adiabatic3D_para(pinter[np]);
               //BC_userdefine3D(p[np]);
            }
            calc_lamb_para(pinter,p[OLD(itime)]);
            drivPF_CA(p[NEW(itime)],pinter);
        } 

		
        if (itime%100  ==0) {
            fileout_p_para(p[OLD(itime)]);
            frac_out_para(p[OLD(itime)]);
            //system("pause");
        }

    }

    for (int np=0;np<nphase;np++) {
        del_buff_d(&p[OLD(itime)][np]);
        del_buff_d(&p[NEW(itime)][np]);
        del_buff_d(&pinter[np]);
    }

    MPI_Finalize();
    //system("pause");
    return 0;
}

void initialize_para(double *p[nphase])
{
        ///////////
        void clean_all(double *p[nphase]);
        void clean_one(double *p);
        void relax_p_para(double *p[nphase], int nrelax);
        void copy_p(double *presult[nphase], double *psource[nphase]);
        void draw_cube_para(double *p,int ix, int iy, int iz, int size);
        void draw_sphere_para(double *p,int ix, int iy, int iz, int diameter);
        void draw_cylinder_para(double *p,int ix, int iy, int iz, char axis, int rad, int height);
        void draw_equilat_tri_prism_para(double *p,int ix, int iy, int iz, char axis, int nlat, int height);
        void draw_rectangular_para(double *p,int ix, int iy, int iz, int sizex, int sizey,int sizez);
        void fill_rectangular_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, int sizex, int sizey, int sizez);
        void fill_cylinder_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int rad, int height);
        void fill_equilat_tri_prism_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int nlat, int height);
        void fill_environ(int phaseindex,double *p[nphase]);
	///////////
	clean_all(p); 
        draw_equilat_tri_prism_para(p[2],nx/2,ny/2,15,'z',42*2,15*2);
        fill_equilat_tri_prism_para(1,p,nx/2,ny/2,15*2+15,'z',42*2,15*2);
        //draw_rectangular_para(p[1],5,nycal/2,nzcal/2 ,10,nycal,nzcal);
        fill_environ(0,p);
        relax_p_para(p, 50);
        //clean_one(p[0]);
        //fill_cylinder_para(2,p,10+10,nycal/2,nzcal/2,'x',20,20);
        //fill_environ(0,p);
}


void BC_adiabatic3D_para(double *p)
{
        //////////////////
        extern int RANK,NPROCS,RANK_L,RANK_U;
	//////////////////
	int l,m,n;
	int lim;
        int lower,upper;
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
                    if (RANK==0) {
                        p[ic(ii(l),ii(m),ii(1)-2)]=p[ic(ii(l),ii(m),ii(2))];
                        p[ic(ii(l),ii(m),ii(1)-1)]=p[ic(ii(l),ii(m),ii(1))];
                    }else if (RANK==NPROCS-1) {
                        p[ic(ii(l),ii(m),ii(nz)+1)]=p[ic(ii(l),ii(m),ii(nz))];
                        p[ic(ii(l),ii(m),ii(nz)+2)]=p[ic(ii(l),ii(m),ii(nz-1))];
                    }
                }
            }
	    //edge
	    //i-direction (i varying)
	    if (l<=nx) {
                if (RANK==0) {
                    p[ic(ii(l),ii(1)-1,ii(1)-1)]=p[ic(ii(l),ii(1),ii(1))];
                    p[ic(ii(l),ii(ny)+1,ii(1)-1)]=p[ic(ii(l),ii(ny),ii(1))];
                }else if (RANK==NPROCS-1) {
                    p[ic(ii(l),ii(ny)+1,ii(nz)+1)]=p[ic(ii(l),ii(ny),ii(nz))];
                    p[ic(ii(l),ii(1)-1,ii(nz)+1)]=p[ic(ii(l),ii(1),ii(nz))];
                }
	    }
            
	    //j-direction (j varying)
	    if (l<=ny) {
                if (RANK==0) {
                    p[ic(ii(1)-1,ii(l),ii(1)-1)]=p[ic(ii(1),ii(l),ii(1))];
                    p[ic(ii(1)-1,ii(l),ii(nz)+1)]=p[ic(ii(1),ii(l),ii(nz))];
                }else if(RANK==NPROCS-1) {
                    p[ic(ii(nx)+1,ii(l),ii(1)-1)]=p[ic(ii(nx),ii(l),ii(1))];
                    p[ic(ii(nx)+1,ii(l),ii(nz)+1)]=p[ic(ii(nx),ii(l),ii(nz))];
                }
	    }
	    //k-direction (k varying)
	    if (l<=nz) {
                p[ic(ii(1)-1,ii(1)-1,ii(l))]=p[ic(ii(1),ii(1),ii(l))];
                p[ic(ii(1)-1,ii(ny)+1,ii(l))]=p[ic(ii(1),ii(ny),ii(l))];
                p[ic(ii(nx)+1,ii(1)-1,ii(l))]=p[ic(ii(nx),ii(1),ii(l))];
                p[ic(ii(nx)+1,ii(ny)+1,ii(l))]=p[ic(ii(nx),ii(ny),ii(l))];
	    }
            
        }
        if (RANK==0) {
            p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1))];
            p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1))];
            p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1))];
            p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1))];
        }
        if (RANK==NPROCS-1) {
            p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz))];
            p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz))];
            p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz))];
            p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz))];
        }
        if (RANK>0 && RANK<NPROCS-1) {
            MPI_Request req1,req2,req3,req4;
            MPI_Status stat1,stat2,stat3,stat4;
            MPI_Request req5,req6,req7,req8;
            MPI_Status stat5,stat6,stat7,stat8;
            req1 = SEND_layer_d(p, ii(1), RANK_L, 1);
            req2 = SEND_layer_d(p, ii(2), RANK_L, 2);
            req3 = RECV_layer_d(p, ii(nz)+1, RANK_U, 1);
            req4 = RECV_layer_d(p, ii(nz)+2, RANK_U, 2);

            req5 = SEND_layer_d(p, ii(nz), RANK_U, 3);
            req6 = SEND_layer_d(p, ii(nz-1), RANK_U, 4);
            req7 = RECV_layer_d(p, ii(1)-1 , RANK_L, 3);
            req8 = RECV_layer_d(p, ii(1)-2 , RANK_L, 4);

            MPI_Wait(&req1,&stat1);
            MPI_Wait(&req2,&stat2);
            MPI_Wait(&req3,&stat3);
            MPI_Wait(&req4,&stat4);
            MPI_Wait(&req5,&stat5);
            MPI_Wait(&req6,&stat6);
            MPI_Wait(&req7,&stat7);
            MPI_Wait(&req8,&stat8);
        }else if (RANK==0) {
            MPI_Request req5,req6,req7,req8;
            MPI_Status stat5,stat6,stat7,stat8;

            req5 = SEND_layer_d(p, ii(nz), RANK_U, 3);
            req6 = SEND_layer_d(p, ii(nz-1), RANK_U, 4);
            req7 = RECV_layer_d(p, ii(nz)+1 , RANK_U, 1);
            req8 = RECV_layer_d(p, ii(nz)+2 , RANK_U, 2);
            MPI_Wait(&req5,&stat5);
            MPI_Wait(&req6,&stat6);
            MPI_Wait(&req7,&stat7);
            MPI_Wait(&req8,&stat8);
        }else if (RANK==NPROCS-1) {
            MPI_Request req1,req2,req3,req4;
            MPI_Status stat1,stat2,stat3,stat4;

            req1 = SEND_layer_d(p, ii(1), RANK_L, 1);
            req2 = SEND_layer_d(p, ii(2), RANK_L, 2);
            req3 = RECV_layer_d(p, ii(1)-1, RANK_L, 3);
            req4 = RECV_layer_d(p, ii(1)-2, RANK_L, 4);
            MPI_Wait(&req1,&stat1);
            MPI_Wait(&req2,&stat2);
            MPI_Wait(&req3,&stat3);
            MPI_Wait(&req4,&stat4);
        }
}

void BC_periodic3D_para(double *p)
{
        //////////////////
        extern int RANK,NPROCS,RANK_L,RANK_U;
	//////////////////
	int l,m,n;
	int lim;
	//////////////////
	
	for (l=1;l<=nmax;l++) {
		//layer
		for (m=1;m<=nmax;m++) {
			//i-layer (j,k varying)
			/*
			if (l<=ny && m <= nz) {
				p[ic(ii(1)-2,ii(l),ii(m))]=p[ic(ii(nx-1),ii(l),ii(m))];
				p[ic(ii(1)-1,ii(l),ii(m))]=p[ic(ii(nx),ii(l),ii(m))];
				p[ic(ii(nx)+1,ii(l),ii(m))]=p[ic(ii(1),ii(l),ii(m))];
				p[ic(ii(nx)+2,ii(l),ii(m))]=p[ic(ii(2),ii(l),ii(m))];
			}*/

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
		/*
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
                */
	}
        /*
	//coners
	p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1)-1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1)-1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1)-1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1)-1)];

	p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz)+1)];
	p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz)+1)];
	p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz)+1)];
        */
        if (1) {
            MPI_Request req1,req2,req3,req4;
            MPI_Status stat1,stat2,stat3,stat4;
            MPI_Request req5,req6,req7,req8;
            MPI_Status stat5,stat6,stat7,stat8;

            req1 = SEND_layer_d(p, ii(1), RANK_L, 1);
            req2 = SEND_layer_d(p, ii(2), RANK_L, 2);
            req3 = RECV_layer_d(p, ii(nx)+1, RANK_U, 1);
            req4 = RECV_layer_d(p, ii(nx)+2, RANK_U, 2);

            req5 = SEND_layer_d(p, ii(nx), RANK_U, 3);
            req6 = SEND_layer_d(p, ii(nx-1), RANK_U, 4);
            req7 = RECV_layer_d(p, ii(1)-1 , RANK_L, 3);
            req8 = RECV_layer_d(p, ii(1)-2 , RANK_L, 4);

            MPI_Wait(&req1,&stat1);
            MPI_Wait(&req2,&stat2);
            MPI_Wait(&req3,&stat3);
            MPI_Wait(&req4,&stat4);
            MPI_Wait(&req5,&stat5);
            MPI_Wait(&req6,&stat6);
            MPI_Wait(&req7,&stat7);
            MPI_Wait(&req8,&stat8);
        }
}

void BC_userdefine3D_para(double *p)
{
        //x: adiabatic
        //y,z: periodic
        //////////////////
        extern int RANK,NPROCS,RANK_L,RANK_U;
	//////////////////
	int l,m,n;
	int lim;
	//////////////////
        for (l=1;l<=nmax;l++) {
                //layer
                for (m=1;m<=nmax;m++) {
                        //i-layer (j,k varying)
                        /*
                        if (l<=ny && m <= nz) {
                                p[ic(ii(1)-2,ii(l),ii(m))]=p[ic(ii(nx-1),ii(l),ii(m))];
                                p[ic(ii(1)-1,ii(l),ii(m))]=p[ic(ii(nx),ii(l),ii(m))];
                                p[ic(ii(nx)+1,ii(l),ii(m))]=p[ic(ii(1),ii(l),ii(m))];
                                p[ic(ii(nx)+2,ii(l),ii(m))]=p[ic(ii(2),ii(l),ii(m))];
                        }*/

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
                /*
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
                */
        }
        /*
        //coners
        p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1)-1)];
        p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1)-1)];
        p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1)-1)];
        p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1)-1)];

        p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz)+1)];
        p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz)+1)];
        p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz)+1)];
        p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz)+1)];
        */
        if (RANK==0) {
            p[ic(ii(1)-1,ii(1)-1,ii(1)-1)]=p[ic(ii(1),ii(1),ii(1))];
            p[ic(ii(1)-1,ii(ny)+1,ii(1)-1)]=p[ic(ii(1),ii(ny),ii(1))];
            p[ic(ii(1)-1,ii(1)-1,ii(nz)+1)]=p[ic(ii(1),ii(1),ii(nz))];
            p[ic(ii(1)-1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(1),ii(ny),ii(nz))];
        }
        if (RANK==NPROCS-1) {
            p[ic(ii(nx)+1,ii(1)-1,ii(1)-1)]=p[ic(ii(nx),ii(1),ii(1))];
            p[ic(ii(nx)+1,ii(ny)+1,ii(1)-1)]=p[ic(ii(nx),ii(ny),ii(1))];
            p[ic(ii(nx)+1,ii(1)-1,ii(nz)+1)]=p[ic(ii(nx),ii(1),ii(nz))];
            p[ic(ii(nx)+1,ii(ny)+1,ii(nz)+1)]=p[ic(ii(nx),ii(ny),ii(nz))];
        }
        
        if (RANK>0 && RANK<NPROCS-1) {
            MPI_Request req1,req2,req3,req4;
            MPI_Status stat1,stat2,stat3,stat4;
            MPI_Request req5,req6,req7,req8;
            MPI_Status stat5,stat6,stat7,stat8;

            req1 = SEND_layer_d(p, ii(1), RANK_L, 1);
            req2 = SEND_layer_d(p, ii(2), RANK_L, 2);
            req3 = RECV_layer_d(p, ii(nx)+1, RANK_U, 1);
            req4 = RECV_layer_d(p, ii(nx)+2, RANK_U, 2);

            req5 = SEND_layer_d(p, ii(nx), RANK_U, 3);
            req6 = SEND_layer_d(p, ii(nx-1), RANK_U, 4);
            req7 = RECV_layer_d(p, ii(1)-1 , RANK_L, 3);
            req8 = RECV_layer_d(p, ii(1)-2 , RANK_L, 4);

            MPI_Wait(&req1,&stat1);
            MPI_Wait(&req2,&stat2);
            MPI_Wait(&req3,&stat3);
            MPI_Wait(&req4,&stat4);
            MPI_Wait(&req5,&stat5);
            MPI_Wait(&req6,&stat6);
            MPI_Wait(&req7,&stat7);
            MPI_Wait(&req8,&stat8);
        }else if (RANK==0) {
            MPI_Request req5,req6,req7,req8;
            MPI_Status stat5,stat6,stat7,stat8;

            req5 = SEND_layer_d(p, ii(nx), RANK_U, 3);
            req6 = SEND_layer_d(p, ii(nx-1), RANK_U, 4);
            req7 = RECV_layer_d(p, ii(1)-1 , RANK_L, 3);
            req8 = RECV_layer_d(p, ii(1)-2 , RANK_L, 4);
            MPI_Wait(&req5,&stat5);
            MPI_Wait(&req6,&stat6);
            MPI_Wait(&req7,&stat7);
            MPI_Wait(&req8,&stat8);
        }else if (RANK==NPROCS-1) {
            MPI_Request req1,req2,req3,req4;
            MPI_Status stat1,stat2,stat3,stat4;

            req1 = SEND_layer_d(p, ii(1), RANK_L, 1);
            req2 = SEND_layer_d(p, ii(2), RANK_L, 2);
            req3 = RECV_layer_d(p, ii(nx)+1, RANK_U, 1);
            req4 = RECV_layer_d(p, ii(nx)+2, RANK_U, 2);
            MPI_Wait(&req1,&stat1);
            MPI_Wait(&req2,&stat2);
            MPI_Wait(&req3,&stat3);
            MPI_Wait(&req4,&stat4);
        }
}

void fileout_p_para(double *p[nphase])
{
    //////////////////
    extern int RANK,NPROCS,RANK_L,RANK_U;
    ///////////
    char fname[100];
    FILE *out;
    int i,j,k;
    int kcal;
    int irank;
    int gosignal=1;
    MPI_Status status1;
    ///////////
    sprintf(fname,"p%010d.plt",itime);

    if (RANK!=0) MPI_Recv(&gosignal,1,MPI_INT,RANK_L,gosignal++,MPI_COMM_WORLD,&status1);
    out = fopen(fname,"a+");
    if (RANK==0) fprintf(out,"zone, i=%d, j=%d, k=%d\n",nxcal,nycal,nzcal);
    for (k=1;k<=nz;k++) {
    for (j=1;j<=ny;j++) {
    for (i=1;i<=nx;i++) {
        kcal=RANK*nz+k;
        fprintf(out,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",i,j,kcal,p[0][ic(ii(i),ii(j),ii(k))],p[1][ic(ii(i),ii(j),ii(k))],p[2][ic(ii(i),ii(j),ii(k))]);
        //fprintf(out,"%d\t%d\t%d\t%lf\t%lf\n",i,j,k,p[0][ii(i)][ii(j)][ii(k)],p[1][ii(i)][ii(j)][ii(k)]);
    }
    }
    }
    fclose(out);
    if (RANK!=NPROCS-1) MPI_Ssend(&gosignal,1,MPI_INT,RANK_U,gosignal,MPI_COMM_WORLD);
}

double calc_frac_p(double *p )
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

void frac_out_para(double *p[nphase])
{
    /////////////////
    double calc_frac_p(double *p);
    /////////////////
    FILE *out;
    double frac_1procs[nphase];
    double frac[nphase];
    int j;
    /////////////////
    for (j=0;j<nphase;j++) {
        frac_1procs[j]=calc_frac_p(p[j]);
        MPI_Allreduce(&frac_1procs[j],&frac[j],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    if (RANK==0) {
        out = fopen("fraction.txt","a+");
        fprintf(out,"%d\t",itime);
        for (j=0;j<nphase;j++) {
            fprintf(out,"%lf%s",frac[j],(j==nphase-1)?"\n":"\t");
        }
        fclose(out);
    }
}

void draw_cube_para(double *p,int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    int kcal;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (abs(i-ix)<=size/2 && abs(j-iy)<=size/2 && abs(kcal-iz)<=size/2) {
            p[ic(ii(i),ii(j),ii(k))]=1.;
        }
    }
    }
    }

}

void draw_sphere_para(double *p,int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    int kcal;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if ((i-ix)*(i-ix)+(j-iy)*(j-iy)+(kcal-iz)*(kcal-iz) <= (size/2)*(size/2)) {
            p[ic(ii(i),ii(j),ii(k))]=1.;
        }
    }
    }
    }

}
void draw_cylinder_para(double *p,int ix, int iy, int iz, char axis, int rad, int height)
{
    //////////////////////////////////
    int i,j,k;
    int kcal;
    int nr;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (axis=='z') {
            nr = (i-ix)*(i-ix)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(kcal-iz)<=height/2) {
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='x'){
            nr = (j-iy)*(j-iy)+(kcal-iz)*(kcal-iz);
            if (nr<=rad*rad && abs(i-ix)<=height/2) {
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='y') {
            nr = (kcal-iz)*(kcal-iz)+(i-ix)*(i-ix);
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

void draw_equilat_tri_prism_para(double *p,int ix, int iy, int iz, char axis, int nlat, int height)
{
    // make a equilateral triangular prism at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    int kcal;
    int nseg1,nseg2,nseg3;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (axis=='z') {
            nseg1 = iy-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(i-ix)+iy+nlat/2);
            nseg3 = (int)( sqrt(3.)*(i-ix)+iy+nlat/2);
            if (j >= nseg1 && j <= nseg2 && j <= nseg3 && abs(kcal-iz)<=height/2){
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='x'){
            nseg1 = iz-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(j-iy)+iz+nlat/2);
            nseg3 = (int)( sqrt(3.)*(j-iy)+iz+nlat/2);
            if (kcal >= nseg1 && kcal <= nseg2 && kcal <= nseg3 && abs(i-ix)<=height/2){
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else if(axis=='y') {
            nseg1 = ix-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(kcal-iz)+ix+nlat/2);
            nseg3 = (int)( sqrt(3.)*(kcal-iz)+ix+nlat/2);
            if (j >= nseg1 && j <= nseg2 && j <= nseg3 && abs(j-iy)<=height/2){
                p[ic(ii(i),ii(j),ii(k))]=1.;
            }
        }else{
            printf("draw_cylinder(): Wrong axis.\n");
        }
    }
    }
    }

}

void draw_rectangular_para(double *p,int ix, int iy, int iz, int sizex, int sizey,int sizez)
{
    // make a cube at position (i,j,k)
    // dimension: sizex*sizey*sizez
    //////////////////////////////////
    int i,j,k;
    int kcal;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (abs(i-ix)<=sizex/2 && abs(j-iy)<=sizey/2 && abs(kcal-iz)<=sizez/2) {
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
void fill_rectangular_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, int sizex, int sizey, int sizez)
{
    ///////////////////////
    int i,j,k,np;
    int kcal;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (abs(i-ix)<=sizex/2 && abs(j-iy)<=sizey/2 && abs(kcal-iz)<=sizez/2) {
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
void fill_cylinder_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int rad, int height)
{
    ///////////////////////
    int i,j,k,np;
    int kcal;
    int nr;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (axis=='z') {
            nr = (i-ix)*(i-ix)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(kcal-iz)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        } else if (axis=='x') {
            nr = (kcal-iz)*(kcal-iz)+(j-iy)*(j-iy);
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
            nr = (i-ix)*(i-ix)+(kcal-iz)*(kcal-iz);
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
void fill_equilat_tri_prism_para(int phaseindex,double *p[nphase],int ix, int iy, int iz, char axis, int nlat, int height)
{
    ///////////////////////
    int i,j,k,np;
    int kcal;
    int nseg1, nseg2, nseg3;
    double sum;
    ///////////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        kcal=RANK*nz+k;
        if (axis=='z') {
            nseg1 = iy-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(i-ix)+iy+nlat/2);
            nseg3 = (int)( sqrt(3.)*(i-ix)+iy+nlat/2);
            if (j >= nseg1 && j <= nseg2 && j <= nseg3 && abs(kcal-iz)<=height/2){
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        }else if(axis=='x'){
            nseg1 = iz-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(j-iy)+iz+nlat/2);
            nseg3 = (int)( sqrt(3.)*(j-iy)+iz+nlat/2);
            if (kcal >= nseg1 && kcal <= nseg2 && kcal <= nseg3 && abs(i-ix)<=height/2){
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        }else if(axis=='y') {
            nseg1 = ix-nlat/4;
            nseg2 = (int)(-sqrt(3.)*(kcal-iz)+ix+nlat/2);
            nseg3 = (int)( sqrt(3.)*(kcal-iz)+ix+nlat/2);
            if (j >= nseg1 && j <= nseg2 && j <= nseg3 && abs(j-iy)<=height/2){
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ic(ii(i),ii(j),ii(k))];
                    }
                }
                p[phaseindex][ic(ii(i),ii(j),ii(k))] = 1.-sum;
            }
        }else{
            printf("draw_cylinder(): Wrong axis.\n");
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






void calc_lamb_para(double *pnew[nphase], double *pold[nphase])
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

    double R_1procs[nphase]={0.};
    double H_1procs[nphase]={0.};
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
            R_1procs[np]+=(pnew[np][ic(iw,jw,kw)])/dt;
            H_1procs[np]+=df(np,pnew,iw,jw,kw);
            if(pfix>=0 && pfix<nphase && fabs(pnew[pfix][ic(iw,jw,kw)])>0.) {
                H_1procs[np]-=df(np,pnew,iw,jw,kw);
            }
        }
        } 
        }
        R_1procs[np]-=V0[np]/dt;
    }
    
    for (np=0;np<nphase;np++) {
        MPI_Allreduce(&R_1procs[np],&R[np],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&H_1procs[np],&H[np],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
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


MPI_Request SEND_layer_d(double *p, int index ,int dest, int tag)
{
    //////////////////
    extern int nx,ny,nz;
    MPI_Request request;
    //////////////////
    //MPI_Isend(&p[ic(index,ii(1)-2,ii(1)-2)],(ny+4)*(nz+4),MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &request);
    MPI_Isend(&p[ic(ii(1)-2,ii(1)-2,index)],(nx+4)*(ny+4),MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &request);
    return request;
}

MPI_Request SEND_layer_i(int *p, int index ,int dest, int tag)
{
    //////////////////
    extern int nx,ny,nz;
    MPI_Request request;
    //////////////////
    //MPI_Isend(&p[ic(index,ii(1)-2,ii(1)-2)],(ny+4)*(nz+4),MPI_INT, dest, tag, MPI_COMM_WORLD, &request);
    MPI_Isend(&p[ic(ii(1)-2,ii(1)-2,index)],(nx+4)*(ny+4),MPI_INT, dest, tag, MPI_COMM_WORLD, &request);
    return request;
}

MPI_Request RECV_layer_d(double *p, int index ,int from, int tag)
{
    //////////////////
    extern int nx,ny,nz;
    MPI_Request request;
    //////////////////
    //MPI_Irecv(&p[ic(index,ii(1)-2,ii(1)-2)],(ny+4)*(nz+4),MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &request);
    MPI_Irecv(&p[ic(ii(1)-2,ii(1)-2,index)],(nx+4)*(ny+4),MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &request);
    return request;
}

MPI_Request RECV_layer_i(int *p, int index ,int from, int tag)
{
    //////////////////
    extern int nx,ny,nz;
    MPI_Request request;
    //////////////////
    //MPI_Irecv(&p[ic(index,ii(1)-2,ii(1)-2)],(ny+4)*(nz+4),MPI_INT, from, tag, MPI_COMM_WORLD, &request);
    MPI_Irecv(&p[ic(ii(1)-2,ii(1)-2,index)],(nx+4)*(ny+4),MPI_INT, from, tag, MPI_COMM_WORLD, &request);
    return request;
}

void relax_p_para(double *p[nphase], int nrelax)
{
    //////////////////
    extern int nx,ny,ny,ntot;
    ///////////////////
    void relaxPF_CH(double *pnew[nphase], double *pold[nphase]);
    void relaxPF_CA(double *pnew[nphase], double *pold[nphase]);
    void copy_p(double *presult[nphase], double *psource[nphase]);
    void update_p(double *pold[nphase], double *pnew[nphase]);
    void BC_periodic3D_para(double *p);
    void BC_adiabatic3D_para(double *p);
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
                        BC_adiabatic3D_para(p[np]);
                }
        }
        
        //relaxPF_CH(pdum,p);
        relaxPF_CA(pdum,p);
    }
    copy_p(p,pdum);
    for (int np=0;np<nphase;np++) {
        del_buff_d(&pdum[np]);
    }
}
