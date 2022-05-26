#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double h(double phi);
double g(double phi1,double phi2);
double dg(double phi1,double phi2);
int max3(int a1, int a2, int a3);
double dmax3(double a1, double a2, double a3); 
double SIG(int np1, int np2);
double W(int np1, int np2);
double EPS(int np1, int np2);
double EMMp(int np1, int np2);
double EMMrlx(int np1, int np2);

#define ii(i) i+1
#define OLD(x) x%2
#define NEW(x) !(x%2)

#define nphase 3

#define nxtot 184
#define nytot 184
#define nztot 114
#define nx 180
#define ny 180
#define nz 110

#define tstable 0.9

#define pi 3.141592

int nmax = max3(nx,ny,nz);

int pfix =  1;

double dxl = 0.5e-10;
double dyl = dxl;
double dzl = dxl;
double XI = 3.*dxl; // half interface thickness

//double SIG = 1.;
double SIGSV = 1.;
double SIGLV = 1.;
double SIGSL = 0.2014;

/*
//polynomial type double well
//double W = 3.*2.94*SIG/XI;
double WSV = 3.*2.94*SIGSV/XI;
double WLV = 3.*2.94*SIGLV/XI;
double WSL = 3.*2.94*SIGSL/XI;
//double EPS = sqrt(6.*SIG*XI/2.94);
double EPSSV = sqrt(6.*SIGSV*XI/2.94);
double EPSLV = sqrt(6.*SIGLV*XI/2.94);
double EPSSL = sqrt(6.*SIGSL*XI/2.94);
*/

//parabolic type double well
double WSV=SIGSV/2./XI;
double WLV=SIGLV/2./XI;
double WSL=SIGSL/2./XI;
double EPSSV=(2./pi)*sqrt(SIGSV*XI);
double EPSLV=(2./pi)*sqrt(SIGLV*XI);
double EPSSL=(2./pi)*sqrt(SIGSL*XI);

double EPSMAX = dmax3(EPSSV,EPSLV,EPSSL);

double EMM=1.;

double dt = pow(dxl,4.)/(12.*EPSMAX*EPSMAX*EMM)*0.1/3.;
double time = 0.;
int itime = 0;

double RR[nphase]={0.};
double V0[nphase]={0.};
double HH[nphase]={0.};
double lamb[nphase]; //interphase driving force coefficient

int main()
{
	/////////////////////
	void fileout_g(void);
	void fileout_dg(void);
        void calc_s(double s[nphase][nxtot][nytot][nztot],double p[nphase][nxtot][nytot][nztot]);
	void initialize(double p[nphase][nxtot][nytot][nztot]);
	void update_p(double pold[nphase][nxtot][nytot][nztot], double pnew[nphase][nxtot][nytot][nztot]);
	void BC_periodic3D(double p[nxtot][nytot][nztot]);
	void BC_adiabatic3D(double p[nxtot][nytot][nztot]);
	void BC_userdefine3D(double p[nxtot][nytot][nztot]);
	void solvPF_CH(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
	void solvPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
	void drivPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
        void calc_lamb(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
	void fileout_p(double p[nphase][nxtot][nytot][nztot]);
	double calc_frac_p(double p[nxtot][nytot][nztot]);
        void frac_out(double p[nphase][nxtot][nytot][nztot]);
        void clean_RR_HH(void);
	/////////////////////
        double p[2][nphase][nxtot][nytot][nztot]; //phase field old/new switchable
	//double pnew[nphase][nxtot][nytot][nztot];
	double pinter[nphase][nxtot][nytot][nztot];
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
		//solvPF_CH(p[NEW(itime)],p[NEW(itime)]);

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
		//solvPF_CA(pnew,pinter);
                //printf("itime:%d, RR[0]=%le, RR[1]=%le\n",itime,RR[0],RR[1]);

                //printf("itime:%d, lamb[0]=%le, lamb[1]=%le, lamb[2]=%le\n",itime,lamb[0],lamb[1],lamb[2]);

		
		if (itime%500  ==0) {
			fileout_p(p[OLD(itime)]);
                        frac_out(p[OLD(itime)]);
			//system("pause");
		}

	}


	//system("pause");
	return 0;
}

void initialize(double p[nphase][nxtot][nytot][nztot])
{
        ///////////
        void clean_all(double p[nphase][nxtot][nytot][nztot]);
        void clean_one(double p[nxtot][nytot][nztot]);
        void relax_p(double p[nphase][nxtot][nytot][nztot], int nrelax);
        void copy_p(double presult[nphase][nxtot][nytot][nztot], double psource[nphase][nxtot][nytot][nztot]);
        void draw_cube(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int size);
        void draw_sphere(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int diameter);
        void draw_cylinder(double p[nxtot][nytot][nztot],int ix, int iy, int iz, char axis, int rad, int height);
        void draw_rectangular(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int sizex, int sizey,int sizez);
        void fill_rectangular(int phaseindex,double p[nphase][nxtot][nytot][nztot],int ix, int iy, int iz, int sizex, int sizey, int sizez);
        void fill_cylinder(int phaseindex,double p[nphase][nxtot][nytot][nztot],int ix, int iy, int iz, char axis, int rad, int height);
        void fill_environ(int phaseindex,double p[nphase][nxtot][nytot][nztot]);
	///////////
	clean_all(p); 
        //draw_cylinder(p[2],nx/2,ny/2,nz/2+6,'z',16,10);
        //fill_cylinder(1,p,nx/2,ny/2,nz/2-6,'z',16,10);
        draw_rectangular(p[1],nx/2,ny/2,10,nx,ny,20);
        fill_environ(0,p);
        relax_p(p, 50);
        clean_one(p[0]);
        fill_cylinder(2,p,nx/2,ny/2,20+40,'z',48,40);
        fill_environ(0,p);
}

void update_p(double pold[nphase][nxtot][nytot][nztot], double pnew[nphase][nxtot][nytot][nztot])
{
	////////////////
	int i,j,k,np;
	////////////////
	for (np=0;np<nphase;np++) {
		for (i=1;i<=nx;i++) {
		for (j=1;j<=ny;j++) {
		for (k=1;k<=nz;k++) {
			pold[np][ii(i)][ii(j)][ii(k)]=pnew[np][ii(i)][ii(j)][ii(k)];
		}
		}
		}
	}
}
void BC_adiabatic3D(double p[nxtot][nytot][nztot])
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
                    p[ii(1)-2][ii(l)][ii(m)]=p[ii(2)][ii(l)][ii(m)];
                    p[ii(1)-1][ii(l)][ii(m)]=p[ii(1)][ii(l)][ii(m)];
                    p[ii(nx)+1][ii(l)][ii(m)]=p[ii(nx)][ii(l)][ii(m)];
                    p[ii(nx)+2][ii(l)][ii(m)]=p[ii(nx-1)][ii(l)][ii(m)];
                }
                //j-layer (k,i varying)
                if (l<=nz && m <= nx) {
                    p[ii(m)][ii(1)-2][ii(l)]=p[ii(m)][ii(2)][ii(l)];
                    p[ii(m)][ii(1)-1][ii(l)]=p[ii(m)][ii(1)][ii(l)];
                    p[ii(m)][ii(ny)+1][ii(l)]=p[ii(m)][ii(ny)][ii(l)];
                    p[ii(m)][ii(ny)+2][ii(l)]=p[ii(m)][ii(ny-1)][ii(l)];
                }
                //k-layer (i,j varying)
                if (l<=nx && m <= ny) {
                    p[ii(l)][ii(m)][ii(1)-2]=p[ii(l)][ii(m)][ii(2)];
                    p[ii(l)][ii(m)][ii(1)-1]=p[ii(l)][ii(m)][ii(1)];
                    p[ii(l)][ii(m)][ii(nz)+1]=p[ii(l)][ii(m)][ii(nz)];
                    p[ii(l)][ii(m)][ii(nz)+2]=p[ii(l)][ii(m)][ii(nz-1)];
                }
            }
	    //edge
	    //i-direction (i varying)
	    if (l<=nx) {
	    	p[ii(l)][ii(1)-1][ii(1)-1]=p[ii(l)][ii(1)][ii(1)];
		p[ii(l)][ii(1)-1][ii(nz)+1]=p[ii(l)][ii(1)][ii(nz)];
		p[ii(l)][ii(ny)+1][ii(1)-1]=p[ii(l)][ii(ny)][ii(1)];
		p[ii(l)][ii(ny)+1][ii(nz)+1]=p[ii(l)][ii(ny)][ii(nz)];
	    }
	    //j-direction (j varying)
	    if (l<=ny) {
		p[ii(1)-1][ii(l)][ii(1)-1]=p[ii(1)][ii(l)][ii(1)];
		p[ii(1)-1][ii(l)][ii(nz)+1]=p[ii(1)][ii(l)][ii(nz)];
		p[ii(nx)+1][ii(l)][ii(1)-1]=p[ii(nx)][ii(l)][ii(1)];
		p[ii(nx)+1][ii(l)][ii(nz)+1]=p[ii(nx)][ii(l)][ii(nz)];
	    }
	    //k-direction (k varying)
	    if (l<=nz) {
	 	p[ii(1)-1][ii(1)-1][ii(l)]=p[ii(1)][ii(1)][ii(l)];
		p[ii(1)-1][ii(ny)+1][ii(l)]=p[ii(1)][ii(ny)][ii(l)];
		p[ii(nx)+1][ii(1)-1][ii(l)]=p[ii(nx)][ii(1)][ii(l)];
		p[ii(nx)+1][ii(ny)+1][ii(l)]=p[ii(nx)][ii(ny)][ii(l)];
	    }
        }
	p[ii(1)-1][ii(1)-1][ii(1)-1]=p[ii(1)][ii(1)][ii(1)];
	p[ii(1)-1][ii(ny)+1][ii(1)-1]=p[ii(1)][ii(ny)][ii(1)];
	p[ii(nx)+1][ii(1)-1][ii(1)-1]=p[ii(nx)][ii(1)][ii(1)];
	p[ii(nx)+1][ii(ny)+1][ii(1)-1]=p[ii(nx)][ii(ny)][ii(1)];

	p[ii(1)-1][ii(1)-1][ii(nz)+1]=p[ii(1)][ii(1)][ii(nz)];
	p[ii(1)-1][ii(ny)+1][ii(nz)+1]=p[ii(1)][ii(ny)][ii(nz)];
	p[ii(nx)+1][ii(1)-1][ii(nz)+1]=p[ii(nx)][ii(1)][ii(nz)];
	p[ii(nx)+1][ii(ny)+1][ii(nz)+1]=p[ii(nx)][ii(ny)][ii(nz)];
}

void BC_periodic3D(double p[nxtot][nytot][nztot])
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
				p[ii(1)-2][ii(l)][ii(m)]=p[ii(nx-1)][ii(l)][ii(m)];
				p[ii(1)-1][ii(l)][ii(m)]=p[ii(nx)][ii(l)][ii(m)];
				p[ii(nx)+1][ii(l)][ii(m)]=p[ii(1)][ii(l)][ii(m)];
				p[ii(nx)+2][ii(l)][ii(m)]=p[ii(2)][ii(l)][ii(m)];
			}
			//j-layer (k,i varying)
			if (l<=nz && m <= nx) {
				p[ii(m)][ii(1)-2][ii(l)]=p[ii(m)][ii(ny-1)][ii(l)];
				p[ii(m)][ii(1)-1][ii(l)]=p[ii(m)][ii(ny)][ii(l)];
				p[ii(m)][ii(ny)+1][ii(l)]=p[ii(m)][ii(1)][ii(l)];
				p[ii(m)][ii(ny)+2][ii(l)]=p[ii(m)][ii(2)][ii(l)];
			}
			//k-layer (i,j varying)
			if (l<=nx && m <= ny) {
				p[ii(l)][ii(m)][ii(1)-2]=p[ii(l)][ii(m)][ii(nz-1)];
				p[ii(l)][ii(m)][ii(1)-1]=p[ii(l)][ii(m)][ii(nz)];
				p[ii(l)][ii(m)][ii(nz)+1]=p[ii(l)][ii(m)][ii(1)];
				p[ii(l)][ii(m)][ii(nz)+2]=p[ii(l)][ii(m)][ii(2)];
			}
		}
		//edge
		//i-direction (i varying)
		if (l<=nx) {
			p[ii(l)][ii(1)-1][ii(1)-1]=p[ii(l)][ii(ny)][ii(nz)];
			p[ii(l)][ii(1)-1][ii(nz)+1]=p[ii(l)][ii(ny)][ii(1)];
			p[ii(l)][ii(ny)+1][ii(1)-1]=p[ii(l)][ii(1)][ii(nz)];
			p[ii(l)][ii(ny)+1][ii(nz)+1]=p[ii(l)][ii(1)][ii(1)];
		}
		//j-direction (j varying)
		if (l<=ny) {
			p[ii(1)-1][ii(l)][ii(1)-1]=p[ii(nx)][ii(l)][ii(nz)];
			p[ii(1)-1][ii(l)][ii(nz)+1]=p[ii(nx)][ii(l)][ii(1)];
			p[ii(nx)+1][ii(l)][ii(1)-1]=p[ii(1)][ii(l)][ii(nz)];
			p[ii(nx)+1][ii(l)][ii(nz)+1]=p[ii(1)][ii(l)][ii(1)];
		}
		//k-direction (k varying)
		if (l<=nz) {
			p[ii(1)-1][ii(1)-1][ii(l)]=p[ii(nx)][ii(ny)][ii(l)];
			p[ii(1)-1][ii(ny)+1][ii(l)]=p[ii(nx)][ii(1)][ii(l)];
			p[ii(nx)+1][ii(1)-1][ii(l)]=p[ii(1)][ii(ny)][ii(l)];
			p[ii(nx)+1][ii(ny)+1][ii(l)]=p[ii(1)][ii(1)][ii(l)];
		}
	}
	//coners
	p[ii(1)-1][ii(1)-1][ii(1)-1]=p[ii(nx)][ii(ny)][ii(1)-1];
	p[ii(1)-1][ii(ny)+1][ii(1)-1]=p[ii(nx)][ii(1)][ii(1)-1];
	p[ii(nx)+1][ii(1)-1][ii(1)-1]=p[ii(1)][ii(ny)][ii(1)-1];
	p[ii(nx)+1][ii(ny)+1][ii(1)-1]=p[ii(1)][ii(1)][ii(1)-1];

	p[ii(1)-1][ii(1)-1][ii(nz)+1]=p[ii(nx)][ii(ny)][ii(nz)+1];
	p[ii(1)-1][ii(ny)+1][ii(nz)+1]=p[ii(nx)][ii(1)][ii(nz)+1];
	p[ii(nx)+1][ii(1)-1][ii(nz)+1]=p[ii(1)][ii(ny)][ii(nz)+1];
	p[ii(nx)+1][ii(ny)+1][ii(nz)+1]=p[ii(1)][ii(1)][ii(nz)+1];
}

void BC_userdefine3D(double p[nxtot][nytot][nztot])
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
				p[ii(1)-2][ii(l)][ii(m)]=p[ii(nx-1)][ii(l)][ii(m)];
				p[ii(1)-1][ii(l)][ii(m)]=p[ii(nx)][ii(l)][ii(m)];
				p[ii(nx)+1][ii(l)][ii(m)]=p[ii(1)][ii(l)][ii(m)];
				p[ii(nx)+2][ii(l)][ii(m)]=p[ii(2)][ii(l)][ii(m)];
			}
			//j-layer (k,i varying)
			if (l<=nz && m <= nx) {
				p[ii(m)][ii(1)-2][ii(l)]=p[ii(m)][ii(ny-1)][ii(l)];
				p[ii(m)][ii(1)-1][ii(l)]=p[ii(m)][ii(ny)][ii(l)];
				p[ii(m)][ii(ny)+1][ii(l)]=p[ii(m)][ii(1)][ii(l)];
				p[ii(m)][ii(ny)+2][ii(l)]=p[ii(m)][ii(2)][ii(l)];
			}
			//k-layer (i,j varying)
			if (l<=nx && m <= ny) {
				p[ii(l)][ii(m)][ii(1)-2]=p[ii(l)][ii(m)][ii(2)];
				p[ii(l)][ii(m)][ii(1)-1]=p[ii(l)][ii(m)][ii(1)];
				p[ii(l)][ii(m)][ii(nz)+1]=p[ii(l)][ii(m)][ii(nz)];
				p[ii(l)][ii(m)][ii(nz)+2]=p[ii(l)][ii(m)][ii(nz-1)];
			}
		}
		//edge
		//i-direction (i varying)
		if (l<=nx) {
			p[ii(l)][ii(1)-1][ii(1)-1]=p[ii(l)][ii(ny)][ii(nz)];
			p[ii(l)][ii(1)-1][ii(nz)+1]=p[ii(l)][ii(ny)][ii(1)];
			p[ii(l)][ii(ny)+1][ii(1)-1]=p[ii(l)][ii(1)][ii(nz)];
			p[ii(l)][ii(ny)+1][ii(nz)+1]=p[ii(l)][ii(1)][ii(1)];
		}
		//j-direction (j varying)
		if (l<=ny) {
			p[ii(1)-1][ii(l)][ii(1)-1]=p[ii(nx)][ii(l)][ii(nz)];
			p[ii(1)-1][ii(l)][ii(nz)+1]=p[ii(nx)][ii(l)][ii(1)];
			p[ii(nx)+1][ii(l)][ii(1)-1]=p[ii(1)][ii(l)][ii(nz)];
			p[ii(nx)+1][ii(l)][ii(nz)+1]=p[ii(1)][ii(l)][ii(1)];
		}
		//k-direction (k varying)
		if (l<=nz) {
			p[ii(1)-1][ii(1)-1][ii(l)]=p[ii(nx)][ii(ny)][ii(l)];
			p[ii(1)-1][ii(ny)+1][ii(l)]=p[ii(nx)][ii(1)][ii(l)];
			p[ii(nx)+1][ii(1)-1][ii(l)]=p[ii(1)][ii(ny)][ii(l)];
			p[ii(nx)+1][ii(ny)+1][ii(l)]=p[ii(1)][ii(1)][ii(l)];
		}
	}
	//coners
	p[ii(1)-1][ii(1)-1][ii(1)-1]=p[ii(nx)][ii(ny)][ii(1)-1];
	p[ii(1)-1][ii(ny)+1][ii(1)-1]=p[ii(nx)][ii(1)][ii(1)-1];
	p[ii(nx)+1][ii(1)-1][ii(1)-1]=p[ii(1)][ii(ny)][ii(1)-1];
	p[ii(nx)+1][ii(ny)+1][ii(1)-1]=p[ii(1)][ii(1)][ii(1)-1];

	p[ii(1)-1][ii(1)-1][ii(nz)+1]=p[ii(nx)][ii(ny)][ii(nz)+1];
	p[ii(1)-1][ii(ny)+1][ii(nz)+1]=p[ii(nx)][ii(1)][ii(nz)+1];
	p[ii(nx)+1][ii(1)-1][ii(nz)+1]=p[ii(1)][ii(ny)][ii(nz)+1];
	p[ii(nx)+1][ii(ny)+1][ii(nz)+1]=p[ii(1)][ii(1)][ii(nz)+1];
}
void solvPF_CH(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
	//////////////////
	void calc_dG(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void calc_EMM(double EM[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
	//////////////////
	int i,j,k;

	double dG[nphase][nphase];
	double dGxp[nphase][nphase];
	double dGxm[nphase][nphase];
	double dGyp[nphase][nphase];
	double dGym[nphase][nphase];
	double dGzp[nphase][nphase];
	double dGzm[nphase][nphase];

        double EMxp[nphase][nphase];
        double EMxm[nphase][nphase];
        double EMyp[nphase][nphase];
        double EMym[nphase][nphase];
        double EMzp[nphase][nphase];
        double EMzm[nphase][nphase];

	double gradGxp[nphase][nphase],gradGxm[nphase][nphase];
	double gradGyp[nphase][nphase],gradGym[nphase][nphase];
	double gradGzp[nphase][nphase],gradGzm[nphase][nphase];

	double dJx[nphase][nphase];
	double dJy[nphase][nphase];
	double dJz[nphase][nphase];

	double divJ[nphase][nphase];

	double delta[nphase]={0.};
	int lp,mp;
	//////////////////

	dt = pow(dxl,4.)/(12.*EPSMAX*EPSMAX*EMM)*0.1/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {

		calc_dG(dG,pold,i,j,k);
		calc_dG(dGxp,pold,i+1,j,k);
		calc_dG(dGxm,pold,i-1,j,k);
		calc_dG(dGyp,pold,i,j+1,k);
		calc_dG(dGym,pold,i,j-1,k);
		calc_dG(dGzp,pold,i,j,k+1);
		calc_dG(dGzm,pold,i,j,k-1);

                calc_EMM(EMxp,pold,i+1,j,k);
                calc_EMM(EMxm,pold,i-1,j,k);
                calc_EMM(EMyp,pold,i,j+1,k);
                calc_EMM(EMym,pold,i,j-1,k);
                calc_EMM(EMzp,pold,i,j,k+1);
                calc_EMM(EMzm,pold,i,j,k-1);

		for (lp=0;lp<nphase;lp++) {
		for (mp=0;mp<nphase;mp++) {
			if (lp!=mp) {
				gradGxp[lp][mp]=(dGxp[lp][mp]-dG[lp][mp])/dxl;
				gradGxm[lp][mp]=(dG[lp][mp]-dGxm[lp][mp])/dxl;
				gradGyp[lp][mp]=(dGyp[lp][mp]-dG[lp][mp])/dyl;
				gradGym[lp][mp]=(dG[lp][mp]-dGym[lp][mp])/dyl;
				gradGzp[lp][mp]=(dGzp[lp][mp]-dG[lp][mp])/dzl;
				gradGzm[lp][mp]=(dG[lp][mp]-dGzm[lp][mp])/dzl;

				dJx[lp][mp]=(EMxp[lp][mp]*gradGxp[lp][mp]-EMxm[lp][mp]*gradGxm[lp][mp])/dxl;
				dJy[lp][mp]=(EMyp[lp][mp]*gradGyp[lp][mp]-EMym[lp][mp]*gradGym[lp][mp])/dyl;
				dJz[lp][mp]=(EMzp[lp][mp]*gradGzp[lp][mp]-EMzm[lp][mp]*gradGzm[lp][mp])/dzl;

				divJ[lp][mp]=dJx[lp][mp]+dJy[lp][mp]+dJz[lp][mp];
			}
		}
		}

		for (lp=0;lp<nphase;lp++) {
			delta[lp]=0.;
			for (mp=0;mp<nphase;mp++) {
				if (lp!=mp) {
					delta[lp] += divJ[lp][mp];
				}
			}
			pnew[lp][ii(i)][ii(j)][ii(k)]=pold[lp][ii(i)][ii(j)][ii(k)]+delta[lp]*dt;
		}

	}
	}
	}

}
void relaxPF_CH(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
	//////////////////
	void calc_dG(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
	//////////////////
	int i,j,k;

	double dG[nphase][nphase];
	double dGxp[nphase][nphase];
	double dGxm[nphase][nphase];
	double dGyp[nphase][nphase];
	double dGym[nphase][nphase];
	double dGzp[nphase][nphase];
	double dGzm[nphase][nphase];

	double gradGxp[nphase][nphase],gradGxm[nphase][nphase];
	double gradGyp[nphase][nphase],gradGym[nphase][nphase];
	double gradGzp[nphase][nphase],gradGzm[nphase][nphase];

	double dJx[nphase][nphase];
	double dJy[nphase][nphase];
	double dJz[nphase][nphase];

	double divJ[nphase][nphase];

	double delta[nphase]={0.};
	int lp,mp;
	//////////////////

	dt = pow(dxl,4.)/(12.*EPSMAX*EPSMAX*EMM)*0.1/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {

		calc_dG(dG,pold,i,j,k);
		calc_dG(dGxp,pold,i+1,j,k);
		calc_dG(dGxm,pold,i-1,j,k);
		calc_dG(dGyp,pold,i,j+1,k);
		calc_dG(dGym,pold,i,j-1,k);
		calc_dG(dGzp,pold,i,j,k+1);
		calc_dG(dGzm,pold,i,j,k-1);

		for (lp=0;lp<nphase;lp++) {
		for (mp=0;mp<nphase;mp++) {
			if (lp!=mp) {
				gradGxp[lp][mp]=(dGxp[lp][mp]-dG[lp][mp])/dxl;
				gradGxm[lp][mp]=(dG[lp][mp]-dGxm[lp][mp])/dxl;
				gradGyp[lp][mp]=(dGyp[lp][mp]-dG[lp][mp])/dyl;
				gradGym[lp][mp]=(dG[lp][mp]-dGym[lp][mp])/dyl;
				gradGzp[lp][mp]=(dGzp[lp][mp]-dG[lp][mp])/dzl;
				gradGzm[lp][mp]=(dG[lp][mp]-dGzm[lp][mp])/dzl;

				dJx[lp][mp]=(EMMrlx(lp,mp)*gradGxp[lp][mp]-EMMrlx(lp,mp)*gradGxm[lp][mp])/dxl;
				dJy[lp][mp]=(EMMrlx(lp,mp)*gradGyp[lp][mp]-EMMrlx(lp,mp)*gradGym[lp][mp])/dyl;
				dJz[lp][mp]=(EMMrlx(lp,mp)*gradGzp[lp][mp]-EMMrlx(lp,mp)*gradGzm[lp][mp])/dzl;

				divJ[lp][mp]=dJx[lp][mp]+dJy[lp][mp]+dJz[lp][mp];
			}
		}
		}

		for (lp=0;lp<nphase;lp++) {
			delta[lp]=0.;
			for (mp=0;mp<nphase;mp++) {
				if (lp!=mp) {
					delta[lp] += divJ[lp][mp];
				}
			}
			pnew[lp][ii(i)][ii(j)][ii(k)]=pold[lp][ii(i)][ii(j)][ii(k)]+delta[lp]*dt;
		}

	}
	}
	}

}

void solvPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
	//////////////////
	void calc_dG(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void calc_EMM(double EM[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void trim_p(double p[nphase][nxtot][nytot][nztot],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

        double EM[nphase][nphase];
	double dG[nphase][nphase];
	double delta[nphase]={0.};
        double p_l,p_m;
	int lp,mp;
        int npw = nphase;
	//////////////////

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
		calc_EMM(EM,pold,i,j,k);
		calc_dG(dG,pold,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ii(i)][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)+1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)-1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)+1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)-1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)+1];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)-1];
		    delta[lp]=0.;
                    if (p_l>0.+0.01 && p_l<7.-0.01) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ii(i)][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)+1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)-1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)+1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)-1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)+1];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)-1];
                            if (p_m>0.+0.01 && p_m<7.-0.01) {
				if (lp!=mp) {
					delta[lp] += -EM[lp][mp]*dG[lp][mp];
				}
                            }
			}
                    }
                    //RR[lp]+=delta[lp]/(-EMM);
		    pnew[lp][ii(i)][ii(j)][ii(k)]=pold[lp][ii(i)][ii(j)][ii(k)]+delta[lp]*dt;
                    //only for parabolic double well potential
		}
                trim_p(pnew,ii(i),ii(j),ii(k));

	}
	}
	}
}

void drivPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
	//////////////////
	void calc_dDRIV(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void calc_EMM(double EM[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void trim_p(double p[nphase][nxtot][nytot][nztot],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

        double EM[nphase][nphase];
	double dG[nphase][nphase];
	double delta[nphase]={0.};
        double p_l,p_m;
	int lp,mp;
        int npw = nphase;
	//////////////////

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;

		calc_EMM(EM,pold,i,j,k);
		calc_dDRIV(dG,pold,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ii(i)][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)+1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)-1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)+1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)-1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)+1];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)-1];
		    delta[lp]=0.;
                    //if (p_l > 0.01 && p_l < 7.- 0.01) {
                    if (1) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ii(i)][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)+1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)-1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)+1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)-1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)+1];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)-1];
                            //if (p_m > 0.01 && p_m < 7.- 0.01) {
                            if (1) {
				if (lp!=mp) {
					delta[lp] += -EM[lp][mp]*dG[lp][mp];
				}
                            }
			}
                    }
		    pnew[lp][ii(i)][ii(j)][ii(k)]=pold[lp][ii(i)][ii(j)][ii(k)]+delta[lp]*dt;
                    //if (pnew[lp][ii(i)][ii(j)][ii(k)]<=0.) pnew[lp][ii(i)][ii(j)][ii(k)]=0.;
                    //if (pnew[lp][ii(i)][ii(j)][ii(k)]>=1.) pnew[lp][ii(i)][ii(j)][ii(k)]=1.;
		}
                //trim_p(pnew,ii(i),ii(j),ii(k));

	}
	}
	}
}
void relaxPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
	//////////////////
	void calc_dG_rlx(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k);
        void trim_p_rlx(double p[nphase][nxtot][nytot][nztot],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

	double dG[nphase][nphase];
	double delta[nphase]={0.};
        double p_l,p_m;
	int lp,mp;
        int npw = nphase;
	//////////////////

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
		calc_dG_rlx(dG,pold,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ii(i)][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)+1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)-1][ii(j)][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)+1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)-1][ii(k)];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)+1];
                    p_l += pold[lp][ii(i)][ii(j)][ii(k)-1];
		    delta[lp]=0.;
                    if (p_l>0.+0.01 && p_l<7.-0.01) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ii(i)][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)+1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)-1][ii(j)][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)+1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)-1][ii(k)];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)+1];
                            p_m += pold[mp][ii(i)][ii(j)][ii(k)-1];
                            if (p_m>0.+0.01 && p_m<7.-0.01) {
				if (lp!=mp) {
					delta[lp] += -EMMrlx(lp,mp)*dG[lp][mp];
				}
                            }
			}
                    }
	            pnew[lp][ii(i)][ii(j)][ii(k)]=pold[lp][ii(i)][ii(j)][ii(k)]+delta[lp]*dt;
		}
                trim_p_rlx(pnew,ii(i),ii(j),ii(k));

	}
	}
	}

}

void fileout_p(double p[nphase][nxtot][nytot][nztot])
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
		fprintf(out,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",i,j,k,p[0][ii(i)][ii(j)][ii(k)],p[1][ii(i)][ii(j)][ii(k)],p[2][ii(i)][ii(j)][ii(k)]);
		//fprintf(out,"%d\t%d\t%d\t%lf\t%lf\n",i,j,k,p[0][ii(i)][ii(j)][ii(k)],p[1][ii(i)][ii(j)][ii(k)]);
	}
	}
	}
	fclose(out);
}

double calc_frac_p(double p[nxtot][nytot][nztot])
{
	/////////////
	double sum=0.;
	int i,j,k;
	/////////////
	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
		sum+=p[ii(i)][ii(j)][ii(k)];
	}
	}
	}

	return sum/(double)(nx*ny*nz);
}

double f(double p1, double p2, double p3)
{
    /////
    double h(double);
    /////
    extern double lamb[nphase];
    double result;
    /////
    //result = lamb[0]*h(p1)+lamb[1]*h(p2)+lamb[2]*h(p3);
    result = lamb[0]*h(p1)+lamb[1]*h(p2);
    return result;
}
double df(int np, double p1[nphase][nxtot][nytot][nztot],int iw,int jw,int kw)
{
    /////
    double dhmn(double phi1,double phi2);
    double dh(double);
    /////
    double result=0.;
    int nn;
    /////
    //for (nn=0;nn<nphase;nn++) {
    //    if (nn!=np) {
    //        result += dhmn(p1[np][iw][jw][kw],p1[nn][iw][jw][kw]);
    //    }
    //} 
    result = dh(p1[np][iw][jw][kw]);
    return result;
}

double h(double phi)
{
    double result;
    if (phi >= 1.) {
        result = 1.;
    }else if (phi<=0.) {
        result = 0.;
    }else{
        result = phi*phi*(3.-2.*phi);
    }
    return result;
}
double dh(double phi)
{
    double result;
    if (phi >= 1.) {
        result = 0.;
    }else if (phi<=0.) {
        result = 0.;
    }else{
        result = 6.*phi*(1.-phi);
    }
    return result;
}

double dhmn(double phi1,double phi2)
{
    //////////////
    double result;
    //////////////
    result = phi1*phi2;

    return result;
}

double g(double phi1,double phi2)
{
    if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0) {
        //return -phi1*phi2;
        return -phi1*phi2;
    }else{
        return phi1*phi2;
    }
}

double dg(double phi1, double phi2)
{
        if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
	    //return -phi2;
	    return 0.;
        }else{
	    return phi2;
        }
}

double tg(double phi1, double phi2, double phi3)
{
    if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
        //return -phi1*phi2*phi3;
        return phi1*phi1*phi2*phi2*phi3*phi3;
    }else{
        return phi1*phi2*phi3;
    }
}

double dtg(double phi1, double phi2, double phi3)
{
    return phi2*phi3;
}

double g_rlx(double phi1,double phi2)
{
    if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
        //return -phi1*phi2;
        return phi1*phi1*phi2*phi2;
    }else{
        return phi1*phi2;
    }
}
double dg_rlx(double phi1, double phi2)
{
        if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
	    //return -phi2;
	    return 0.;
        }else{
	    return phi2;
        }
}

int max3(int a1, int a2, int a3) 
{
	int result = a1;
	if (a2 >= result) {
		result = a2;
	}
	if (a3 >= result) {
		result = a3;
	}
	return result;
}

void calc_EMM(double EM[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k)
{
	/////////////
	int iw,jw,kw;
	int lp,mp,np;
        int npw = nphase;
	/////////////
	iw=ii(i);
	jw=ii(j);
	kw=ii(k);
        //p[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
	for (lp=0;lp<nphase;lp++) {
        for (mp=0;mp<nphase;mp++) {
             if (lp!=mp) {
                 //EM[lp][mp]=EMMp(0,1)*p1 + EMMp(1,0)*p2 + EMMp(0,2)*p3;
                 EM[lp][mp]=EMMp(lp,mp);
             }      
        }
	}
}
void calc_dG(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k)
{
        /////////////
        int makephaselist(int plist[nphase],double p[nphase][nxtot][nytot][nztot],int iw, int jw, int kw);
        /////////////
        extern double RR[nphase];
        extern double HH[nphase];
	/////////////
	int iw,jw,kw;
	double pxx[nphase],pyy[nphase],pzz[nphase];
	double lapl[nphase];
        double dFdphi[nphase];
	int lp,mp,np;
        int npw=nphase;
        int plist[nphase];
	/////////////
	iw=ii(i);
	jw=ii(j);
	kw=ii(k);
        //p[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
        npw=makephaselist(plist,p,iw,jw,kw);
	for (lp=0;lp<nphase;lp++) {
		pxx[lp] = (p[lp][iw+1][jw][kw]-2.*p[lp][iw][jw][kw]+p[lp][iw-1][jw][kw])/(dxl*dxl);
		pyy[lp] = (p[lp][iw][jw+1][kw]-2.*p[lp][iw][jw][kw]+p[lp][iw][jw-1][kw])/(dyl*dyl);
		pzz[lp] = (p[lp][iw][jw][kw+1]-2.*p[lp][iw][jw][kw]+p[lp][iw][jw][kw-1])/(dzl*dzl);
		lapl[lp] = pxx[lp] + pyy[lp] + pzz[lp];
        }
	for (lp=0;lp<nphase;lp++) {
                dFdphi[lp] = 0.;
                for (mp=0;mp<nphase;mp++) {
                    if (lp!=mp) {
                        dFdphi[lp] += EPS(lp,mp)*EPS(lp,mp)*lapl[mp]/2.;
                        dFdphi[lp] += W(lp,mp)*dg(p[lp][iw][jw][kw],p[mp][iw][jw][kw]);
                    }
                    for (np=0;np<nphase;np++) {
                       if (lp!=mp && mp!=np && np!=lp) {
                           dFdphi[lp] += 3.*WSV*dtg(p[lp][iw][jw][kw],p[mp][iw][jw][kw],p[np][iw][jw][kw]);
                       }
                    }
                    
                }
                dFdphi[lp]*=(double)plist[lp];
                //RR[lp]+=dFdphi[lp];
                HH[lp]+=df(lp,p,iw,jw,kw);
                HH[lp]*=(double)plist[lp];
	}
	for (lp=0;lp<nphase;lp++) {
	for (mp=0;mp<nphase;mp++) {
		if (lp!=mp) {
			dG[lp][mp]=dFdphi[lp]-dFdphi[mp];
			dG[lp][mp]*=(1./(double)npw);
                        //dG[lp][mp]*=s[lp][iw][jw][kw]*s[mp][iw][jw][kw];
		}
	}
	}
}
void calc_dDRIV(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k)
{
        int makephaselist(int plist[nphase],double p[nphase][nxtot][nytot][nztot],int iw, int jw, int kw);
	/////////////
        extern double lamb[nphase];

	int iw,jw,kw;
        double dFdphi[nphase];
	int lp,mp,np;
        int npw = nphase;
        int p_list[nphase];
	/////////////
	iw=ii(i);
	jw=ii(j);
	kw=ii(k);
        npw=makephaselist(p_list,p,iw,jw,kw);
	for (lp=0;lp<nphase;lp++) {
            dFdphi[lp] = lamb[lp]*df(lp,p,iw,jw,kw);
            dFdphi[lp] *= (double)p_list[lp];
	}
	for (lp=0;lp<nphase;lp++) {
	for (mp=0;mp<nphase;mp++) {
		if (lp!=mp) {
			dG[lp][mp]=dFdphi[lp]-dFdphi[mp];
			dG[lp][mp]*=(1./(double)npw);
                        //dG[lp][mp]*=s[lp][iw][jw][kw]*s[mp][iw][jw][kw];
		}
	}
	}
}

void calc_dG_rlx(double dG[nphase][nphase], double p[nphase][nxtot][nytot][nztot],int i,int j,int k)
{
        /////////////
        int makephaselist(int plist[nphase],double p[nphase][nxtot][nytot][nztot],int iw, int jw, int kw);
        double g_rlx(double, double);
        double dg_rlx(double,double);
	/////////////
	int iw,jw,kw;
	double pxx[nphase],pyy[nphase],pzz[nphase];
	double lapl[nphase];
        double dFdphi[nphase];
	int lp,mp,np;
        int npw = nphase;
        int plist[nphase];

        double SIGrlx=1.;
        double Wrlx=SIGrlx/2./XI;
        double EPSrlx=(2./pi)*sqrt(SIGrlx*XI);
	/////////////
	iw=ii(i);
	jw=ii(j);
	kw=ii(k);
        //p[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
        npw=makephaselist(plist,p,iw,jw,kw);
	for (lp=0;lp<nphase;lp++) {
		pxx[lp] = (p[lp][iw+1][jw][kw]-2.*p[lp][iw][jw][kw]+p[lp][iw-1][jw][kw])/(dxl*dxl);
		pyy[lp] = (p[lp][iw][jw+1][kw]-2.*p[lp][iw][jw][kw]+p[lp][iw][jw-1][kw])/(dyl*dyl);
		pzz[lp] = (p[lp][iw][jw][kw+1]-2.*p[lp][iw][jw][kw]+p[lp][iw][jw][kw-1])/(dzl*dzl);
		lapl[lp] = pxx[lp] + pyy[lp] + pzz[lp];
        }
	for (lp=0;lp<nphase;lp++) {
                dFdphi[lp] = 0.;
                for (mp=0;mp<nphase;mp++) {
                    if (lp!=mp) {
                        dFdphi[lp] += EPSrlx*EPSrlx*lapl[mp]/2.+Wrlx*dg_rlx(p[lp][iw][jw][kw],p[mp][iw][jw][kw]);
                    }
                }
	}
        for (lp=0;lp<nphase;lp++) {
            for (mp=0;mp<nphase;mp++) {
            for (np=0;np<nphase;np++) {
                if (lp<mp && mp<np) {
                    //dFdphi[lp] += 6.*W(mp,np)*dtg(p[lp][iw][jw][kw],p[mp][iw][jw][kw],p[np][iw][jw][kw]);
                }
            }
            }
            dFdphi[lp]*=(double)plist[lp];
        }
	for (lp=0;lp<nphase;lp++) {
	for (mp=0;mp<nphase;mp++) {
		if (lp!=mp) {
			dG[lp][mp]=dFdphi[lp]-dFdphi[mp];
			dG[lp][mp]*=(1./(double)npw);
                        //dG[lp][mp]*=s[lp][iw][jw][kw]*s[mp][iw][jw][kw];
		}
	}
	}
}

void calc_s(double s[nphase][nxtot][nytot][nztot],double p[nphase][nxtot][nytot][nztot])
{
    /////////////////
    int i,j,k,np;
    /////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        for (np=0;np<nphase;np++) {
            if (p[np][ii(i)][ii(j)][ii(k)]>0.) {
            //if (np!=2) {
                s[np][ii(i)][ii(j)][ii(k)]=1.;
            }else{
                s[np][ii(i)][ii(j)][ii(k)]=0.;
            }
        }
    }
    }
    }
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

void frac_out(double p[nphase][nxtot][nytot][nztot])
{
    /////////////////
    double calc_frac_p(double p[nxtot][nytot][nztot]);
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

void draw_cube(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if (abs(i-ix)<=size/2 && abs(j-iy)<=size/2 && abs(k-iz)<=size/2) {
            p[ii(i)][ii(j)][ii(k)]=1.;
        }
    }
    }
    }

}

void draw_sphere(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int size)
{
    // make a cube at position (i,j,k)
    //////////////////////////////////
    int i,j,k;
    //////////////////////////////////
    
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        if ((i-ix)*(i-ix)+(j-iy)*(j-iy)+(k-iz)*(k-iz) <= (size/2)*(size/2)) {
            p[ii(i)][ii(j)][ii(k)]=1.;
        }
    }
    }
    }

}
void draw_cylinder(double p[nxtot][nytot][nztot],int ix, int iy, int iz, char axis, int rad, int height)
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
                p[ii(i)][ii(j)][ii(k)]=1.;
            }
        }else if(axis=='x'){
            nr = (j-iy)*(j-iy)+(k-iz)*(k-iz);
            if (nr<=rad*rad && abs(i-ix)<=height/2) {
                p[ii(i)][ii(j)][ii(k)]=1.;
            }
        }else if(axis=='y') {
            nr = (k-iz)*(k-iz)+(i-ix)*(i-ix);
            if (nr<=rad*rad && abs(j-iy)<=height/2) {
                p[ii(i)][ii(j)][ii(k)]=1.;
            }
        }else{
            printf("draw_cylinder(): Wrong axis.\n");
        }
    }
    }
    }

}

void draw_rectangular(double p[nxtot][nytot][nztot],int ix, int iy, int iz, int sizex, int sizey,int sizez)
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
            p[ii(i)][ii(j)][ii(k)]=1.;
        }
    }
    }
    }

}

void fill_environ(int phaseindex,double p[nphase][nxtot][nytot][nztot])
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
                 sum += p[np][ii(i)][ii(j)][ii(k)];
            }
        }
        p[phaseindex][ii(i)][ii(j)][ii(k)] = 1.-sum;
    }
    }
    }
}
void fill_rectangular(int phaseindex,double p[nphase][nxtot][nytot][nztot],int ix, int iy, int iz, int sizex, int sizey, int sizez)
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
                   sum += p[np][ii(i)][ii(j)][ii(k)];
                }
            }
            p[phaseindex][ii(i)][ii(j)][ii(k)] = 1.-sum;
        }
    }
    }
    }
}
void fill_cylinder(int phaseindex,double p[nphase][nxtot][nytot][nztot],int ix, int iy, int iz, char axis, int rad, int height)
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
                       sum += p[np][ii(i)][ii(j)][ii(k)];
                    }
                }
                p[phaseindex][ii(i)][ii(j)][ii(k)] = 1.-sum;
            }
        } else if (axis=='x') {
            nr = (k-iz)*(k-iz)+(j-iy)*(j-iy);
            if (nr<=rad*rad && abs(i-ix)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ii(i)][ii(j)][ii(k)];
                    }
                }
                p[phaseindex][ii(i)][ii(j)][ii(k)] = 1.-sum;
            }
        } else if (axis=='y') {
            nr = (i-ix)*(i-ix)+(k-iz)*(k-iz);
            if (nr<=rad*rad && abs(j-iy)<=height/2) {
                sum=0.;
                for (np=0;np<nphase;np++) {
                    if (np!=phaseindex) {
                       sum += p[np][ii(i)][ii(j)][ii(k)];
                    }
                }
                p[phaseindex][ii(i)][ii(j)][ii(k)] = 1.-sum;
            }
        }else{
            printf("fill_cylinder(): Wrong axis.\n");
        }
    }
    }
    }
}
void clean_all(double p[nphase][nxtot][nytot][nztot])
{
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        for (np=0;np<nphase;np++) {
            p[np][ii(i)][ii(j)][ii(k)]=0.;
        }
    }
    }
    }
    
}

void clean_one(double p[nxtot][nytot][nztot])
{
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        p[ii(i)][ii(j)][ii(k)]=0.;
    }
    }
    }
    
}

void relax_p(double p[nphase][nxtot][nytot][nztot], int nrelax)
{
    ///////////////////
    void relaxPF_CH(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
    void relaxPF_CA(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot]);
    void copy_p(double presult[nphase][nxtot][nytot][nztot], double psource[nphase][nxtot][nytot][nztot]);
    void update_p(double pold[nphase][nxtot][nytot][nztot], double pnew[nphase][nxtot][nytot][nztot]);
    void BC_periodic3D(double p[nxtot][nytot][nztot]);
    void BC_adiabatic3D(double p[nxtot][nytot][nztot]);
    ///////////////////
    double pdum[nphase][nxtot][nytot][nztot];
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

double EMMp(int np1, int np2)
{
    ///////////////////
    extern int pfix;
    //
    double result;
    //double p1,p2,p3;
    ///////////////////
    
    if (np1==pfix || np2==pfix) {
        result = 0.;
    }else{
        result = EMM;
    }

    //result = EMM; 
    return result;
}

double EMMrlx(int np1, int np2)
{
    ///////////////////
    double result;
    ///////////////////
    result = EMM;
    return result;
}
void copy_p(double presult[nphase][nxtot][nytot][nztot], double psource[nphase][nxtot][nytot][nztot])
{
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        for (np=0;np<nphase;np++) {
            presult[np][ii(i)][ii(j)][ii(k)]=psource[np][ii(i)][ii(j)][ii(k)];
        }
    }
    }
    }
}

double SIG(int np1, int np2)
{
    //////////////
    double result;
    //////////////
    if ((np1==0 && np2==1) || (np1==1 && np2==0)) {
        result = SIGLV;
    }
    if ((np1==0 && np2==2) || (np1==2 && np2==0)) {
        result = SIGSV;
    }
    if ((np1==1 && np2==2) || (np1==2 && np2==1)) {
        result = SIGSL;
    }
    if (np1==np2) {
        result = 0.;
    }
    
    return result;
}

double EPS(int np1, int np2)
{
    //////////////
    double result;
    //////////////
    if ((np1==0 && np2==1) || (np1==1 && np2==0)) {
        result = EPSLV;
    }
    if ((np1==0 && np2==2) || (np1==2 && np2==0)) {
        result = EPSSV;
    }
    if ((np1==1 && np2==2) || (np1==2 && np2==1)) {
        result = EPSSL;
    }
    if (np1==np2) {
        result = 0.;
    }
    
    return result;
}

double W(int np1, int np2)
{
    //////////////
    double result;
    //////////////
    if ((np1==0 && np2==1) || (np1==1 && np2==0)) {
        result = WLV;
    }
    if ((np1==0 && np2==2) || (np1==2 && np2==0)) {
        result = WSV;
    }
    if ((np1==1 && np2==2) || (np1==2 && np2==1)) {
        result = WSL;
    }
    if (np1==np2) {
        result = 0.;
    }
    
    return result;
}
double dmax3(double a1, double a2, double a3) 
{
	double result = a1;
	if (a2 >= result) {
		result = a2;
	}
	if (a3 >= result) {
		result = a3;
	}
	return result;
}

void calc_lamb(double pnew[nphase][nxtot][nytot][nztot], double pold[nphase][nxtot][nytot][nztot])
{
    //////////////
    double df(int,double p[nphase][nxtot][nytot][nztot],int,int,int);
    //////////////
    extern double lamb[nphase];
    extern double V0[nphase];
    extern double RR[nphase];
    extern double HH[nphase];
    extern double dxl,dyl,dzl;
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
            if (itime==0)V0[np]+=pold[np][iw][jw][kw];
            R[np]+=(pnew[np][iw][jw][kw])/dt;
            H[np]+=df(np,pnew,iw,jw,kw);
            if(pfix>=0 && pfix<nphase && fabs(pnew[pfix][iw][jw][kw])>0.) {
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

int makephaselist(int plist[nphase],double p[nphase][nxtot][nytot][nztot],int iw, int jw, int kw)
{
    ////////////
    extern double dxl, dyl, dzl;
    ////////////
    int npw=0;
    int np;
    double lap_p[nphase];
    ////////////
    for (np=0;np<nphase;np++) {
        lap_p[np]=(p[np][iw+1][jw][kw]-2.*p[np][iw][jw][kw]+p[np][iw-1][jw][kw])/dxl/dxl;
        lap_p[np]+=(p[np][iw][jw+1][kw]-2.*p[np][iw][jw][kw]+p[np][iw][jw-1][kw])/dyl/dyl;
        lap_p[np]+=(p[np][iw][jw][kw+1]-2.*p[np][iw][jw][kw]+p[np][iw][jw][kw-1])/dzl/dzl;
        if (lap_p[np]!=0.||p[np][iw][jw][kw]>0.) {
            plist[np]=1;
        }else{
            plist[np]=0;
        }
        npw += plist[np];
    }
    //if(npw==3)printf("three npw: %d %d %d %d\n",itime,iw,jw,kw);
    return npw;
}

void trim_p(double p[nphase][nxtot][nytot][nztot],int iw,int jw,int kw)
{
    /////////
    extern int pfix;    

    double psum=0.;
    int np;
    int npw=nphase; 
    int plist[nphase];
    /////////
    npw=makephaselist(plist,p,iw,jw,kw);
    for (np=0;np<nphase;np++) {
        if(p[np][iw][jw][kw]<=0.)p[np][iw][jw][kw]=0.;
        if(p[np][iw][jw][kw]>=1.)p[np][iw][jw][kw]=1.;
        psum+=p[np][iw][jw][kw];
    }
    if(pfix>=0 && pfix<nphase) { 
        psum -= p[pfix][iw][jw][kw];
    }
    for (np=0;np<nphase;np++) {
        if (pfix>=0 && pfix<nphase) {
            if(np!=pfix && psum!=0.)p[np][iw][jw][kw]*=(1.-p[pfix][iw][jw][kw])/psum;
        }else{
            p[np][iw][jw][kw]/=psum;
        }
    }
    
}

void trim_p_rlx(double p[nphase][nxtot][nytot][nztot],int iw,int jw,int kw)
{
    /////////
    double psum=0.;
    int np;
    int npw=nphase; 
    int plist[nphase];
    /////////
    npw=makephaselist(plist,p,iw,jw,kw);
    for (np=0;np<nphase;np++) {
        if(p[np][iw][jw][kw]<=0.)p[np][iw][jw][kw]=0.;
        if(p[np][iw][jw][kw]>=1.)p[np][iw][jw][kw]=1.;
        psum+=p[np][iw][jw][kw];
    }
    for (np=0;np<nphase;np++) {
        p[np][iw][jw][kw]/=psum;
    }
    
}
