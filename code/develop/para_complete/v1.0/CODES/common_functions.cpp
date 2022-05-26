#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "defined_values.h"
#include "common_functions.h"

extern int nx,ny,nz,ntot,nmax;

int pfix =  2;

double dxl = 1.0e-8;
double dyl = dxl;
double dzl = dxl;
double XI = 2.5*dxl; // half interface thickness

//double SIG = 1.;
double SIGSV = 1.;
double SIGLV = 1.;
double SIGSL = 0.2014;
//double SIGSL = 1.;

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
double dtime;
double time = 0.;
int itime = 0;
int dtkind;

double RR[nphase]={0.};
double V0[nphase]={0.};
double HH[nphase]={0.};
double lamb[nphase]; //interphase driving force coefficient


double calc_BoltzmannF(double *BF, double *p, double *volt, double z, double TEMP)

{
    //////////////
    extern int nx,ny,nz;
    extern double dxl,dyl,dzl;
    //////////////
    int i,j,k;
    int iw,jw,kw;
    double total=0.;
    //////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        iw=ii(i);
        jw=ii(j);
        kw=ii(k);
        BF[ic(iw,jw,kw)]=exp(-z*FARAD/(GasConst*TEMP)*volt[ic(iw,jw,kw)]);
        total+=p[ic(iw,jw,kw)]*BF[ic(iw,jw,kw)];
    }
    }
    }
    //total*=dxl*dyl*dzl;
    return total;
}

double df(int np, double *p1[nphase],int iw,int jw,int kw)
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
    result = dh(p1[np][ic(iw,jw,kw)]);
    return result;
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
double dg_rlx(double phi1, double phi2)
{
        if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
	    //return -phi2;
	    return 0.;
        }else{
	    return phi2;
        }
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
double dtg(double phi1, double phi2, double phi3)
{
    return phi2*phi3;
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
double g(double phi1,double phi2)
{
    if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0) {
        //return -phi1*phi2;
        return -phi1*phi2;
    }else{
        return phi1*phi2;
    }
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
double tg(double phi1, double phi2, double phi3)
{
    if (phi1>=1.||phi1<=0. || phi2>=1.||phi2<=0.) {
        //return -phi1*phi2*phi3;
        return phi1*phi1*phi2*phi2*phi3*phi3;
    }else{
        return phi1*phi2*phi3;
    }
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

int ic(int i, int j, int k)
{
    //coordinate mapping only for fields in real space [nx+4 by ny+4 by nz+4]
    //z-outermost order
    //same index order with FORTRAN 3d static arrays
    //advantageous on parallel data output
    extern int nx,ny,nz;
    //return (i)*(ny+4)*(nz+4)+(j)*(nz+4)+(k); //x-outermost
    return (i)+(nx+4)*(j)+(nx+4)*(ny+4)*k; //z-outermost
}
int icf(int i, int j, int k)
{
    //coordinate mapping only for field in reciprocal space(e.g. Fourier space) [nx by ny by nz]
    //z-outermost order
    extern int nx,ny,nz;
    //return (i)*(ny)*(nz)+(j)*(nz)+(k); //x-outermost
    return (i)+(nx)*(j)+(nx)*(ny)*k; //z-outermost
}
int makephaselist(int plist[nphase],double *p[nphase],int iw, int jw, int kw)
{
    ////////////
    extern double dxl, dyl, dzl;
    ////////////
    int npw=0;
    int np;
    double lap_p[nphase];
    ////////////
    for (np=0;np<nphase;np++) {
        lap_p[np]=(p[np][ic(iw+1,jw,kw)]-2.*p[np][ic(iw,jw,kw)]+p[np][ic(iw-1,jw,kw)])/dxl/dxl;
        lap_p[np]+=(p[np][ic(iw,jw+1,kw)]-2.*p[np][ic(iw,jw,kw)]+p[np][ic(iw,jw-1,kw)])/dyl/dyl;
        lap_p[np]+=(p[np][ic(iw,jw,kw+1)]-2.*p[np][ic(iw,jw,kw)]+p[np][ic(iw,jw,kw-1)])/dzl/dzl;
        if (lap_p[np]!=0.||p[np][ic(iw,jw,kw)]>0.) {
            plist[np]=1;
        }else{
            plist[np]=0;
        }
        npw += plist[np];
    }
    //if(npw==3)printf("three npw: %d %d %d %d\n",itime,iw,jw,kw);
    return npw;
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


void calc_dDRIV(double dG[nphase][nphase], double *p[nphase],double *volt, double *CHG,int i,int j,int k)
{
    int makephaselist(int plist[nphase],double *p[nphase],int iw, int jw, int kw);
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
        if(lp==1)dFdphi[lp] += 0.5*CHG[ic(iw,jw,kw)]*volt[ic(iw,jw,kw)]*df(1,p,iw,jw,kw);
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
void calc_dG(double dG[nphase][nphase], double *p[nphase],int i,int j,int k)
{
    /////////////
    int makephaselist(int plist[nphase],double *p[nphase],int iw, int jw, int kw);
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
		pxx[lp] = (p[lp][ic(iw+1,jw,kw)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw-1,jw,kw)])/(dxl*dxl);
		pyy[lp] = (p[lp][ic(iw,jw+1,kw)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw,jw-1,kw)])/(dyl*dyl);
		pzz[lp] = (p[lp][ic(iw,jw,kw+1)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw,jw,kw-1)])/(dzl*dzl);
		lapl[lp] = pxx[lp] + pyy[lp] + pzz[lp];
    }
	for (lp=0;lp<nphase;lp++) {
        dFdphi[lp] = 0.;
        for (mp=0;mp<nphase;mp++) {
            if (lp!=mp) {
                dFdphi[lp] += EPS(lp,mp)*EPS(lp,mp)*lapl[mp]/2.;
                dFdphi[lp] += W(lp,mp)*dg(p[lp][ic(iw,jw,kw)],p[mp][ic(iw,jw,kw)]);
            }
            for (np=0;np<nphase;np++) {
                if (lp!=mp && mp!=np && np!=lp) {
                    //dFdphi[lp] += 6.*WSV*dtg(p[lp][ic(iw,jw,kw)],p[mp][ic(iw,jw,kw)],p[np][ic(iw,jw,kw)]);
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
void calc_dG_rlx(double dG[nphase][nphase], double *p[nphase],int i,int j,int k)
{
    /////////////
    int makephaselist(int plist[nphase],double *p[nphase],int iw, int jw, int kw);
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
		pxx[lp] = (p[lp][ic(iw+1,jw,kw)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw-1,jw,kw)])/(dxl*dxl);
		pyy[lp] = (p[lp][ic(iw,jw+1,kw)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw,jw-1,kw)])/(dyl*dyl);
		pzz[lp] = (p[lp][ic(iw,jw,kw+1)]-2.*p[lp][ic(iw,jw,kw)]+p[lp][ic(iw,jw,kw-1)])/(dzl*dzl);
		lapl[lp] = pxx[lp] + pyy[lp] + pzz[lp];
    }
	for (lp=0;lp<nphase;lp++) {
        dFdphi[lp] = 0.;
        for (mp=0;mp<nphase;mp++) {
            if (lp!=mp) {
                dFdphi[lp] += EPSrlx*EPSrlx*lapl[mp]/2.+Wrlx*dg_rlx(p[lp][ic(iw,jw,kw)],p[mp][ic(iw,jw,kw)]);
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
void calc_EMM(double EM[nphase][nphase], double *p[nphase],int i,int j,int k)
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

void copy_array_d(double *out, double *in)
{
    //////////////////
    extern int nx,ny,ny,ntot;
    ///////////////////
    int i,j,k;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        out[ic(ii(i),ii(j),ii(k))]=in[ic(ii(i),ii(j),ii(k))];
    }
    }
    }
}
void copy_array_d_scale(double *out, double *in, double mag)
{
    //////////////////
    extern int nx,ny,ny,ntot;
    ///////////////////
    int i,j,k;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        out[ic(ii(i),ii(j),ii(k))]=in[ic(ii(i),ii(j),ii(k))]*mag;
    }
    }
    }
}

void copy_p(double *presult[nphase], double *psource[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
    ///////////////////
    int i,j,k,np;
    ///////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        for (np=0;np<nphase;np++) {
            presult[np][ic(ii(i),ii(j),ii(k))]=psource[np][ic(ii(i),ii(j),ii(k))];
        }
    }
    }
    }
}
void del_buff_cx(double complex **p)
{
    delete[] *p;
    *p = NULL;
}
void del_buff_d(double **p)
{
    delete[] *p;
    *p = NULL;
}

void del_buff_i(int **p)
{
    delete[] *p;
    *p = NULL;
}

void drivPF_CA(double *pnew[nphase], double *pold[nphase], double *volt, double *CHG)
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dDRIV(double dG[nphase][nphase], double *p[nphase],double *volt, double *CHG, int i,int j,int k);
        void calc_EMM(double EM[nphase][nphase], double *p[nphase],int i,int j,int k);
        void trim_p(double *p[nphase],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

        double EM[nphase][nphase];
	double dG[nphase][nphase];
	double delta[nphase]={0.};
        double p_l,p_m;
	int lp,mp;
        int npw = nphase;
	//////////////////

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;

		calc_EMM(EM,pold,i,j,k);
		calc_dDRIV(dG,pold,volt,CHG,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ic(ii(i),ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)+1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)-1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)+1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)-1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)+1)];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)-1)];
		    delta[lp]=0.;
                    //if (p_l > 0.01 && p_l < 7.- 0.01) {
                    if (1) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ic(ii(i),ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)+1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)-1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)+1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)-1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)+1)];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)-1)];
                            //if (p_m > 0.01 && p_m < 7.- 0.01) {
                            if (1) {
				if (lp!=mp) {
					delta[lp] += -EM[lp][mp]*dG[lp][mp];
				}
                            }
			}
                    }
		    pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dtime;
                    //if (pnew[lp][ii(i)][ii(j)][ii(k)]<=0.) pnew[lp][ii(i)][ii(j)][ii(k)]=0.;
                    //if (pnew[lp][ii(i)][ii(j)][ii(k)]>=1.) pnew[lp][ii(i)][ii(j)][ii(k)]=1.;
		}
                //trim_p(pnew,ii(i),ii(j),ii(k));

	}
	}
	}
}
void relaxPF_CA(double *pnew[nphase], double *pold[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dG_rlx(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
    void trim_p_rlx(double *p[nphase],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

	double dG[nphase][nphase];
	double delta[nphase]={0.};
        double p_l,p_m;
	int lp,mp;
        int npw = nphase;
        double dtRLX;
	//////////////////

        dtRLX = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;
 
	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;
		calc_dG_rlx(dG,pold,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ic(ii(i),ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)+1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)-1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)+1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)-1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)+1)];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)-1)];
		    delta[lp]=0.;
                    if (p_l>0.+0.01 && p_l<7.-0.01) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ic(ii(i),ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)+1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)-1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)+1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)-1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)+1)];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)-1)];
                            if (p_m>0.+0.01 && p_m<7.-0.01) {
				if (lp!=mp) {
					delta[lp] += -EMMrlx(lp,mp)*dG[lp][mp];
				}
                            }
			}
                    }
	            pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dtRLX;
		}
                trim_p_rlx(pnew,ii(i),ii(j),ii(k));

	}
	}
	}

}
void relaxPF_CH(double *pnew[nphase], double *pold[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dG(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
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
        double dtRLX;
	//////////////////

        dtRLX = pow(dxl,4.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;
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
			pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dtRLX;
		}

	}
	}
	}

}
#include <stdio.h>
double calc_Ce(double *Cnew, double *C, double *p[nphase], double *v[3])
{
        //////////////////
        extern int nx,ny,ny,ntot;
        extern double dxl,dyl,dzl;
        extern double z[ncomp];
        extern double DIFFUSIVITY;
        extern double dtime;
	//////////////////
	int i,j,k;
        int iw,jw,kw;
        double C0;
        double Cxp, Cyp, Czp;
        double Cxm, Cym, Czm;
        double Cxpyp,Cypzp,Czpxp;
        double Cxmym,Cymzm,Czmxm;
        double Cxpym,Cxmyp;
        double Cypzm,Cymzp;
        double Czpxm,Czmxp;
        double p0,p1,p2,p3,p4,p5,p6;
        double p13,p14,p15,p16;
        double p35,p36;
        double p23,p24,p25,p26;
        double p45,p46;

        double pxp,pyp,pzp;
        double pxm,pym,pzm;
        double pxpyp,pypzp,pzpxp;
        double pxmym,pymzm,pzmxm;
        double pxpym,pxmyp;
        double pypzm,pymzp;
        double pzpxm,pzmxp;
        double dpmag0;
        double dpmagxp,dpmagyp,dpmagzp;
        double dpmagxm,dpmagym,dpmagzm;
        double nvecxp_x,nvecxp_y,nvecxp_z;
        double nvecxm_x,nvecxm_y,nvecxm_z;
        double nvecyp_x,nvecyp_y,nvecyp_z;
        double nvecym_x,nvecym_y,nvecym_z;
        double nveczp_x,nveczp_y,nveczp_z;
        double nveczm_x,nveczm_y,nveczm_z;
        double nvecmagxp,nvecmagxm;
        double nvecmagyp,nvecmagym;
        double nvecmagzp,nvecmagzm;
        double vx0,vy0,vz0;
        double vxxp,vxyp,vxzp;
        double vxxm,vxym,vxzm;
        double vyxp,vyyp,vyzp;
        double vyxm,vyym,vyzm;
        double vzxp,vzyp,vzzp;
        double vzxm,vzym,vzzm;

        double Jxp, Jyp, Jzp;
        double Jxm, Jym, Jzm;
        double J0xp_x=0.,J0xp_y=0.,J0xp_z=0.;
        double J0xm_x=0.,J0xm_y=0.,J0xm_z=0.;
        double J0yp_x=0.,J0yp_y=0.,J0yp_z=0.;
        double J0ym_x=0.,J0ym_y=0.,J0ym_z=0.;
        double J0zp_x=0.,J0zp_y=0.,J0zp_z=0.;
        double J0zm_x=0.,J0zm_y=0.,J0zm_z=0.;
        double J0dotnxp, J0dotnyp, J0dotnzp;
        double J0dotnxm, J0dotnym, J0dotnzm;
        double JBxp=0.,JByp=0.,JBzp=0.;
        double JBxm=0.,JBym=0.,JBzm=0.;
	double divJ,divJ0,divJB;

        double diffu0 = DIFFUSIVITY;

        double Ctot=0.;
	//////////////////

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            iw=ii(i);
            jw=ii(j);
            kw=ii(k);
            C0 = C[ic(iw,jw,kw)];
            Cxp =0.5*(C[ic(iw+1,jw,kw)]+C[ic(iw,jw,kw)]);
            Cxm =0.5*(C[ic(iw,jw,kw)]+C[ic(iw-1,jw,kw)]);
            Cyp =0.5*(C[ic(iw,jw+1,kw)]+C[ic(iw,jw,kw)]);
            Cym =0.5*(C[ic(iw,jw,kw)]+C[ic(iw,jw-1,kw)]);
            Czp =0.5*(C[ic(iw,jw,kw+1)]+C[ic(iw,jw,kw)]);
            Czm =0.5*(C[ic(iw,jw,kw)]+C[ic(iw,jw,kw-1)]);

            Cxpyp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw+1,jw,kw)]+C[ic(iw,jw+1,kw)]+C[ic(iw+1,jw+1,kw)]);
            Cypzp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw+1,kw)]+C[ic(iw,jw,kw+1)]+C[ic(iw,jw+1,kw+1)]);
            Czpxp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw,kw+1)]+C[ic(iw+1,jw,kw)]+C[ic(iw+1,jw,kw+1)]);
            Cxmym=0.25*(C[ic(iw,jw,kw)]+C[ic(iw-1,jw,kw)]+C[ic(iw,jw-1,kw)]+C[ic(iw-1,jw-1,kw)]);
            Cymzm=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw-1,kw)]+C[ic(iw,jw,kw-1)]+C[ic(iw,jw-1,kw-1)]);
            Czmxm=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw,kw-1)]+C[ic(iw-1,jw,kw)]+C[ic(iw-1,jw,kw-1)]);
            Cxpym=0.25*(C[ic(iw,jw,kw)]+C[ic(iw+1,jw,kw)]+C[ic(iw,jw-1,kw)]+C[ic(iw+1,jw-1,kw)]);
            Cxmyp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw-1,jw,kw)]+C[ic(iw,jw+1,kw)]+C[ic(iw-1,jw+1,kw)]);
            Cypzm=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw+1,kw)]+C[ic(iw,jw,kw-1)]+C[ic(iw,jw+1,kw-1)]);
            Cymzp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw-1,kw)]+C[ic(iw,jw,kw+1)]+C[ic(iw,jw-1,kw+1)]);
            Czpxm=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw,kw+1)]+C[ic(iw-1,jw,kw)]+C[ic(iw-1,jw,kw+1)]);
            Czmxp=0.25*(C[ic(iw,jw,kw)]+C[ic(iw,jw,kw-1)]+C[ic(iw+1,jw,kw)]+C[ic(iw+1,jw,kw-1)]);
            vxxp =0.5*(v[0][ic(iw+1,jw,kw)]+v[0][ic(iw,jw,kw)]);
            vxxm =0.5*(v[0][ic(iw,jw,kw)]+v[0][ic(iw-1,jw,kw)]);
            vxyp =0.5*(v[0][ic(iw,jw+1,kw)]+v[0][ic(iw,jw,kw)]);
            vxym =0.5*(v[0][ic(iw,jw,kw)]+v[0][ic(iw,jw-1,kw)]);
            vxzp =0.5*(v[0][ic(iw,jw,kw+1)]+v[0][ic(iw,jw,kw)]);
            vxzm =0.5*(v[0][ic(iw,jw,kw)]+v[0][ic(iw,jw,kw-1)]);
            vyxp =0.5*(v[1][ic(iw+1,jw,kw)]+v[1][ic(iw,jw,kw)]);
            vyxm =0.5*(v[1][ic(iw,jw,kw)]+v[1][ic(iw-1,jw,kw)]);
            vyyp =0.5*(v[1][ic(iw,jw+1,kw)]+v[1][ic(iw,jw,kw)]);
            vyym =0.5*(v[1][ic(iw,jw,kw)]+v[1][ic(iw,jw-1,kw)]);
            vyzp =0.5*(v[1][ic(iw,jw,kw+1)]+v[1][ic(iw,jw,kw)]);
            vyzm =0.5*(v[1][ic(iw,jw,kw)]+v[1][ic(iw,jw,kw-1)]);
            vzxp =0.5*(v[2][ic(iw+1,jw,kw)]+v[2][ic(iw,jw,kw)]);
            vzxm =0.5*(v[2][ic(iw,jw,kw)]+v[2][ic(iw-1,jw,kw)]);
            vzyp =0.5*(v[2][ic(iw,jw+1,kw)]+v[2][ic(iw,jw,kw)]);
            vzym =0.5*(v[2][ic(iw,jw,kw)]+v[2][ic(iw,jw-1,kw)]);
            vzzp =0.5*(v[2][ic(iw,jw,kw+1)]+v[2][ic(iw,jw,kw)]);
            vzzm =0.5*(v[2][ic(iw,jw,kw)]+v[2][ic(iw,jw,kw-1)]);


            p0 = 0.;
            p1 = 0.;
            p2 = 0.;
            p3 = 0.;
            p4 = 0.;
            p5 = 0.;
            p6 = 0.;

            if(p[1][ic(iw,jw,kw)]>=1.0) p0=1.;
            if(p[1][ic(iw+1,jw,kw)]>=1.0) p1=1.;
            if(p[1][ic(iw-1,jw,kw)]>=1.0) p2=1.;
            if(p[1][ic(iw,jw+1,kw)]>=1.0) p3=1.;
            if(p[1][ic(iw,jw-1,kw)]>=1.0) p4=1.;
            if(p[1][ic(iw,jw,kw+1)]>=1.0) p5=1.;
            if(p[1][ic(iw,jw,kw-1)]>=1.0) p6=1.;

            p13 = 0.;
            p14 = 0.;
            p15 = 0.;
            p16 = 0.;
            p35 = 0.;
            p36 = 0.;

            if(p[1][ic(iw+1,jw+1,kw)]>=1.0) p13=1.;
            if(p[1][ic(iw+1,jw-1,kw)]>=1.0) p14=1.;
            if(p[1][ic(iw+1,jw,kw+1)]>=1.0) p15=1.;
            if(p[1][ic(iw+1,jw,kw-1)]>=1.0) p16=1.;
            if(p[1][ic(iw,jw+1,kw+1)]>=1.0) p35=1.;
            if(p[1][ic(iw,jw+1,kw-1)]>=1.0) p36=1.;

            p23 = 0.;
            p24 = 0.;
            p25 = 0.;
            p26 = 0.;
            p45 = 0.;
            p46 = 0.;

            if(p[1][ic(iw-1,jw+1,kw)]>=1.0) p23=1.;
            if(p[1][ic(iw-1,jw-1,kw)]>=1.0) p24=1.;
            if(p[1][ic(iw-1,jw,kw+1)]>=1.0) p25=1.;
            if(p[1][ic(iw-1,jw,kw-1)]>=1.0) p26=1.;
            if(p[1][ic(iw,jw-1,kw+1)]>=1.0) p45=1.;
            if(p[1][ic(iw,jw-1,kw-1)]>=1.0) p46=1.;

            pxp =0.5*(p1+p0);
            pxm =0.5*(p0+p2);
            pyp =0.5*(p3+p0);
            pym =0.5*(p0+p4);
            pzp =0.5*(p5+p0);
            pzm =0.5*(p0+p6);
            pxpyp=0.25*(p0+p1+p3+p13);
            pypzp=0.25*(p0+p3+p5+p35);
            pzpxp=0.25*(p0+p5+p1+p15);
            pxmym=0.25*(p0+p2+p4+p24);
            pymzm=0.25*(p0+p4+p6+p46);
            pzmxm=0.25*(p0+p6+p2+p26);
            pxpym=0.25*(p0+p1+p4+p14);
            pxmyp=0.25*(p0+p2+p3+p23);
            pypzm=0.25*(p0+p3+p6+p36);
            pymzp=0.25*(p0+p4+p5+p45);
            pzpxm=0.25*(p0+p5+p2+p26);
            pzmxp=0.25*(p0+p6+p1+p16);


            dpmag0=sqrt(pow((pxp-pxm)/dxl,2.)+pow((pyp-pym)/dyl,2.)+pow((pzp-pzm)/dzl,2.));
            dpmagxp=sqrt(pow(2.*(pxp-p0)/dxl,2.)+pow((pxpyp-pxpym)/dyl,2.)+pow((pzpxp-pzmxp)/dzl,2.));
            dpmagxm=sqrt(pow(2.*(p0-pxm)/dxl,2.)+pow((pxmyp-pxmym)/dyl,2.)+pow((pzpxm-pzmxm)/dzl,2.));
            dpmagyp=sqrt(pow((pxpyp-pxmyp)/dxl,2.)+pow(2.*(pyp-p0)/dyl,2.)+pow((pypzp-pypzm)/dzl,2.));
            dpmagym=sqrt(pow((pxpym-pxmym)/dxl,2.)+pow(2.*(p0-pym)/dyl,2.)+pow((pymzp-pymzm)/dzl,2.));
            dpmagzp=sqrt(pow((pzpxp-pzpxm)/dxl,2.)+pow((pypzp-pymzp)/dyl,2.)+pow(2.*(pzp-p0)/dzl,2.));
            dpmagzm=sqrt(pow((pzmxp-pzmxm)/dxl,2.)+pow((pypzm-pymzm)/dyl,2.)+pow(2.*(p0-pzm)/dzl,2.));
            nvecxp_x=0.; nvecxp_y=0.; nvecxp_z=0.;
            nvecxm_x=0.; nvecxm_y=0.; nvecxm_z=0.;
            nvecyp_x=0.; nvecyp_y=0.; nvecyp_z=0.;
            nvecym_x=0.; nvecym_y=0.; nvecym_z=0.;
            nveczp_x=0.; nveczp_y=0.; nveczp_z=0.;
            nveczm_x=0.; nveczm_y=0.; nveczm_z=0.;
            //if (fabs(dpmagxp)>1.e-2/dxl) {
            if (fabs(dpmagxp)>0.) {
                nvecxp_x = -2.*(pxp-p0)/dxl/dpmagxp;
                nvecxp_y = -(pxpyp-pxpym)/dyl/dpmagxp;
                nvecxp_z = -(pzpxp-pzmxp)/dzl/dpmagxp;
            }
            //if (fabs(dpmagxm)>1.e-2/dxl) {
            if (fabs(dpmagxm)>0.) {
                nvecxm_x = -2.*(p0-pxm)/dxl/dpmagxm;
                nvecxm_y = -(pxmyp-pxmym)/dyl/dpmagxm;
                nvecxm_z = -(pzpxm-pzmxm)/dzl/dpmagxm;
            }
            //if (fabs(dpmagyp)>1.e-2/dyl) {
            if (fabs(dpmagyp)>0.) {
                nvecyp_x = -(pxpyp-pxmyp)/dxl/dpmagyp;
                nvecyp_y = -2.*(pyp-p0)/dyl/dpmagyp;
                nvecyp_z = -(pypzp-pypzm)/dzl/dpmagyp;
            }
            //if (fabs(dpmagym)>1.e-2/dyl) {
            if (fabs(dpmagym)>0.) {
                nvecym_x = -(pxpym-pxmym)/dxl/dpmagym;
                nvecym_y = -2.*(p0-pym)/dyl/dpmagym;
                nvecym_z = -(pymzp-pymzm)/dzl/dpmagym;
            }
            //if (fabs(dpmagzp)>1.e-2/dzl) {
            if (fabs(dpmagzp)>0.) {
                nveczp_x = -(pzpxp-pzpxm)/dxl/dpmagzp;
                nveczp_y = -(pypzp-pymzp)/dyl/dpmagzp;
                nveczp_z = -2.*(pzp-p0)/dzl/dpmagzp;
            }
            //if (fabs(dpmagzm)>1.e-2/dzl) {
            if (fabs(dpmagzm)>0.) {
                nveczm_x = -(pzmxp-pzmxm)/dxl/dpmagzm;
                nveczm_y = -(pypzm-pymzm)/dyl/dpmagzm;
                nveczm_z = -2.*(p0-pzm)/dzl/dpmagzm;
            }
            nvecmagxp = sqrt(pow(nvecxp_x,2.)+pow(nvecxp_y,2.)+pow(nvecxp_z,2.));
            nvecmagxm = sqrt(pow(nvecxm_x,2.)+pow(nvecxm_y,2.)+pow(nvecxm_z,2.));
            nvecmagyp = sqrt(pow(nvecyp_x,2.)+pow(nvecyp_y,2.)+pow(nvecyp_z,2.));
            nvecmagym = sqrt(pow(nvecym_x,2.)+pow(nvecym_y,2.)+pow(nvecym_z,2.));
            nvecmagzp = sqrt(pow(nveczp_x,2.)+pow(nveczp_y,2.)+pow(nveczp_z,2.));
            nvecmagzm = sqrt(pow(nveczm_x,2.)+pow(nveczm_y,2.)+pow(nveczm_z,2.));
            //if (p0>=0.999) {
            if (1) {
            J0xp_x = -diffu0*2.*(Cxp-C0)/dxl + Cxp*vxxp;
            J0xp_y = -diffu0*(Cxpyp-Cxpym)/dyl + Cxp*vyxp;
            J0xp_z = -diffu0*(Czpxp-Czmxp)/dzl + Cxp*vzxp;
            J0xm_x = -diffu0*2.*(C0-Cxm)/dxl + Cxm*vxxm;
            J0xm_y = -diffu0*(Cxmyp-Cxmym)/dyl + Cxm*vyxm;
            J0xm_z = -diffu0*(Czpxm-Czmxm)/dzl + Cxm*vzxm;
            J0yp_x = -diffu0*(Cxpyp-Cxmyp)/dxl + Cyp*vxyp;
            J0yp_y = -diffu0*2.*(Cyp-C0)/dyl + Cyp*vyyp;
            J0yp_z = -diffu0*(Cypzp-Cypzm)/dzl + Cyp*vzyp;
            J0ym_x = -diffu0*(Cxpym-Cxmym)/dxl + Cym*vxym;
            J0ym_y = -diffu0*2.*(C0-Cym)/dyl + Cym*vyym;
            J0ym_z = -diffu0*(Cymzp-Cymzm)/dzl + Cym*vzym;
            J0zp_x = -diffu0*(Czpxp-Czpxm)/dxl + Czp*vxzp;
            J0zp_y = -diffu0*(Cypzp-Cymzp)/dyl + Czp*vyzp;
            J0zp_z = -diffu0*2.*(Czp-C0)/dzl + Czp*vzzp;
            J0zm_x = -diffu0*(Czmxp-Czmxm)/dxl + Czm*vxzm;
            J0zm_y = -diffu0*(Cypzm-Cymzm)/dyl + Czm*vyzm;
            J0zm_z = -diffu0*2.*(C0-Czm)/dzl + Czm*vzzm;
            }

            J0dotnxp = J0xp_x*nvecxp_x+J0xp_y*nvecxp_y+J0xp_z*nvecxp_z;
            J0dotnxm = J0xm_x*nvecxm_x+J0xm_y*nvecxm_y+J0xm_z*nvecxm_z;
            J0dotnyp = J0yp_x*nvecyp_x+J0yp_y*nvecyp_y+J0yp_z*nvecyp_z;
            J0dotnym = J0ym_x*nvecym_x+J0ym_y*nvecym_y+J0ym_z*nvecym_z;
            J0dotnzp = J0zp_x*nveczp_x+J0zp_y*nveczp_y+J0zp_z*nveczp_z;
            J0dotnzm = J0zm_x*nveczm_x+J0zm_y*nveczm_y+J0zm_z*nveczm_z;
            /* 
            JBxp = -1.*J0dotnxp*nvecxp_x;
            JBxm = -1.*J0dotnxm*nvecxm_x;
            JByp = -1.*J0dotnyp*nvecyp_y;
            JBym = -1.*J0dotnym*nvecym_y;
            JBzp = -1.*J0dotnzp*nveczp_z;
            JBzm = -1.*J0dotnzm*nveczm_z;
           */
            
            JBxp =  - J0xp_x*nvecmagxp;
            JBxm =  - J0xm_x*nvecmagxm;
            JByp =  - J0yp_y*nvecmagyp;
            JBym =  - J0ym_y*nvecmagym;
            JBzp =  - J0zp_z*nvecmagzp;
            JBzm =  - J0zm_z*nvecmagzm;
            
            divJ0 =  (J0xp_x-J0xm_x)/(dxl);
            divJ0 +=  (J0yp_y-J0ym_y)/(dyl);
            divJ0 +=  (J0zp_z-J0zm_z)/(dzl);
            divJB =  (JBxp-JBxm)/(dxl);
            divJB +=  (JByp-JBym)/(dyl);
            divJB +=  (JBzp-JBzm)/(dzl);
            divJ = divJ0 + divJB;
            //divJ = divJ0;
            if(1)Cnew[ic(iw,jw,kw)] = C[ic(iw,jw,kw)] + dtime*(-divJ);
            //if (Cnew[ic(iw,jw,kw)]<0.) Cnew[ic(iw,jw,kw)]=0.; 
            Ctot+=Cnew[ic(iw,jw,kw)];
	}
	}
	}

        return Ctot;
}

double calc_Ce2(double *Cnew, double *C, double *p[nphase], double *V)
{
        //////////////////
        double calc_diffusivity(double p);
        //////////////////
        extern int nx,ny,ny,ntot;
        extern double dxl,dyl,dzl;
        extern double THERM;
        extern double z[ncomp];
        extern double DIFFUSIVITY;
        extern double dtime;
	//////////////////
	int i,j,k;
        int iw,jw,kw;

        //Concentration
        double Cc;
        double Cip, Cjp, Ckp;
        double Cim, Cjm, Ckm;

        //Electric field
        double EXxp,EYyp,EZzp;
        double EXxm,EYym,EZzm;
        //Misc. values for electro-convection term
        double mole_Coulomb = z[0]*FARAD;

        //Phase field
        double pc;
        double pip,pjp,pkp;
        double pim,pjm,pkm;

        //Fluxes
        double Jxp, Jyp, Jzp;
        double Jxm, Jym, Jzm;
        double JDxp,JDyp,JDzp;
        double JDxm,JDym,JDzm;
        double JExp,JEyp,JEzp;
        double JExm,JEym,JEzm;

        //Divergences of fluxes
	double divJ,divJD,divJE;

        //Diffusivity
        double Dc;
        double Dip,Djp,Dkp;
        double Dim,Djm,Dkm;
       
        //Interpolated diffusivity
        double Dxp,Dyp,Dzp;
        double Dxm,Dym,Dzm;

        //Coefficients for electro-convection term
        double Ac;
        double Aip,Ajp,Akp;
        double Aim,Ajm,Akm;

        //Interpolated coefficients for electro-convection term
        double Axp,Ayp,Azp;
        double Axm,Aym,Azm;

        double Ctot=0.;
	//////////////////

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            iw=ii(i);
            jw=ii(j);
            kw=ii(k);

            Cc = C[ic(iw,jw,kw)];
            Cip = C[ic(iw+1,jw,kw)];
            Cim = C[ic(iw-1,jw,kw)];
            Cjp = C[ic(iw,jw+1,kw)];
            Cjm = C[ic(iw,jw-1,kw)];
            Ckp = C[ic(iw,jw,kw+1)];
            Ckm = C[ic(iw,jw,kw-1)];

            EXxp = -(V[ic(iw+1,jw,kw)]-V[ic(iw,jw,kw)])/dxl;
            EXxm = -(V[ic(iw,jw,kw)]-V[ic(iw-1,jw,kw)])/dxl;
            EYyp = -(V[ic(iw,jw+1,kw)]-V[ic(iw,jw,kw)])/dyl;
            EYym = -(V[ic(iw,jw,kw)]-V[ic(iw,jw-1,kw)])/dyl;
            EZzp = -(V[ic(iw,jw,kw+1)]-V[ic(iw,jw,kw)])/dzl;
            EZzm = -(V[ic(iw,jw,kw)]-V[ic(iw,jw,kw-1)])/dzl;
    
            pc = p[1][ic(iw,jw,kw)];
            pip =p[1][ic(iw+1,jw,kw)];
            pim =p[1][ic(iw-1,jw,kw)];
            pjp =p[1][ic(iw,jw+1,kw)];
            pjm =p[1][ic(iw,jw-1,kw)];
            pkp =p[1][ic(iw,jw,kw+1)];
            pkm =p[1][ic(iw,jw,kw-1)];

            Dc = calc_diffusivity(pc);
            Dip = calc_diffusivity(pip);
            Dim = calc_diffusivity(pim);
            Djp = calc_diffusivity(pjp);
            Djm = calc_diffusivity(pjm);
            Dkp = calc_diffusivity(pkp);
            Dkm = calc_diffusivity(pkm);

            Dxp = 2.*1./(1./Dip + 1./Dc);
            Dxm = 2.*1./(1./Dim + 1./Dc);
            Dyp = 2.*1./(1./Djp + 1./Dc);
            Dym = 2.*1./(1./Djm + 1./Dc);
            Dzp = 2.*1./(1./Dkp + 1./Dc);
            Dzm = 2.*1./(1./Dkm + 1./Dc);

            Ac = Dc*Cc*mole_Coulomb/THERM;
            Aip = Dip*Cip*mole_Coulomb/THERM;
            Aim = Dim*Cim*mole_Coulomb/THERM;
            Ajp = Djp*Cjp*mole_Coulomb/THERM;
            Ajm = Djm*Cjm*mole_Coulomb/THERM;
            Akp = Dkp*Ckp*mole_Coulomb/THERM;
            Akm = Dkm*Ckm*mole_Coulomb/THERM;

            Axp = 2.*1./(1./Aip + 1./Ac);
            Axm = 2.*1./(1./Aim + 1./Ac);
            Ayp = 2.*1./(1./Ajp + 1./Ac);
            Aym = 2.*1./(1./Ajm + 1./Ac);
            Azp = 2.*1./(1./Akp + 1./Ac);
            Azm = 2.*1./(1./Akm + 1./Ac);

            JDxp = -Dxp*(Cip-Cc)/dxl;
            JDxm = -Dxm*(Cc-Cim)/dxl;
            JDyp = -Dyp*(Cjp-Cc)/dyl;
            JDym = -Dym*(Cc-Cjm)/dyl;
            JDzp = -Dzp*(Ckp-Cc)/dzl;
            JDzm = -Dzm*(Cc-Ckm)/dzl;

            JExp = Axp*EXxp;
            JExm = Axm*EXxm;
            JEyp = Ayp*EYyp;
            JEym = Aym*EYym;
            JEzp = Azp*EZzp;
            JEzm = Azm*EZzm;
           
            divJD =  (JDxp-JDxm)/(dxl);
            divJD +=  (JDyp-JDym)/(dyl);
            divJD +=  (JDzp-JDzm)/(dzl);
            divJE =  (JExp-JExm)/(dxl);
            divJE +=  (JEyp-JEym)/(dyl);
            divJE +=  (JEzp-JEzm)/(dzl);
            divJ = divJD + divJE;
            //divJ = divJD;
            if(1)Cnew[ic(iw,jw,kw)] = C[ic(iw,jw,kw)] + dtime*(-divJ);
            //if (Cnew[ic(iw,jw,kw)]<0.) Cnew[ic(iw,jw,kw)]=0.; 
            Ctot+=Cnew[ic(iw,jw,kw)];
	}
	}
	}

        return Ctot;
}

double calc_diffusivity(double p)
{
    ///////////////
    extern double DIFFUSIVITY;
    ///////////////
    double result=0.;
    ///////////////
    if(p > 0.0){
        result = DIFFUSIVITY;
    }
    //result = DIFFUSIVITY;
    return result;
}

void solvPF_CA(double *pnew[nphase], double *pold[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dG(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
    void calc_EMM(double EM[nphase][nphase], double *p[nphase],int i,int j,int k);
    void trim_p(double *p[nphase],int iw,int jw,int kw);
	//////////////////
	int i,j,k;

    double EM[nphase][nphase];
	double dG[nphase][nphase];
	double delta[nphase]={0.};
    double p_l,p_m;
	int lp,mp;
    int npw = nphase;
	//////////////////

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
		calc_EMM(EM,pold,i,j,k);
		calc_dG(dG,pold,i,j,k);

		for (lp=0;lp<nphase;lp++) {
                    p_l = pold[lp][ic(ii(i),ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)+1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i)-1,ii(j),ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)+1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j)-1,ii(k))];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)+1)];
                    p_l += pold[lp][ic(ii(i),ii(j),ii(k)-1)];
		    delta[lp]=0.;
                    if (p_l>0.+0.01 && p_l<7.-0.01) {
			for (mp=0;mp<nphase;mp++) {
                            p_m = pold[mp][ic(ii(i),ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)+1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i)-1,ii(j),ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)+1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j)-1,ii(k))];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)+1)];
                            p_m += pold[mp][ic(ii(i),ii(j),ii(k)-1)];
                            if (p_m>0.+0.01 && p_m<7.-0.01) {
				if (lp!=mp) {
					delta[lp] += -EM[lp][mp]*dG[lp][mp];
				}
                            }
			}
                    }
                    //RR[lp]+=delta[lp]/(-EMM);
		    pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dtime;
                    //only for parabolic double well potential
		}
                trim_p(pnew,ii(i),ii(j),ii(k));

	}
	}
	}
}
void solvPF_CH(double *pnew[nphase], double *pold[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dG(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
    void calc_EMM(double EM[nphase][nphase], double *p[nphase],int i,int j,int k);
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
			pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dtime;
		}

	}
	}
	}

}
void trim_p(double *p[nphase],int iw,int jw,int kw)
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
        if(p[np][ic(iw,jw,kw)]<=0.)p[np][ic(iw,jw,kw)]=0.;
        if(p[np][ic(iw,jw,kw)]>=1.)p[np][ic(iw,jw,kw)]=1.;
        psum+=p[np][ic(iw,jw,kw)];
    }
    if(pfix>=0 && pfix<nphase) { 
        psum -= p[pfix][ic(iw,jw,kw)];
    }
    for (np=0;np<nphase;np++) {
        if (pfix>=0 && pfix<nphase) {
            if(np!=pfix && psum!=0.)p[np][ic(iw,jw,kw)]*=(1.-p[pfix][ic(iw,jw,kw)])/psum;
        }else{
            p[np][ic(iw,jw,kw)]/=psum;
        }
    }
    
}
void trim_p_rlx(double *p[nphase],int iw,int jw,int kw)
{
    /////////
    double psum=0.;
    int np;
    int npw=nphase; 
    int plist[nphase];
    /////////
    npw=makephaselist(plist,p,iw,jw,kw);
    for (np=0;np<nphase;np++) {
        if(p[np][ic(iw,jw,kw)]<=0.)p[np][ic(iw,jw,kw)]=0.;
        if(p[np][ic(iw,jw,kw)]>=1.)p[np][ic(iw,jw,kw)]=1.;
        psum+=p[np][ic(iw,jw,kw)];
    }
    for (np=0;np<nphase;np++) {
        p[np][ic(iw,jw,kw)]/=psum;
    }
    
}
void update_p(double *pold[nphase], double *pnew[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	////////////////
	int i,j,k,np;
        int iw,jw,kw;
	////////////////
	for (np=0;np<nphase;np++) {
	    for (i=1;i<=nx;i++) {
	    for (j=1;j<=ny;j++) {
	    for (k=1;k<=nz;k++) {
                iw=ii(i);
                jw=ii(j);
                kw=ii(k);
	        pold[np][ic(iw,jw,kw)]=pnew[np][ic(iw,jw,kw)];
	    }
	    }
	    }
	}
}

//double calc_v_fromEpot(double *v[3], double *Epot, double *p[nphase])
double calc_v_fromEpot(double *v[3], double *Epot)
{
    //////////////////
    extern int nx,ny,ny,ntot;
    extern double DIFFUSIVITY;
    extern double THERM;
    extern double z[ncomp];
    ////////////////
    int i,j,k,np;
    int iw,jw,kw;
    double mole_Coulomb=z[0]*FARAD;
    double velmax=0.,speed;
    double diffu;
    ////////////////
    for (i=1;i<=nx;i++) {
    for (j=1;j<=ny;j++) {
    for (k=1;k<=nz;k++) {
        iw=ii(i);
        jw=ii(j);
        kw=ii(k);
        //diffu = calc_diffusivity(p[1][ic(iw,jw,kw)]);
        diffu = DIFFUSIVITY;
        v[0][ic(iw,jw,kw)] = -(Epot[ic(iw+1,jw,kw)]-Epot[ic(iw-1,jw,kw)])/(2.*dxl)*diffu*mole_Coulomb/THERM;
        v[1][ic(iw,jw,kw)]  = -(Epot[ic(iw,jw+1,kw)]-Epot[ic(iw,jw-1,kw)])/(2.*dyl)*diffu*mole_Coulomb/THERM;
        v[2][ic(iw,jw,kw)]  = -(Epot[ic(iw,jw,kw+1)]-Epot[ic(iw,jw,kw-1)])/(2.*dzl)*diffu*mole_Coulomb/THERM;
        speed = sqrt(v[0][ic(iw,jw,kw)] *v[0][ic(iw,jw,kw)] +v[1][ic(iw,jw,kw)] *v[1][ic(iw,jw,kw)] +v[2][ic(iw,jw,kw)] *v[2][ic(iw,jw,kw)] );
        if (speed >= velmax) velmax = speed;
    }
    }
    }
    return velmax;
}



double find_timestep(int PFkind,int DIFFUkind)
{
    /*
    PFkind: What kind of equation solving for phase field is it?
        1: Cahn-Allen type
        2: Cahn-Hilliard type
    DIFFUkind: What kind of equation solving for diffusion equation is it?
        1: Diffusion only
        2: Convection-diffusion equation
    */ 
    ////////////////////
    extern double vmax;
    extern double DIFFUSIVITY;
    extern double maxPeC;
    extern double dxl,dyl,dzl;
    extern double EPSMAX, EMM;
    extern double dt, dtime;
    ////////////////////
    double mindtfind(double,double,double,int,int,int);
    ////////////////////
    double dtkind1, dtkind2, dtkind3;
    double dtkind1CA,dtkind1CH;
    double dtkind2DIFF,dtkind2CONVDIFF;
    double dtkindCOURANT;
    double dtcurr,dtmin=9999999.;
    ////////////////////
    maxPeC = fabs(vmax*dxl/DIFFUSIVITY);
    dtkind1CH = pow(dxl,4.)/(12.*EPSMAX*EPSMAX*EMM)/3.;
    dtkind1CA = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)/3.;
    dtkind2DIFF = 0.5*pow(dxl,2.)/DIFFUSIVITY/3.;
    if (maxPeC!=0.) {
      //dtkind2CONVDIFF = pow(dxl,2.)/DIFFUSIVITY*(1./fabs(2.*maxPeC))/3.;
      dtkind2CONVDIFF = pow(dxl,2.)/DIFFUSIVITY*(1./fabs(1.+maxPeC))/3.;
    }else{
      dtkind2CONVDIFF = 99999.;
    }
    if (vmax!=0.) {
        dtkindCOURANT = dxl/vmax/3.;
    }else{
        dtkindCOURANT = 999999.;
    }
    dtkind1 = dtkind1CA;
    if (PFkind==2) dtkind1 = dtkind1CH;
    dtkind2 = dtkind2DIFF;
    if (DIFFUkind==2) dtkind2 = dtkind2CONVDIFF;
    dtkind3 = dtkindCOURANT;
    //dtcurr = mindtfind(dtkind1,dtkind2,dtkind3,1,2,3);
    dtcurr = mindtfind(9999.,dtkind2,dtkind3,1,2,3);
    if(dtcurr <= dtmin) dtmin = dtcurr;
    set_dtime(dtmin*tstable);

    return maxPeC;
}


double mindtfind(double a, double b, double c,int kind1,int kind2,int kind3)
{
    ///////////////
    double result=a;
    ///////////////
    set_dtkind(kind1);
    if (result >  b) {
        result = b;
        set_dtkind(kind2);
    }
    if (result >  c) {
        result = c;
        set_dtkind(kind3);
    }
    return result;
}

void timestepwrite(void)
{
    ///////////////////
    extern double maxPeC;
    ///////////////////
    FILE *outdt;
    ///////////////////
    outdt = fopen("timesteps.txt","a");
    fprintf(outdt,"%d\t%le\t%le\t%d\t%le\n",itime, time, dtime, dtkind, maxPeC);
    fclose(outdt);
}

void set_dtime(double val)
{
    /////////////////////
    extern double dtime;
    /////////////////////
    dtime = val;
}
void set_dtkind(int val)
{
    dtkind = val;
}
