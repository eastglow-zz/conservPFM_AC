#include <stdlib.h>
#include <math.h>
#include "defined_values.h"
#include "common_functions.h"

extern int nx,ny,nz,ntot,nmax;

int pfix =  2;

double dxl = 0.5e-6;
double dyl = dxl;
double dzl = dxl;
double XI = 3.*dxl; // half interface thickness

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
double time = 0.;
int itime = 0;

double RR[nphase]={0.};
double V0[nphase]={0.};
double HH[nphase]={0.};
double lamb[nphase]; //interphase driving force coefficient

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
    extern int nx,ny,nz;
    //return (i)*(ny+4)*(nz+4)+(j)*(nz+4)+(k); //x-outermost
    return (i)+(nx+4)*(j)+(nx+4)*(ny+4)*k; //z-outermost
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


void calc_dDRIV(double dG[nphase][nphase], double *p[nphase],int i,int j,int k)
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

void drivPF_CA(double *pnew[nphase], double *pold[nphase])
{
    //////////////////
    extern int nx,ny,ny,ntot;
	//////////////////
	void calc_dDRIV(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
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

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

	for (i=1;i<=nx;i++) {
	for (j=1;j<=ny;j++) {
	for (k=1;k<=nz;k++) {
            //pold[2][ii(i)][ii(j)][ii(k)]<=0 ? npw=nphase-1 : npw=nphase;

		calc_EMM(EM,pold,i,j,k);
		calc_dDRIV(dG,pold,i,j,k);

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
		    pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dt;
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
	//////////////////

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

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
	            pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dt;
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
			pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dt;
		}

	}
	}
	}

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

	dt = pow(dxl,2.)/(2.*EPSMAX*EPSMAX*EMM)*tstable/3.;

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
		    pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dt;
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
			pnew[lp][ic(ii(i),ii(j),ii(k))]=pold[lp][ic(ii(i),ii(j),ii(k))]+delta[lp]*dt;
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


