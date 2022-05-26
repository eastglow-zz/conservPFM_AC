double df(int np, double *p1[nphase],int iw,int jw,int kw);
double dg(double phi1, double phi2);
double dg_rlx(double phi1, double phi2);
double dh(double phi);
double dhmn(double phi1,double phi2);
double dmax3(double a1, double a2, double a3);
double dtg(double phi1, double phi2, double phi3);
double f(double p1, double p2, double p3);
double g(double phi1,double phi2);
double g_rlx(double phi1,double phi2);
double h(double phi);
double tg(double phi1, double phi2, double phi3);

double EMMp(int np1, int np2);
double EMMrlx(int np1, int np2);
double EPS(int np1, int np2);
double SIG(int np1, int np2);
double W(int np1, int np2);

int ic(int i, int j, int k);
int makephaselist(int plist[nphase],double *p[nphase],int iw, int jw, int kw);
int max3(int a1, int a2, int a3);

void calc_dDRIV(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
void calc_dG(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
void calc_dG_rlx(double dG[nphase][nphase], double *p[nphase],int i,int j,int k);
void calc_EMM(double EM[nphase][nphase], double *p[nphase],int i,int j,int k);
void del_buff_d(double **p);
void del_buff_i(int **p);
void drivPF_CA(double *pnew[nphase], double *pold[nphase]);
void relaxPF_CA(double *pnew[nphase], double *pold[nphase]);
void relaxPF_CH(double *pnew[nphase], double *pold[nphase]);
void solvPF_CA(double *pnew[nphase], double *pold[nphase]);
void solvPF_CH(double *pnew[nphase], double *pold[nphase]);
void trim_p(double *p[nphase],int iw,int jw,int kw);
void trim_p_rlx(double *p[nphase],int iw,int jw,int kw);
