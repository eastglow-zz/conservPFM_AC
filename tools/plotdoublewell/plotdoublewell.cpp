#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dbl(double phi1, double phi2);
double triple(double phi1, double phi2, double phi3);
int main()
{
  ///////////
  double x,y,dx,dy;
  double phi1, phi2, phi3;
  double dphi1, dphi2, dphi3;
  double f;
  double alpha;
  FILE *fout;
  char fname[50];
  int numphi1,numphi2;
  ///////////
  dphi1 = 0.01;
  dphi2 = 0.01;
  printf("alpha?:");
  scanf("%lf",&alpha);
  sprintf(fname,"f_%lf.plt",alpha);
  fout = fopen(fname,"w");
  fprintf(fout,"zone, i=, j=, f=point\n");
  for (phi1=0.,numphi1=1;phi1<1.+dphi1;phi1+=dphi1,numphi1++) {
    for (phi2=1.-phi1,numphi2=1;phi1+phi2>0.-dphi2 && phi2>0.;phi2-=dphi2,numphi2++) {
       phi3 = 1.-phi1-phi2;
       if (phi3>=0. && phi3 <=1.) {
         x = phi1;
         y = phi2;
         f = dbl(phi1,phi2) + alpha*triple(phi1,phi2,phi3);
         fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,phi1,phi2,phi3,f);
       }
    }
  }
  printf("numphi1 = %d, numphi2 = %d\n",numphi1,numphi2);

  fclose(fout);
  return 0;
}

double dbl(double phi1, double phi2)
{
  ////////////
  double result;
  ////////////
  result = phi1*phi1*phi2*phi2;
  return result;
}
double triple(double phi1, double phi2, double phi3)
{
  /////////////
  double result;
  /////////////
  result =pow(phi1*phi2*phi3,2.);
  return result;
}
