/**********************************************************************
 *                                                                    *
 * HILL_STAB.C                                                        *
 *                                                                    *
 * Author: Rory Barnes (rory@astro.washington.edu)                    *
 *                                                                    *
 * To compile: gcc -o hillstab hill_stab.c -lm                        *
 *                                                                    *
 * This code calculates the relative proximity of a three-body system *
 * to "Hill stability." A Hill stable system is one for which the     *
 * the ordering of the planets remains constant, i.e. the most        *
 * distent body may escape to infinity and the system would still be  *
 * Hill stable. Hill stability may be calculated for any system, but  *
 * this code is optimized for planetary systems. This code was used   *
 * in R. Barnes & R. Greenberg, 2006, ApJ, 674, L163-L166 to          *
 * demonstrate that for a system of a solar-type star and two         *
 * Jupiter-ish mass planets that Hill stability approximates          *
 * "Lagrange stability," which means no swapping && no ejections.     *
 *                                                                    *
 * The user inputs the central mass in solar units and the orbiters   *
 * parameters are entered into the next two lines with the format:    *
 * Mass SemiMajorAxis Eccentricity Inclination ArgumentPericenter     *
 * LongAscNode MeanAnomaly. The units are Jupiter masses, AU,         *
 * and degrees. The final line must state either "bodycentric" or     *
 * "barycentric" to indicate the coordinate system of the orbital     *
 * elements. There are no command line options.                       *
 *                                                                    *
 * This code will output 2 numbers: Exact and Approx. Both numbers    *
 * represent relative proximity to the boundary, with unity on the    *
 * boundary, values < 1 are Hill unstable, and > 1 are Hill stable.   * 
 * Exact is the value of beta/beta_{crit} from BG06, which is         *
 * computed by calculating the energy and angular momentum including  *
 * the primary. Approx is the value of delta/delta_{crit} from BG06,  *
 * which is computed assuming the central mass dominates, see Gladman *
 * (1993). As this option assumes the central body's center is very   *
 * close to the system's center-of-mass, Approx uses the input        *
 * elements.                                                          *
 *                                                                    *
 * Thanks to Russell Deitrick for catching a bug in read_init()!      *
 *                                                                    *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#define dot(a,b)        (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#define BIGG	6.672e-8
#define AUCM	1.49598e13
#define MSUN	1.98892e33
#define MJUP    1.8987e30

typedef struct {
  double a,e,i,lasc,aper,mean_an;
} ELEMS;

char *sLower(char cString[]) {
  int iPos;
  for (iPos=0;cString[iPos];iPos++)
    cString[iPos] = tolower(cString[iPos]);
    
  return cString;
}

void hel_bar(double **hel,double **bar,double *ms,double *m,int n) {
   int i,p;

   for(i=1;i<=3;i++)
      bar[0][i] = 0;
   for(p=1;p<=n;p++) {
      for(i=1;i<=3;i++)
         bar[0][i] -= m[p]/ms[n]*hel[p][i];
   }
   for(p=1;p<=n;p++) {
      for(i=1;i<=3;i++)
         bar[p][i] = hel[p][i]+bar[0][i];
   }
}

void read_init(char *infile, ELEMS *p1, ELEMS *p2,double *m,int *origin) {
  char c_origin[24];
  FILE *fp;

  fp=fopen(infile,"r");

  fscanf(fp,"%lf",&m[0]);
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&m[1],&p1->a,&p1->e,&p1->i,&p1->aper,&p1->lasc,&p1->mean_an);
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&m[2],&p2->a,&p2->e,&p2->i,&p2->aper,&p2->lasc,&p2->mean_an);
  fscanf(fp,"%s",c_origin);
  if (memcmp(sLower(c_origin),"ba",2) == 0) 
    *origin = 0;
  else if (memcmp(sLower(c_origin),"bo",2) == 0) 
    *origin = 1;
  else {
    fprintf(stderr,"ERROR: Unknown coordinate system. Options are barycentric or bodycentric.\n");
    exit(1);
  }

  /* Convert to cgs, rads */
  m[0] *= MSUN;
  m[1] *= MJUP;
  m[2] *= MJUP;
  p1->a *= AUCM;
  p1->aper *= M_PI/180;
  p1->i *= M_PI/180;
  p1->lasc *= M_PI/180;
  p2->a *= AUCM;
  p2->aper *= M_PI/180;
  p2->i *= M_PI/180;
  p2->lasc *= M_PI/180;
}

void cartes(double *x, double *v, ELEMS elem_ptr,double mu) {
  double a,e,m,cosi,sini,cos_lasc,sin_lasc,cos_aper,sin_aper;
  double es,ec,w,wp,wpp,wppp,ecc,dx,lo,up,next;int iter;
  double sin_ecc,cos_ecc,l1,m1,n1,l2,m2,n2;
  double xi,eta,vel_scl;

  a = elem_ptr.a;
  e = elem_ptr.e;
  m = elem_ptr.mean_an;
  cosi = cos(elem_ptr.i);
  sini = sin(elem_ptr.i);
  cos_lasc = cos(elem_ptr.lasc);
  sin_lasc = sin(elem_ptr.lasc);
  cos_aper = cos(elem_ptr.aper);
  sin_aper = sin(elem_ptr.aper);

  /*
   * Reduce mean anomoly to [0, 2*PI)
   */
  m -= ((int)(m/(2*M_PI)))*2*M_PI;
  /*
   * Solve kepler's equation.
   */
  if (sin(m)>0)
    ecc = m+0.85*e;
  else 
    ecc = m-0.85*e;
  lo = -2*M_PI;
  up = 2*M_PI;
  for(iter=1;iter<=32;iter++) {
    es = e*sin(ecc);
    ec = e*cos(ecc);
    w = ecc-es-m;
    wp = 1-ec;
    wpp = es;
    wppp = ec;
    if (w>0)
      up = ecc;
    else 
      lo = ecc;
    dx = -w/wp;
    dx = -w/(wp+dx*wpp/2);
    dx = -w/(wp+dx*wpp/2+dx*dx*wppp/6);
    next = ecc+dx;
    if (ecc==next) 
      break;
    if ((next>lo) && (next<up))
      ecc= next;
    else ecc= (lo+up)/2;
    if((ecc==lo) || (ecc==up))
      break;
    if (iter>30)
      printf("%4d %23.20f %e\n",iter,ecc,up-lo);
  }
  if(iter>32) {
    fprintf(stderr,"ERROR: Kepler solultion failed.\n");
    exit(1);
  }

  cos_ecc = cos(ecc);
  sin_ecc = sin(ecc);

  l1 = cos_lasc*cos_aper-sin_lasc*sin_aper*cosi;
  m1 = sin_lasc*cos_aper+cos_lasc*sin_aper*cosi;
  n1 = sin_aper*sini;
  l2 = -cos_lasc*sin_aper-sin_lasc*cos_aper*cosi;
  m2 = -sin_lasc*sin_aper+cos_lasc*cos_aper*cosi;
  n2 = cos_aper*sini;

  xi = a*(cos_ecc-e);
  eta = a*sqrt(1-e*e)*sin_ecc;
  x[0] = l1*xi+l2*eta;
  x[1] = m1*xi+m2*eta;
  x[2] = n1*xi+n2*eta;
  vel_scl = sqrt((mu*a)/dot(x,x));
  xi = -vel_scl*sin_ecc;
  eta = vel_scl*sqrt(1-e*e)*cos_ecc;
  v[0] = l1*xi+l2*eta;
  v[1] = m1*xi+m2*eta;
  v[2] = n1*xi+n2*eta;
}

double GetExact(double **r,double **v,double *m) {
  double c,h,ke;
  double br[3][3],bv[3][3];
  double c1[3][3];
  double ctot[3];
  double p_a,crit;
  double mm,mtot;
  double r01,r02,r12;
  int i,p;
 
  for (i=0;i<3;i++) {
    for (p=0;p<3;p++) {
      br[p][i]=0;
      bv[p][i]=0;
    }
  }

  mtot=m[0]+m[1]+m[2];
  mm=m[0]*m[1] + m[0]*m[2] + m[1]*m[2];
  r01=sqrt(pow((r[0][0]-r[1][0]),2) + pow((r[0][1]-r[1][1]),2) +
        pow((r[0][2]-r[1][2]),2));
  r12=sqrt(pow((r[1][0]-r[2][0]),2) + pow((r[1][1]-r[2][1]),2) +
        pow((r[1][2]-r[2][2]),2));
  r02=sqrt(pow((r[0][0]-r[2][0]),2) + pow((r[0][1]-r[2][1]),2) +
        pow((r[0][2]-r[2][2]),2));

  /* Convert to barycentric coordinates */
  for (p=1;p<=2;p++) {
    for (i=0;i<3;i++) {
       br[0][i] -= m[p]/mtot*r[p][i];
       bv[0][i] -= m[p]/mtot*v[p][i];
    }
  }
  for (p=1;p<=2;p++) {
    for (i=0;i<3;i++) {
       br[p][i] = r[p][i]+br[0][i];
       bv[p][i] = v[p][i]+bv[0][i];
    }
  }
  /* Total energy and angular momentum */
  ke=0;
  for (p=0;p<3;p++) {
    c1[p][0]=m[p]*(br[p][1]*bv[p][2]-br[p][2]*bv[p][1]);
    c1[p][1]=m[p]*(br[p][2]*bv[p][0]-br[p][0]*bv[p][2]);
    c1[p][2]=m[p]*(br[p][0]*bv[p][1]-br[p][1]*bv[p][0]);
    for (i=0;i<3;i++) 
      ke += 0.5*m[p]*v[p][i]*v[p][i];
  }

  for (i=0;i<3;i++) {
      ctot[i]=0;
      for (p=0;p<3;p++)
	  ctot[i] += c1[p][i];
  }

  h=ke-BIGG*(m[0]*m[1]/r01 + m[0]*m[2]/r02 + m[1]*m[2]/r12);
  c=sqrt(ctot[0]*ctot[0] + ctot[1]*ctot[1] + ctot[2]*ctot[2]);

  p_a = -2*c*c*h*mtot/(BIGG*BIGG*pow(mm,3));

  /* Note that M+B and G93 use m3 as the central mass, but we use m[0] */
  crit = 1 + pow(3,(4.0/3))*(m[1]*m[2])/
       (pow(m[0],(2.0/3))*(pow((m[1]+m[2]),(4.0/3))))
       - (m[1]*m[2]*(11*m[1]+7*m[2]))/(3*m[0]*(m[1]+m[2])*(m[1]+m[2]));

  return p_a/crit;
}

double GetApprox(ELEMS p1,ELEMS p2,double *m) {
  double mu[2],zeta,gamma[2],lambda,dTotMass;
  double p_a,crit;

  dTotMass = m[0] + m[1] + m[2];
  gamma[0] = sqrt(1 - p1.e*p1.e);
  gamma[1] = sqrt(1 - p2.e*p2.e);
  lambda = sqrt(p2.a/p1.a);
  mu[0] = m[1]/m[0];
  mu[1] = m[2]/m[0];
  zeta = mu[0]+mu[1];

  p_a = pow(zeta,-3)*(mu[0]+mu[1]/(lambda*lambda))*pow((mu[0]*gamma[0] + mu[1]*gamma[1]*lambda),2);
  crit = 1 + pow(3,(4./3))*mu[0]*mu[1]/pow(zeta,4./3);
 
  return p_a/crit;
}
 
int main(int argc, char *argv[]) {
  double **x,**v,*m;
  double **bax,**bav;
  double *ms;
  double ratio,e1,e2;
  ELEMS p1,p2;
  double exact,approx;
  int i,j,origin;

  if (argc != 2) {
    fprintf(stderr,"Usage: %s inputfile\n",argv[0]);
    exit(1);
  }

  m=malloc(3*sizeof(double));
  x=malloc(3*sizeof(double*));
  v=malloc(3*sizeof(double*));
  bax=malloc(3*sizeof(double*));
  bav=malloc(3*sizeof(double*));
  ms=malloc(3*sizeof(double));
  for (i=0;i<3;i++) {
    x[i]=malloc(3*sizeof(double));
    v[i]=malloc(3*sizeof(double));
    bax[i]=malloc(3*sizeof(double));
    bav[i]=malloc(3*sizeof(double));
  }

  read_init(argv[1],&p1,&p2,m,&origin);
  /* bodycentric coordinates */
  for (i=0;i<3;i++) {
    x[0][i]=0.0;
    v[0][i]=0.0;
  }

  cartes(x[1],v[1],p1,m[0]*BIGG);
  cartes(x[2],v[2],p2,m[0]*BIGG);
  if (origin) {
    ms[0]=m[0];
    for (i=1;i<3;i++)
      ms[i] = ms[i-1] + m[i];
    hel_bar(x,bax,ms,m,2);
  }
  exact=GetExact(x,v,m);
  approx=GetApprox(p1,p2,m);

  if (exact < 1e-4 || exact > 1e4) 
    printf("Exact: %.5e\n",exact);
  else
    printf("Exact: %.5lf\n",exact);

  if (approx < 1e-4 || approx > 1e4) 
    printf("Approx: %.5e\n",approx);
  else
    printf("Approx: %.5lf\n",approx);

  return 0;
}
