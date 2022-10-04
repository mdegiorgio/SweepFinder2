#include "my_rand.h"

#define PI 3.141592654



double gammLn(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947,-86.50532032942,
			24.01409824083,-1.231739572450,
			0.1208650973866e-2,-0.5395239385e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -=(x+0.5)*log(tmp);
  ser=1.00000000019015;
  for (j=0; j<5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}



int static z_rndu=137, x_rndu=11, y_rndu=23;
static int tempz=-1, tempx, tempy;

void SetSeed(int seed)
{
  z_rndu = 170*(seed%178) + 137;
  x_rndu = 11;
  y_rndu=23;
}


void saveseed() {
  tempz = z_rndu;
  tempx = x_rndu;
  tempy = y_rndu;
}


void loadseed() {
  if (tempz==-1) 
    fprintf(stderr, "no seed saved; cannot load\n");
  else {
   z_rndu = tempz;
   y_rndu = tempy;
   x_rndu = tempx;
  }
}

double uniform()
{
  /*
    
    U(0,1): AS 183: Appl. Stat. 31:188-190
    Wichmann BA & Hill ID.  1982.  An efficient and portable
    pseudo-random number generator.  Appl. Stat. 31:188-190
    
    x, y, z are any numbers in the range 1-30000.  Integer operation up
    to 30323 required.
    
    Suggested to me by Z. Yang who also provided me with
    the source code used here.
  */
  
  double r;
  
  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177); 
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176); 
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178); 
  if (x_rndu<0) x_rndu+=30269; 
  if (y_rndu<0) y_rndu+=30307; 
  if (z_rndu<0) z_rndu+=30323; 
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0; 
  return (r - (int)r); 
}



int rand_int(int N)
{
  return uniform()*N;
}


int Binomial(int n, double p) {
  int i, k=0;
  for (i=0; i<n; i++) 
    if (uniform() < p) k++;
  return k;
}


//not tested or used ever
int Geometric(double p) {
  int i=0;
  double u;
  do {
    u=uniform();
    i++;
  }
  while (u < p);
  return i;
}


/* Generates an exponential random variable*/
double expo(double c)
{
  return - (1.0/c)*log(uniform());
}


/*generating a gamma when a<1*/
double gamsmall(double a)
{
  double c, d, e, z, x;
  c=1.0/a;
  d=pow(a, a/(1.0-a))*(1.0-a);
  do {
    z=expo(1.0);
    e=expo(1.0);
    x=pow(z,c);
  }
  while (z+e<=d+x);
  return x;
}


/*generating a gamma when a>1*/
double gambig(double a)
{
  double b, c, w, y, x, u, v, z, acc=0;
  b = a-1.0;
  c = 3.0*a-0.75;
  do {
    u = uniform();
    v = uniform();
    w = u*(1-u);
    y = sqrt((c/w))*(u-0.5);
    x=b+y;
    if (x>0.0){
      z=64.0*w*w*w*v*v;
      z=64.0*w*w*w*v*v;
      if (z<=1-2*y*y/x)
        acc=1;
      else if (log(z)<=2*(b*log(x/b)-y))
        acc=1;
    }
  }
  while (acc==0);
  return x;
}


double Rgamma(double a, double lambda)
{
  if (a<1.0) return gamsmall(a)/lambda;
  else if (a>1.0) return gambig(a)/lambda;
  else return expo(1.0);
}


/*--------------------------------------*/

/* Dirichlet random generator
   a and b are arrays of length k, containing doubles.
   a is the array of parameters
   b is the output array, where b ~ Dirichlet(a)  
*/

void RDirichlet(const double * a, const int k, double * b)
{
  int i;
  double sum=0.0;
  for(i=0;i<k;i++)
    {
      b[i]=Rgamma(a[i],1);
      sum += b[i];
    }
  for(i=0;i<k;i++)
    {
      b[i] /= sum;
    }
}


double RBeta(const double alpha, const double beta) {
  double in[2], out[2];
  in[0]=alpha;
  in[1]=beta;
  RDirichlet(in, 2, out);
  return out[0];
}


/*---------------------------------------*/


double Poisson(double xm)
{
  double gammLn(double xx);
  static double sq,alxm,g,oldm=(-1.0);
  double em, t, y;
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      ++em;
      t *=uniform();
    } while (t>g);
  } 
  else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammLn(xm+1.0);
    }
    do {
      do {
	y=tan(PI*uniform());
	em=sq*y+xm;
      } while (em< 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammLn(em+1.0)-g);
    } while (uniform()>t);
  }
  return em;
}


