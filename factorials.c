#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "factorials.h"

static double *ln_fac=NULL;
static int max_n=-1;

void get_factorials(int n) {
  int i, start;
  if (n > 1000000) {
    double *val=NULL;
    fprintf(stderr, "get_factorials: n=%i\n", n);
    val[100]=3.0;  //crash on purpose so we get a core dump
    exit(-1);
  }
  
  if (n==-1) {
    if (ln_fac!=NULL) {
      free(ln_fac);
      ln_fac=NULL;
    }
    return;
  }
  if (n <= max_n) return;
  if (max_n==-1) {
    if (n<3) n=3;
    ln_fac = malloc((n+1)*sizeof(double));
    ln_fac[0] = 0;
    start=1;
  }
  else {
    ln_fac = realloc(ln_fac, (n+1)*sizeof(double));
    start = max_n+1;
  }
  for (i=start; i<=n; i++)
    ln_fac[i] = ln_fac[i-1]+log((double)i);
  max_n = n;
}

double factorial(int n) {
  if (n <0) {
    fprintf(stderr, "factorial got %i\n", n);
    exit(-1);
  }
  if (n > max_n) get_factorials(n);
  return exp(ln_fac[n]);
}


double ln_factorial(int n) {
  if (n < 0) {
    fprintf(stderr, "ln_factorial: n=%i\n", n);
    exit(-1);
  }
  if (n > max_n)
    get_factorials(n);
  return ln_fac[n];
}

double xchoosey(int x, int y) {
  if (y < 0) {
    fprintf(stderr, "xchoosey got y=%i\n", y);
    exit(-1);
  }
  if (x > max_n) get_factorials(x);
  if (y > max_n) get_factorials(y);
  if (x<y) return 0;
  if (y==x || y==0) return 1;
  return exp(ln_fac[x]-ln_fac[x-y]-ln_fac[y]);
}

double ln_xchoosey(int x, int y) {
  if (x < 0 || y < 0) {
    fprintf(stderr, "ln_xchoosey: %i, %i\n", x, y);
    exit(-1);
  }
  if (x > max_n) get_factorials(x);
  if (y > max_n) get_factorials(y);
  if (y==x || y==0) return 0;
  if (y > x || x==0) {
    double *x1=NULL;
    fprintf(stderr, "ln_xchoosey: x=%i y=%i\n", x, y); 
    fprintf(stderr, "%f\n", x1[100000]);
    exit(-1);
  }
  return ln_fac[x]-ln_fac[x-y]-ln_fac[y];
}


double multinom(int x, int n, ...) {
  int i, y;
  va_list ap;
  double ret_val=0.0;
  
  va_start(ap, n);
  for (i=0; i<n; i++) {
    y = va_arg(ap, int);
    if (x < y) return 0.0;
    ret_val += ln_xchoosey(x, y);
    x -= y;
    if (x < 0) return 0.0;
  }
  va_end(ap);
  return exp(ret_val);
}
  
