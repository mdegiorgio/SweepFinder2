#include <math.h>
#include "SweepFinder2.h"
#include "bfgs.h"


static double *upbound=NULL, *lowbound=NULL;
static int numpar;

struct rdatatype {
  int x, n, folded, count;
  struct rdatatype *next;
};
static struct rdatatype *rdata=NULL;
static int rdatasize=0;


void free_rdata() {
  struct rdatatype *curr, *next;
  for (curr=rdata; curr!=NULL; curr=next) {
    next=curr->next;
    free(curr);
  }
  rdata=NULL;
  rdatasize=0;
}

void get_rnlim() {
  struct rdatatype *curr=rdata;
  nmax=nmin=0;
  if (curr==NULL) {
    printf("get_rnlim: curr is NULL!\n"); fflush(stdout);
    return;
  }
  nmax=nmin=curr->n;
  for (curr=curr->next; curr!=NULL; curr=curr->next) {
    if (curr->n > nmax) nmax=curr->n;
    if (curr->n < nmin) nmin=curr->n;
  }
}

int rdata_match(struct datatype data, struct rdatatype rdata) {
  if (data.x!=rdata.x) return data.x-rdata.x;
  if (data.n!=rdata.n) return data.n-rdata.n;
  if (data.folded!=rdata.folded) return data.folded-rdata.folded;
  return 0;
}


void add_rdata(struct datatype data) {
  struct rdatatype *curr, *newnode, *prev=NULL;
  for (curr=rdata; curr!=NULL; curr=curr->next) {
    if (rdata_match(data, *curr) >=0) break;
    prev = curr;
  }
  if (curr!=NULL) {
    if (rdata_match(data, *curr)==0) {
      (curr->count)++;
      return;
    }
  }
  newnode = malloc(sizeof(struct rdatatype));
  newnode->x = data.x;
  newnode->n = data.n;
  newnode->folded = data.folded;
  newnode->count=1;
  newnode->next=curr;
  if (curr==rdata)
    rdata=newnode;
  else prev->next=newnode;
  rdatasize++;
}


void remove_rdata(struct datatype data) {
  struct rdatatype *curr, *prev=NULL;
  for (curr=rdata; curr!=NULL; curr=curr->next) {
    if (rdata_match(data, *curr)==0) break;
    prev=curr;
  }
  if (curr==NULL) {
    fprintf(stderr, "remove_rdata: could not find data!\n"); 
    fprintf(stderr, "x=%i n=%i folded=%i\n",
	    data.x, data.n, data.folded);
    exit(-1);}
  if (curr->count > 1) {
    (curr->count)--;
    return;
  }
  if (curr==rdata) rdata=curr->next;
  else prev->next = curr->next;
  free(curr);
  rdatasize--;
}
  

void print_rdata(char *outfn) {
  FILE *outfile=my_fopen(outfn, "w");
  struct rdatatype *curr;
  for (curr=rdata; curr!=NULL; curr=curr->next) 
    fprintf(outfile, "%i\t%i\t%i\t%i\n",
	    curr->x, curr->n, curr->folded, curr->count);
  fclose(outfile);
}


void get_rdata(char *outfn) {
  int i;
  free_rdata();
  rdata=NULL;
  for (i=0; i<datasize; i++) 
    add_rdata(data[i]);
  if (outfn!=NULL)
    print_rdata(outfn);
  fflush(stdout);
}
	

void load_rdata(char *fn) {
  FILE *infile;
  struct rdatatype *newnode, *prev=NULL;
  int n, x, f, c;
  rdatasize=0;
  infile = my_fopen(fn, "r");
  printf("loading %s\n", fn); fflush(stdout);
  free_rdata();
  while (EOF!=fscanf(infile, "%i %i %i %i", &x, &n, &f, &c)) {
    if (rdatasize==0) {
      while ('\n'!=(c=fgetc(infile))) if (!isspace(c)) break;
    }
    newnode = malloc(sizeof(struct rdatatype));
    newnode->x=x;
    newnode->n=n;
    newnode->folded=f;
    newnode->count=c;
    newnode->next=NULL;
    if (rdata==NULL) rdata=newnode;
    else prev->next=newnode;
    prev=newnode;
    rdatasize++;
  }
  fflush(stdout);
  fclose(infile);
  printf("done loaddata rdatasize=%i\n", rdatasize);
} 
      

void makeguess(double *p, char *fn) {
  int i, mini=1;
  double sum=0.0;
  struct rdatatype *curr;

  for (i=0; i<xmax; i++) p[i]=0.0;
  for (curr=rdata; curr!=NULL; curr=curr->next)
    p[curr->x] += curr->count;
  if (invar==2) mini=0;
  for (i=mini; i<xmax; i++) {
    if (p[i]==0.0) p[i]=1.0;
    sum += p[i];
  }
  for (i=0; i<xmax; i++) p[i] /= sum;
  return;
}


void get_p(double *p, const double *invec) {
  double sum=0.0;
  int i, pos;
  if (invar==2) pos=0;
  else pos=1;
  for (i=0; i<xmax; i++) p[i]=0.0;
  for (i=0; i<numpar; i++) 
    sum += (p[i+pos+1] = exp(invec[i]));
  sum++;
  p[pos]=1.0;
  for (i=0; i<xmax; i++) p[i] /= sum;
  fflush(stdout);
}


void get_invec(double *invec, const double *p) {
  int pos, i;
  double norm;
  if (invar==2) pos=0;
  else pos=1;
  if (p[pos]<=0.0) {
    fprintf(stderr, "get_invec pos=%i p[%i]=%e\n", pos, pos, p[pos]);
    exit(-1);
  }
  norm = log(p[pos]);
  for (i=0; i<numpar; i++) {
    if (p[pos+1+i]<=0.0) {
      fprintf(stderr, "get_invec pos=%i i=%i p[%i]=%e\n", pos, i,
	      pos+i+1, p[pos+i+1]);
      exit(-1);
    }
    invec[i] = log(p[pos+1+i])-norm;
  }
}


double likelihood_freq_onesnp(int x, int n, int folded, 
			      int count, double *p) {
  int first=1, minx, maxx, xi;
  double pr=0.0;

 likelihood_freq_calc:
  minx = x;
  maxx = x + nmax-n;
  for (xi=minx; xi<=maxx; xi++) {
    pr += p[xi]*exp(ln_xchoosey(xi, x)+
		    ln_xchoosey(nmax-xi, n-x)-
		    ln_xchoosey(nmax, n));
  }
  if (folded && first) {
    first=0;
    x = n-x;
    if (x!=n-x)
      goto likelihood_freq_calc;
  }
  return log(pr)*count;
}


double likelihood_freq_one_rdata(struct rdatatype *curr, double *p) {
  return likelihood_freq_onesnp(curr->x, curr->n, curr->folded, 
				curr->count, p);
}

       
double likelihood_freq(double *p) {
  double like=0.0;
  struct rdatatype *curr;

  for (curr=rdata; curr!=NULL; curr=curr->next) 
    like += likelihood_freq_one_rdata(curr, p);
  return like;
}


double p_likelihood_freq(const double *invec) {
  double like, *p = malloc(xmax*sizeof(double));
  get_p(p, invec);
  like = likelihood_freq(p);
  free(p);
  //  if (rank==0) printf("like=%f\n", like);
  return -like;
}


void d_likelihood_freq(const double *invec, double *outvec) {
  static int *ng = NULL;
  static int curr_numpar=-1;
  int i;
  if (curr_numpar!=numpar) {
    if (ng!=NULL) free(ng);
    ng = malloc(numpar*sizeof(int));
    for (i=0; i<numpar; i++) ng[i]=1;
    curr_numpar=numpar;
  }
  assert(lowbound!=NULL && upbound!=NULL);
  getgradient(numpar, invec, ng, outvec, p_likelihood_freq, 
	      lowbound, upbound);
}

//assumes rdata already loaded and lims already set
void getfreq(char *outfn) {
  double *p, *invec;
  int i;
  FILE *outfile;

  get_rdata(NULL);
  printf("getfreq rdatasize=%i\n", rdatasize);
  if (rdatasize<=0) {
    fprintf(stderr, "getfreq_nompi rdataisze=%i\n", rdatasize);
    exit(-1);
  }
  p=malloc(xmax*sizeof(double));
  makeguess(p, NULL);
  if (invar==2) numpar=xmax-1;
  else numpar = xmax-2;
  lowbound=malloc(numpar*sizeof(double));
  upbound=malloc(numpar*sizeof(double));
  for (i=0; i<numpar; i++) {
    lowbound[i]=-15;
    upbound[i]=10;
  }
  invec=malloc(numpar*sizeof(double));
  get_invec(invec, p);
  for (i=0; i<numpar; i++) {
    if (invec[i] < lowbound[i]) invec[i]=lowbound[i];
    if (invec[i] > upbound[i]) invec[i]=upbound[i];
  }
  //  doNRinits(numpar);
  printf("calling findmax numpar=%i\n", numpar);
  findmax_bfgs(numpar, invec, p_likelihood_freq, d_likelihood_freq,
	       lowbound, upbound, NULL, 10);
  //  findmax(invec, numpar, p_likelihood_freq, d_likelihood_freq,
  //	  lowbound, upbound, 1);
  get_p(p, invec);
  outfile = my_fopen(outfn, "w");
  for (i=0; i<xmax; i++)
    fprintf(outfile, "%i\t%e\n", i, p[i]);
  fclose(outfile);
  free(invec);
  free(upbound);
  free(lowbound);
  free(p);
  //  doNRfree(numpar);
  lowbound=upbound=NULL;
  printf("frequency spectrum written to %s\n", outfn);
  free_rdata();
  return;
}


double *loadfreq(char *infn) {
  FILE *infile=my_fopen(infn, "r");
  double *p = malloc(xmax*sizeof(double)), val;
  int i;
  for (i=0; i<xmax; i++) p[i]=-1.0;
  p[0]=0.0;
  while (EOF!=fscanf(infile, "%i %le", &i, &val)) {
    if (val<0.0 || i<0 || i>=xmax) {
      fprintf(stderr, "error loading frequency spectrum from %s\n", infn);
      exit(-1);
    }
    p[i] = val;
  }
  for (i=0; i<xmax; i++) assert(p[i] >=0.0);
  fclose(infile);
  return p;
}
