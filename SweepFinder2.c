#include <stdarg.h>
#include <math.h>
#include "SweepFinder2.h"
#include "my_rand.h"
#include "freq.h"
#include "sort.h"
//#include "matrix.h"

//stores the probabilities over a grid of values of alpha*d
//indexed by prob[n][x][2][gridsize]
//gridsize is hard-coded, not same as argument to program
//prob[n][x][0][g] holds probabilities, 
//prob[n][x][1][g] holds vector that spline function needs
double ****prob=NULL;
double *ads=NULL; //holds vector of alpha*d values used to get prob

//stores the probabilities over a grid of values of alpha*d
//indexed by prob_bvalue[n][x][gridsize_bvalue][2][gridsize]
//gridsize_bvalue is hard-coded, 100 values uniformly between 0 and 1
//gridsize is hard-coded, not same as argument to program
//prob_bvalue[n][x][b][0][g] holds probabilities,
//prob_bvalue[n][x][b][1][g] holds vector that spline function needs
double *****prob_bvalue=NULL;

int invar=0;
int datasize=0, sweep_width, nmax, nmin, xmax;
double lowbound[MAXPAR], upbound[MAXPAR];
double minFreqPos=-1.0, maxFreqPos=-1.0;
struct datatype *data=NULL;

int gridsize_bvalue = 101; // Sets size of B-value grid
double *bvalue_grid; // B-value grid

double N_anc = -1.0; // IN DIPLOID INDIVIUALS (READ FROM INPUT)
double N_curr = -1.0; // IN DIPLOID INDIVIUALS (READ FROM INPUT)
double split_time = -1.0; // IN GENERATIONS (READ FROM INPUT)

void getlims() {
  int i;
  invar=0;
  if (datasize==0) {
    nmin=nmax=0;
    return;
  }
  nmin = nmax = data[0].n;
  minFreqPos = data[0].loc;
  maxFreqPos = minFreqPos;
  for (i=0; i<datasize; i++) {
    if(data[i].loc < minFreqPos) minFreqPos = data[i].loc;
    if(maxFreqPos < data[i].loc) maxFreqPos = data[i].loc;
    if (data[i].n>nmax) nmax=data[i].n;
    if (data[i].n<nmin) nmin=data[i].n;
    if (invar==0 && data[i].x==data[i].n) invar=1;
    if (invar!=2 && data[i].x==0 && data[i].n > 0) invar=2;
    if (invar!=2 && data[i].x==data[i].n && data[i].n > 0 && data[i].folded) invar=2;
  }
  if (invar) xmax=nmax+1;
  else xmax=nmax;
}

// variables for storing user-defined grid
int flagUserDefinedGrid = 0;
int datasize_grid=0;
double gridposmin, gridposmax;
double *data_grid=NULL;

// obtains the limits for the user-defined grid
void getlims_grid()
{
	int i;

	if(datasize_grid==0) {
		gridposmin=gridposmax=0.0;

		return;
	}

	gridposmin = data_grid[0];
	gridposmax = data_grid[0];

	for(i=0; i<datasize_grid; i++) {
		if(data_grid[i] > gridposmax) {
			gridposmax=data_grid[i];
		}

		if(data_grid[i]<gridposmin) {
			gridposmin=data_grid[i];
		}
	}
}

// variables for storing recombination map
int datasize_rec=0;
double recratemin, recratemax;
struct datatype_rec *data_rec=NULL;

// obtains the limits for the recombination map
void getlims_rec()
{
	int i;

	if(datasize_rec==0) {
		recratemin=recratemax=0.0;

		return;
	}

	recratemin = data_rec[0].rate;
	recratemax = data_rec[0].rate;

	for(i=0; i<datasize_rec; i++) {
		if(data_rec[i].rate>recratemax) {
			recratemax=data_rec[i].rate;
		}

		if(data_rec[i].rate<recratemin) {
			recratemin=data_rec[i].rate;
		}
	}
}

// variables for storing the b-value map
int datasize_bvalue=0;
double bvaluemin, bvaluemax;
struct datatype_bvalue *data_bvalue=NULL;

// obtains the limits for the recombination map
void getlims_bvalue()
{
	int i;

	if(datasize_bvalue==0) {
		bvaluemin=bvaluemax=0.0;

		return;
	}

	bvaluemin = data_bvalue[0].bvalue;
	bvaluemax = data_bvalue[0].bvalue;

	for(i=0; i<datasize_bvalue; i++) {
		if(data_bvalue[i].bvalue>bvaluemax) {
			bvaluemax=data_bvalue[i].bvalue;
		}

		if(data_bvalue[i].bvalue<bvaluemin) {
			bvaluemin=data_bvalue[i].bvalue;
		}
	}
}


FILE *my_fopen(char *fn, char *mode) {
  FILE *rv = fopen(fn, mode);
  if (rv==NULL) {
    fprintf(stderr, "error opening %s mode \"%s\"\n", fn, mode);
    exit(-1);
  }
  return rv;
}


void get_nlim() {
  int i;
  if (datasize==0) {
    nmax=nmin=0;
    return;
  }
  nmax=nmin=data[0].n;
  for (i=1; i<datasize; i++) {
    if (data[i].n > nmax) nmax=data[i].n;
    if (data[i].n < nmin) nmin=data[i].n;
  }
}


void get_lims() {
  fflush(stdout);
  get_nlim();
  if (invar) xmax = nmax+1;
  else xmax = nmax;
}

#define NCOLUMN 4
#define LOC 0
#define X 1
#define N 2
#define FOLDED 3

//all column headers except folded are required.
//if folded is not present, all SNPs will be assumed unfolded
char colname[NCOLUMN][10]={"position", "x", "n", "folded"};

//get rid of some data so we can test missing data
void losedata(double per) {
  struct datatype *newdata = malloc(datasize*sizeof(struct datatype));
  int newdatasize=0;
  int i,j , x, n;
  for (i=0; i<datasize; i++) {
    x = data[i].x;
    n = data[i].n;
    for (j=0; j<data[i].n; j++) {
      if (uniform()  < per) {
	if (j < data[i].x) x--;
	n--;
      }
    }
    if (x!=0 && x!=n) {
      newdata[newdatasize].x = x;
      newdata[newdatasize].n = n;
      newdata[newdatasize].folded=data[i].folded;
      newdata[newdatasize].loc = data[i].loc;
      newdatasize++;
    }
  }
  if (datasize > 0) free(data);
  datasize = newdatasize;
  data = newdata;
  data = realloc(newdata, datasize*sizeof(struct datatype));
  printf("done losedata datasize=%i\n", datasize);
}

void readms_error(char *infn) {
  printf("error reading data from %s\n", infn);
  exit(-1);
}


int readsnps_ms(char *infn) {
  static FILE *infile=NULL;
  static char *curr_infile=NULL;
  static int numchr=0, rep=0;
  char str[1000], c;
  int i, j;

  if (infn!=curr_infile) {
    if (infile!=NULL) 
      fclose(infile);
    printf("opening %s\n", infn);
    rep=0;
    infile=fopen(infn, "r");
    curr_infile=infn;
    if (EOF==fscanf(infile, "%s %i", str, &numchr)) readms_error(infn);
    if (strstr(str, "ms")==NULL || numchr <=0) readms_error(infn);
    printf("numchr=%i\n", numchr);
  }
  while (1) {
    if (EOF==fscanf(infile, "%s", str)) {
      fclose(infile);
      infile=NULL;
      curr_infile=NULL;
      return 0;
    }
    if (strcmp(str, "segsites:")==0) break;
  }
  if(EOF==fscanf(infile, "%i",&datasize)) readms_error(infn);
  printf("reading %s rep %i\n", infn, ++rep);
  printf("datasize=%i\n", datasize);
  if (datasize<=0) {
    printf("readsnps_ms: segsites must be greater than 0!\n");
    readms_error(infn);
  }
  data = malloc(datasize*sizeof(struct datatype));
  if(EOF==fscanf(infile, "%s", str)) readms_error(infn);
  if(strcmp(str, "positions:")!=0) readms_error(infn);
  for (i=0; i<datasize; i++) {
    if (EOF==fscanf(infile, "%lf", &data[i].loc))
      readms_error(infn);
    data[i].folded=0;
    data[i].x=0;
    data[i].n=numchr;
  }
  for (i=0; i<numchr; i++) {
    while ('\n'!=(c=fgetc(infile))) if (isspace(c)==0) readms_error(infn);
    for (j=0; j<datasize; j++) {
      c=fgetc(infile);
      if (c!='0' && c!='1') readms_error(infn);
      data[j].x += (c=='1');
    }
  }
  for (i=0; i<datasize; i++) {
    if (data[i].x==0 || data[i].x>=data[i].n) {
      printf("invalid data in msfile %s data[%i].x=%i n=%i\n", infn,
	     i, data[i].x, data[i].n);
      readms_error(infn);
    }
  }
  getlims();
  return 1;
}

  
  
      
      


//this is the format the old version of yuseob's program outputs
//it doesn't output n anywhere so you have to know what it is
void readsnps_yuseob(char *infn, int rep, int n) {
  FILE *infile = my_fopen(infn, "r");
  char str[1000], c;
  int count=0, nsites, i;
  if (datasize > 0) free(data);
  datasize=0;
  while (1) {
    while (EOF!=fscanf(infile, "%s", str))
      if (strcmp(str, "sites")==0) break;
    if (strcmp(str, "sites")!=0) break;
    while (':'!=(c=fgetc(infile))) assert(isspace(c));
    fscanf(infile, "%i", &nsites);
    if (rep < 0 || count==rep) {
      if (datasize==0) data = malloc(nsites*sizeof(struct datatype));
      else data = realloc(data, (datasize+nsites)*sizeof(struct datatype));
      for (i=0; i<nsites; i++) {
	assert(EOF!=fscanf(infile, "%lf %i", &data[datasize+i].loc,
			   &data[datasize+i].x));
	data[datasize+i].folded=0;
	data[datasize+i].n=n;
      }
      datasize += nsites;
      if (count==rep) break;
    }
    count++;
  }
  fclose(infile);
  losedata(0.1);
  getlims();
  printf("done getlims_yuseob datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n",
	 datasize, nmax, nmin, xmax, invar);
}


void printdata(char *outfn) {
  int i;
  FILE *outfile = fopen(outfn, "w");
  fprintf(outfile, "position\tx\tn\tfolded\n");
  for (i=0; i<datasize; i++)
    fprintf(outfile, "%f\t%i\t%i\t%i\n", data[i].loc,
	    data[i].x, data[i].n, data[i].folded);
  fclose(outfile);
}

// prints recombination map
void printdata_rec(char *outfn)
{
	int i;
	FILE *outfile = fopen(outfn, "w");

	fprintf(outfile, "position\trate\n");

	for(i=0; i<datasize_rec; i++) {
		fprintf(outfile, "%f\t%le\n", data_rec[i].loc,data_rec[i].rate);
	}

	fclose(outfile);
}

// prints b-value map
void printdata_bvalue(char *outfn)
{
	int i;
	FILE *outfile = fopen(outfn, "w");

	fprintf(outfile, "position\tbvalue\n");

	for(i=0; i<datasize_bvalue; i++) {
		fprintf(outfile, "%f\t%le\n", data_bvalue[i].loc,data_bvalue[i].bvalue);
	}

	fclose(outfile);
}


void readsnps_error(char *infn) {
  printf("error reading %s\n", infn);
  exit(-1);
}


void readsnps(char *infn) {
  FILE *infile=my_fopen(infn, "r");
  int colpos[NCOLUMN], pos=0, col=0, i, j;
  char c, str[1000];
  if (datasize > 0) free(data);
  datasize=0;
  c=fgetc(infile);
  for (i=0; i<NCOLUMN; i++) colpos[i]=-1;
  while (c!='\n' && c!=EOF) {
    str[0]=c;
    pos=1;
    while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF)
      str[pos++]=c;
    str[pos]='\0';
    for (i=0; i<NCOLUMN; i++)
      if (strcmp(str, colname[i])==0) {
	colpos[i]=col;
	break;
      }
    col++;
    if (c=='\n' || c==EOF) break;
    c=fgetc(infile);
  }
  if (colpos[LOC]<0 || colpos[X]<0 || colpos[N]<0) {
    fprintf(stderr, "readsnps: infile should have columns named position, x, and n (and optionally folded\n");
    exit(-1);
  }
  while (EOF!=(c=fgetc(infile)))
    if (c=='\n') datasize++;
  fclose(infile);
  infile=my_fopen(infn, "r");
  while ('\n'!=(c=fgetc(infile)) && c!=EOF);
  data = malloc(datasize*sizeof(struct datatype));
  for (i=0; i<datasize; i++) {
    for (j=0; j<col; j++) {
      if (colpos[LOC]==j) 
	{if(EOF==fscanf(infile, "%lf", &data[i].loc)) readsnps_error(infn);}
      else if (colpos[X]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].x)) readsnps_error(infn);}
      else if (colpos[N]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].n)) readsnps_error(infn);}
      else if (colpos[FOLDED]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].folded)) readsnps_error(infn);}
      else {
	while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF);
	if (c=='\n' || c==EOF) if(j!=col-1) readsnps_error(infn);
      }
    }
  }
  fclose(infile);
  if (colpos[FOLDED] < 0)
    for (i=0; i<datasize; i++) 
      data[i].folded=0;
  getlims();
  printf("done readsnps datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n", datasize, nmax, nmin, xmax, invar);
}

void readgrid(char *infn) { // Okay, we have read in the corret data file
  FILE *infile=my_fopen(infn, "r");
  int i;
  char c;
  if (datasize_grid > 0) free(data_grid);
  datasize_grid=0;
  c=fgetc(infile);
  while (EOF!=(c=fgetc(infile)))
    if (c=='\n') datasize_grid++;
  fclose(infile);
  infile=my_fopen(infn, "r");
  data_grid = malloc(datasize_grid*sizeof(double));
  for (i=0; i<datasize_grid; i++) {
    fscanf(infile, "%lf", &data_grid[i]);
  }
  fclose(infile);
  getlims_grid();

  printf("done readgrid datasize_grid=%i gridposmin=%le gridposmax=%le\n", datasize_grid, gridposmin, gridposmax);
}

// sets constants for recombination map
#define NCOLUMN_REC 2
#define LOC_REC 0
#define RATE_REC 1

// all recombination map column headers are required.
char colname_rec[NCOLUMN_REC][10]={"position", "rate"};

// reads the recombination map
void readrecmap(char *infn)
{
	FILE *infile=my_fopen(infn, "r");
	int colpos[NCOLUMN_REC], pos=0, col=0, i, j;
	char c, str[1000];

	if(datasize_rec > 0) {
		free(data_rec);
	}

	datasize_rec=0;
	c=fgetc(infile);

	for(i=0; i<NCOLUMN_REC; i++) {
		colpos[i]=-1;
	}

	while(c!='\n' && c!=EOF) {
		str[0]=c;
		pos=1;

		while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
			str[pos++]=c;
		}

		str[pos]='\0';

		for(i=0; i<NCOLUMN_REC; i++) {
			if(strcmp(str, colname_rec[i])==0) {
				colpos[i]=col;
				break;
			}
		}

		col++;

		if(c=='\n' || c==EOF) {
			break;
		}

		c=fgetc(infile);
	}

	if(colpos[LOC_REC]<0 || colpos[RATE_REC]<0) {
		fprintf(stderr, "readrecmap: infile should have columns named position and rate\n");
		exit(-1);
	}

	while(EOF!=(c=fgetc(infile))) {
		if(c=='\n') {
			datasize_rec++;
		}
	}

	fclose(infile);
	infile=my_fopen(infn, "r");

	while('\n'!=(c=fgetc(infile)) && c!=EOF) {
		;
	}

	data_rec = malloc(datasize_rec*sizeof(struct datatype_rec));

	for(i=0; i<datasize; i++) {
//	for(i=0; i<datasize_rec; i++) {
		for(j=0; j<col; j++) {
			if(colpos[LOC_REC]==j) {
				if(EOF==fscanf(infile, "%lf", &data_rec[i].loc)) {
					readsnps_error(infn);
				}
			}
			else if(colpos[RATE_REC]==j) {
				if(EOF==fscanf(infile, "%le", &data_rec[i].rate)) {
					readsnps_error(infn);
				}

				// If not 0, begin constructing cumulative rate
				if(i > 0) {
					// Impose a minimum recombination value 10^{-6}
			//		if(data_rec[i].rate == 0.0) {
			//			data_rec[i].rate = 1.0e-6;		
			//		}

					data_rec[i].rate = data_rec[i].rate + data_rec[i - 1].rate;
				}
			}
			else {
				while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
					;
				}

				if(c=='\n' || c==EOF) {
					if(j!=col-1) {
						readsnps_error(infn);
					}
				}
			}
		}
	}

	for(i=0; i<datasize_rec; i++) {
		data_rec[i].rate = data_rec[i].rate * 1e6;
	}

	fclose(infile);

	getlims_rec();

	printf("done readrecmap datasize_rec=%i recratemin=%le recratemax=%le\n", datasize_rec, recratemin, recratemax);
}

// sets constants for recombination map
#define NCOLUMN_BVALUE 2
#define LOC_BVALUE 0
#define VALUE_BVALUE 1

// all B value column headers are required.
char colname_bvalue[NCOLUMN_BVALUE][10]={"position", "bvalue"};

// reads the B-value map
void readbvalues(char *infn)
{
	FILE *infile=my_fopen(infn, "r");
	int colpos[NCOLUMN_BVALUE], pos=0, col=0, i, j;
	char c, str[1000];

	if(datasize_bvalue > 0) {
		free(data_bvalue);
	}

	datasize_bvalue=0;
	c=fgetc(infile);

	for(i=0; i<NCOLUMN_BVALUE; i++) {
		colpos[i]=-1;
	}

	while(c!='\n' && c!=EOF) {
		str[0]=c;
		pos=1;

		while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
			str[pos++]=c;
		}

		str[pos]='\0';

		for(i=0; i<NCOLUMN_BVALUE; i++) {
			if(strcmp(str, colname_bvalue[i])==0) {
				colpos[i]=col;
				break;
			}
		}

		col++;

		if(c=='\n' || c==EOF) {
			break;
		}

		c=fgetc(infile);
	}

	if(colpos[LOC_BVALUE]<0 || colpos[VALUE_BVALUE]<0) {
		fprintf(stderr, "readbvalues: infile should have columns named position and bvalue\n");
		exit(-1);
	}

	while(EOF!=(c=fgetc(infile))) {
		if(c=='\n') {
			datasize_bvalue++;
		}
	}

	fclose(infile);
	infile=my_fopen(infn, "r");

	while('\n'!=(c=fgetc(infile)) && c!=EOF) {
		;
	}

	data_bvalue = malloc(datasize_bvalue*sizeof(struct datatype_bvalue));

	for(i=0; i<datasize_bvalue; i++) {
		for(j=0; j<col; j++) {
			if(colpos[LOC_BVALUE]==j) {
				if(EOF==fscanf(infile, "%lf", &data_bvalue[i].loc)) {
					readsnps_error(infn);
				}
			}
			else if(colpos[VALUE_BVALUE]==j) {
				if(EOF==fscanf(infile, "%le", &data_bvalue[i].bvalue)) {
					readsnps_error(infn);
				}
			}
			else {
				while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
					;
				}

				if(c=='\n' || c==EOF) {
					if(j!=col-1) {
						readsnps_error(infn);
					}
				}
			}
		}
	}

	fclose(infile);

	getlims_bvalue();

	printf("done readbvalues datasize_bvalue=%i bvaluemin=%le bvaluemax=%le\n", datasize_bvalue, bvaluemin, bvaluemax);
}


//probability of choosing j from sample size s given frequency 
//spectrum p
double p_j_s(int j, int s, double *p) { // Calculation of pj,H
  int i;
  double rv=0.0;
  if (j < 0 || j > s) return 0.0;
  for (i=j; i<=nmax; i++) // ADVISED: Formula in paper goes to n-1, this formula goes to n
    if (s-j >=0 && i<xmax)
      rv += p[i]*xchoosey(i,j)*xchoosey(nmax-i,s-j)/
	xchoosey(nmax,s);
  return rv;
}

	 
double p_kescapes(int k, double ad) { // Calculation of Pe(k)
  double pe = exp(-ad); // Strangely, this is just the opposite of what is given in the paper
  return xchoosey(nmax, k)*pow((1.0-pe),k)*pow(pe,(nmax-k)); // But this reverses the opposite and makes it correct
}


double get_pstar_sub(int b0, double* p, double ad) { // Calculation of pB*
  int k, b;
  static double lastad=-1.0;
  static double *rv=NULL;

  if (lastad!=ad) {
    if (rv==NULL) rv = malloc((nmax+1)*sizeof(double));
    for (b=0; b<=nmax; b++) {
      rv[b] = 0.0;
      if (b < xmax)
	rv[b] = p[b]*p_kescapes(nmax, ad); // First term of pB* calculation
      for (k=0; k<nmax; k++) { // Last set of terms for pB* calculation
	rv[b] += p_kescapes(k, ad)*
	  (p_j_s(b+1-nmax+k, k+1, p)*(b+1-nmax+k)/(k+1) +
	   p_j_s(b, k+1, p)*(k+1-b)/(k+1));
      }
    }
    lastad = ad;
  }
  return rv[b0];
}


double get_pstar(int i, double* p, double ad) {
  double pr, sum;

  ad = ad + 1E-7; // Calculation becomes unstable when too close to the test site

  sum=1.0; // If sum < 1, then conditioning on only variable sites
  if (invar==0 || invar==1)
    sum -= get_pstar_sub(0, p, ad); // Subtract on invariable sites (except fixed differences)
  if (invar==0) 
    sum -= get_pstar_sub(nmax, p, ad); // Subtract out fixed differences
  assert(sum <= 1.0 && sum > 0.0);
  pr = get_pstar_sub(i, p, ad)/sum; // Calculate the pB*, with B=i, conditioning on the types of sites observed (variable or invariable)
  return pr;
}

double splint(double* xvec, double *yvec, double *yvec2, int n, 
	      double x)
{
  int lowpos,hipos,pos;
  double diff,b,a;
  double y;
  
  if (x < xvec[0]) return xvec[0];
  lowpos=0;
  hipos=n-1;
  while (hipos-lowpos > 1) {
    pos=(hipos+lowpos)/2;
    if (xvec[pos] > x) hipos=pos;
    else lowpos=pos;
  }
  diff=xvec[hipos]-xvec[lowpos];
  if (diff == 0.0) {fprintf(stderr, "splint error\n"); exit(-1);}
  a=(xvec[hipos]-x)/diff;
  b=(x-xvec[lowpos])/diff;
  y=a*yvec[lowpos]+b*yvec[hipos]+((a*a*a-a)*yvec2[lowpos]+(b*b*b-b)*yvec2[hipos])*(diff*diff)/6.0;
  if (y<=0.0 || y>=1.0) {
    y=(yvec[lowpos]*a+b*yvec[hipos]);
  }
  if (y <0.0 || y>1.0) {
    printf("y=%e lowpos=%i hipos=%i a=%e b=%e yvec=%e, %e\n", 
	   y, lowpos, hipos, a, b, yvec[lowpos], yvec[hipos]);
    y = (yvec[lowpos]+yvec[hipos])/2.0;
  }
  return y;
}


void spline(double *x, double *y ,int n, double yp1, double ypn,
	    double *y2)
{
  int i,k;
  double p,qn,sig,un,*u;
  
  /*	for (i=1; i<=n; i++) {
	
  printf("%i\t%e\t%e\t%e\n", i, x[i], y[i], y2[i]);
  }
  exit(-1);*/
  
  u=malloc(n*sizeof(double));
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    if (x[1]-x[0]==0.0) {
      printf("spline: x[1]=x[0]=%e\n", x[1]); exit(-1);}
    
  }
  for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    if (x[i+1]==x[i-1]) {
      printf("spline2: x[%i]=x[%i]=%e\n", i+1, i-1, x[i+1]);
      exit(-1);}
    p=sig*y2[i-1]+2.0;
    if (p==0.0) {
      printf("spline: p=%e\n", p); exit(-1);}
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    if (x[i+1]==x[i]) {
      printf("spline3: x[%i]=x[%i]=%e\n", i+1, i, x[i+1]); exit(-1);}
    if (x[i]==x[i-1]) {
      printf("spline4: x[%i]=x[%i]=%e\n", i, i-1, x[i]); exit(-1);}
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    printf("spline: Else!\n"); exit(-1);
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  if (qn*y2[n-2]+1.0==0.0) {
    printf("spline: qn=%e y2[%i]=%e\n", qn, n-2, y2[n-2]); exit(-1);}
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}


void calcprobs(double *p, int gridsize) {
  double pr, maxval, minval, interval;
  int i, n, x, j, k, minx, maxx, maxxi, xi, numlog, numrest;
  printf("calcprobs nmax=%i nmin=%i xmax=%i invar=%i\n",
	 nmax, nmin, xmax, invar);
  fflush(stdout);

  ads = malloc(gridsize*sizeof(double));
  numlog = gridsize*0.5;
  numrest = gridsize-numlog;
  //this is a bit ad-hoc but i'm choosing values of ad that cover
  //the entire range of minval to maxval, but concentrate on the higher
  //values which are usually more relevant
  //so first 50% of values are chosen on a log scale, other rest on linear
  minval=1.0e-8;
  maxval=BIG_ENOUGH;
  interval=log(maxval/minval)/(double)(numlog-1);
  for (i=0; i<numlog; i++)
    ads[i]=exp(log(minval)+i*interval);
  minval=0.1;
  maxval=3.0;
  interval=(maxval-minval)/(double)(numrest-1);
  for (i=0; i<numrest; i++)
    ads[i+numlog]=minval+interval*i;
  dsort(gridsize, ads-1);
  for (i=1; i<gridsize; i++)
    if (ads[i]==ads[i-1]) {
      printf("ads[%i]=ads[%i]=%e\n\n\n", i, i-1, ads[i]);
      exit(-1);
    }

  prob = malloc((nmax+1)*sizeof(double***));
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    prob[n] = malloc(maxx*sizeof(double**));
    for (x=minx; x<maxx; x++) {
      prob[n][x] = malloc(2*sizeof(double*));
      for (j=0; j<2; j++) {
	prob[n][x][j] = malloc(gridsize*sizeof(double));
	for (k=0; k<gridsize; k++) 
	  prob[n][x][j][k] = 0.0;
      }
    }
  }
  for (k=0; k<gridsize; k++) {
    fflush(stdout);
    for (n=nmin; n<=nmax; n++) {
      if (invar==0) {minx=1; maxx=n;}
      else if (invar==1) {minx=1; maxx=n+1;}
      else {minx=0; maxx=n+1;}
      for (x=minx; x<maxx; x++) {
	maxxi = x + nmax - n;
	for (xi=x; xi<=maxxi; xi++) {
	  pr = exp(ln_xchoosey(xi, x)+
		   ln_xchoosey(nmax-xi, n-x)-
		   ln_xchoosey(nmax, n));
	  prob[n][x][0][k] +=
	    pr*get_pstar(xi, p, ads[k]);
	}
      }
    }
  }
  
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    for (x=minx; x<maxx; x++) {
      spline(ads, prob[n][x][0], gridsize, 1.0e30, 
	     1.0e30, prob[n][x][1]);
    }
  }
  printf("done calcprob\n");
  fflush(stdout);
}

// Modifies the frequency spectrum using B-values
void ModifySpectrumWithBvalues(double *p, int n, double bvalue, double *p_modified)
{
	int i;
	//double N_anc = 10000.0;//250.0;	 ///// TEMPORARY PARAMETERS, MAKE SURE TO DELETE THIS AND INSTEAD READ IN VALUES
	//double N_curr = 10000.0;//250.0;
	//double split_time = 250000.0;//2000.0;	
	double sum = 0.0;


	
	if(invar == 0) {
		for(i = 0; i < xmax; i++) {	
			p_modified[i] = p[i]; // All polymorphic, so just use original spectrum		
		}
	}	
	else {
		if(invar == 1) {
			p_modified[0] = p[0];
		}
		else {
			sum = p[0];
			p_modified[0] = p[0];
		}

		// Polymorphism component
		for(i = 1; i < xmax - 1; i++) {
			p_modified[i] = bvalue * p[i];
			sum = sum + p_modified[i];
		}

		// Divergence component
		p_modified[xmax - 1] = p[xmax - 1] * (split_time+2.0*bvalue*(N_anc-N_curr*(1.0-1.0/nmax))) / (split_time+2.0*(N_anc-N_curr*(1.0-1.0/nmax)));
		//p_modified[xmax - 1] = p[xmax - 1] * (split_time + 2.0 * bvalue * (N_anc - N_curr*(1.0 - 1.0/n))) / (split_time + 2.0 * (N_anc - N_curr*(1.0 - 1.0/n)));
		sum = sum + p_modified[xmax - 1];

		assert(sum > 0.0);

		// Normalize
		for(i = 0; i < xmax; i++) {
			p_modified[i] = p_modified[i] / sum;
			assert(p_modified[i] <= 1.0);		
		}
	}
}

void calcprobs_with_bvalues(double *p, int gridsize) {
  double pr, maxval, minval, interval;
  int i, n, x, j, k, minx, maxx, maxxi, xi, numlog, numrest, bindex; // Added bindex
  double *p_modified = malloc(xmax*sizeof(double)); // make new allele frequency spectrum

  // Initialize to modified frequency spectrum
  for(i = 0; i < xmax; i++) {
    //p_modified[i]=p[i];
    p_modified[i] = 0.0;
  }

  printf("calcprobs_with_bvalues nmax=%i nmin=%i xmax=%i invar=%i\n",
	 nmax, nmin, xmax, invar);
  fflush(stdout);

  ads = malloc(gridsize*sizeof(double));
  numlog = gridsize*0.5;
  numrest = gridsize-numlog;
  //this is a bit ad-hoc but i'm choosing values of ad that cover
  //the entire range of minval to maxval, but concentrate on the higher
  //values which are usually more relevant
  //so first 50% of values are chosen on a log scale, other rest on linear
  minval=1.0e-8;
  maxval=BIG_ENOUGH;
  interval=log(maxval/minval)/(double)(numlog-1);
  for (i=0; i<numlog; i++)
    ads[i]=exp(log(minval)+i*interval);
  minval=0.1;
  maxval=3.0;
  interval=(maxval-minval)/(double)(numrest-1);
  for (i=0; i<numrest; i++)
    ads[i+numlog]=minval+interval*i;
  dsort(gridsize, ads-1);
  for (i=1; i<gridsize; i++)
    if (ads[i]==ads[i-1]) {
      printf("ads[%i]=ads[%i]=%e\n\n\n", i, i-1, ads[i]);
      exit(-1);
    }

  prob_bvalue = malloc((nmax+1)*sizeof(double****)); // Changed to prob_bvalue and inserted an extra pointer
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    prob_bvalue[n] = malloc(maxx*sizeof(double***)); // Changed to prob_bvalue and inserted an extra pointer
    for (x=minx; x<maxx; x++) {
      prob_bvalue[n][x] = malloc(gridsize_bvalue*sizeof(double**)); // Changed to prob_bvalue, inserted gridsize_bvalue and an extra pointer
      for(bindex = 0; bindex < gridsize_bvalue; bindex++) { // Added this loop
        prob_bvalue[n][x][bindex] = malloc(2*sizeof(double*)); // Add this memory allocation
        for (j=0; j<2; j++) {
	  prob_bvalue[n][x][bindex][j] = malloc(gridsize*sizeof(double)); // Changed to prob_bvalue and inserted [bindex]
	  for (k=0; k<gridsize; k++) 
	    prob_bvalue[n][x][bindex][j][k] = 0.0; // Changed to prob_bvalue and inserted [bindex]
        }
      }
    }
  }
  printf("Populating background spectra with B-values\n");
  for(bindex = 0; bindex < gridsize_bvalue; bindex++) { // Added this loop
    //printf("considering bvalue = %lf\n", bvalue_grid[bindex]); // Print B-values to consider
    ModifySpectrumWithBvalues(p, n, bvalue_grid[bindex], p_modified); // MODIFY FREQUENCY SPECTRUM TO INCLUDE BVALUES
    for (k=0; k<gridsize; k++) {
      fflush(stdout);
      for (n=nmin; n<=nmax; n++) {
        if (invar==0) {minx=1; maxx=n;}
        else if (invar==1) {minx=1; maxx=n+1;}
        else {minx=0; maxx=n+1;}
        for (x=minx; x<maxx; x++) {
	  maxxi = x + nmax - n;
	  for (xi=x; xi<=maxxi; xi++) {
	    pr = exp(ln_xchoosey(xi, x)+
		   ln_xchoosey(nmax-xi, n-x)-
		   ln_xchoosey(nmax, n));
	    prob_bvalue[n][x][bindex][0][k] += pr*get_pstar(xi, p_modified, ads[k]); // CALCULATE PROBABILITIES WITH MODIFIED SPECTRUM 
	  }
        }
      }
    }
  }
  
  for(bindex = 0; bindex < gridsize_bvalue; bindex++) { // Added this loop
    for (n=nmin; n<=nmax; n++) {
      if (invar==0) {minx=1; maxx=n;}
      else if (invar==1) {minx=1; maxx=n+1;}
      else {minx=0; maxx=n+1;}
      for (x=minx; x<maxx; x++) {
        spline(ads, prob_bvalue[n][x][bindex][0], gridsize, 1.0e30, 
	     1.0e30, prob_bvalue[n][x][bindex][1]); // Changed to prob_bvalue and inserted [bindex]
      }
    }
  }
  free(p_modified); // Free temporary Bvalue-corrected frequency spectrum
  printf("done calcprob_with_bvalues\n"); // Added _with_bvalues
  fflush(stdout);
}


//actually returns LR now
double ln_likelihood(double *p, double alpha, double sweep) {
  int s, n, x, gridsize=500, first;
  double like=0.0, dist, ad, pr;

  if (prob==NULL)
    calcprobs(p, gridsize);

  sweep_width=0;
  for (s=0; s<datasize; s++) {
    dist = fabs(data[s].loc-sweep);
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) continue;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0;
    first=1;
  calc:
    if (ad < ads[0]) 
      pr += prob[n][x][0][0];
    else 
      pr += splint(ads, prob[n][x][0],
		   prob[n][x][1], 
		   gridsize, ad);
    
    if (pr >= 0.0 && pr < 1.0);
    else {
      //printf("1: val=%e\n", pr);  // FOR DEBUGGING
      //printf("%i %i %e\n", nmax-n, x, ad);  // FOR DEBUGGING
      //printf("%e %e\t%e %e\n", ads[gridsize-2], 
	  //   prob[n][x][0][gridsize-2],
	  //   ads[gridsize-1],
	  //   prob[n][x][0][gridsize-1]); // FOR DEBUGGING
      //      debug_splint=1;
      splint(ads, prob[n][x][0],
	     prob[n][x][1],
	     gridsize, ad);
      
    }    
    if (pr < 0.0 || pr >= 1.00000001) {
        if(pr >= 1.00000001) { // Mod for errors in numerical precision
          pr = 1.0;
        } 
        else {
          printf("val=%e\n", pr); exit(-1);
        }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc;
    }
    like += log(pr) - data[s].baselike;
  }
  return like;
}

// For a given b-value, find an appropriate bin (grid index)
int GetIndexForBvalue(double bvalue)
{
	int i;

	for(i = gridsize_bvalue - 1; i >= 0; i--) {
		if(bvalue >= bvalue_grid[i]) {
			return i;
		}
	} 
	printf("B-value found that outside bvalue_grid range\n");//// MAKE THIS MESSAGE A LITTLE BETTER
	exit(-1);
}

//actually returns LR now
// Uses B-values
double ln_likelihood_with_bvalues(double *p, double alpha, double sweep) {
  int s, n, x, gridsize=500, first, b; // Added variable b
  double like=0.0, dist, ad, pr;

  if (prob_bvalue==NULL)
    calcprobs_with_bvalues(p, gridsize); // Changed to calcprobs_with_bvalues

  sweep_width=0;
  for (s=0; s<datasize; s++) {
    dist = fabs(data[s].loc-sweep);
    // Assuming 1 cM/Mb
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) continue;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0;
    first=1;

    b = GetIndexForBvalue(data_bvalue[s].bvalue); // Find B-value index

  calc:
    if (ad < ads[0]) 
      pr += prob_bvalue[n][x][b][0][0]; // Changed to prob_bvalue and inserted [b]
    else 
      pr += splint(ads, prob_bvalue[n][x][b][0],
		   prob_bvalue[n][x][b][1], 
		   gridsize, ad); // Changed to prob_bvalue and inserted [b] 
    
    if (pr >= 0.0 && pr < 1.0);
    else {
   //   printf("1: val=%e\n", pr);  // FOR DEBUGGING
   //   printf("%i %i %e\n", nmax-n, x, ad);  // FOR DEBUGGING
   //   printf("%e %e\t%e %e\n", ads[gridsize-2], 
//	     prob_bvalue[n][x][b][0][gridsize-2],
	//     ads[gridsize-1],
	//     prob_bvalue[n][x][b][0][gridsize-1]); // Changed to prob_bvalue and inserted [b]
      //      debug_splint=1;
      splint(ads, prob_bvalue[n][x][b][0],
	     prob_bvalue[n][x][b][1],
	     gridsize, ad); // Changed to prob_bvalue and inserted [b]
      
    }    
    if (pr < 0.0 || pr >= 1.00000001) {
      if(pr >= 1.00000001) { // Mod for errors in numerical precision
        pr = 1.0;
      }
      else {
        printf("val=%e\n", pr); exit(-1);
      }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc;
    }
    like += log(pr) - data[s].baselike;
  }
  return like;
}

//actually returns LR now
//Uses genetic distance instead of physical location
double ln_likelihood_with_map(double *p, double alpha, double sweep_genetic_distance) { // modified to sweep_genetic_distance
  int s, n, x, gridsize=500, first;
  double like=0.0, dist, ad, pr;

  if (prob==NULL)
    calcprobs(p, gridsize);

  sweep_width=0;
  for (s=0; s<datasize; s++) {
    dist = fabs(data_rec[s].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    // Assuming 1 cM/Mb
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) continue;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0;
    first=1;
  calc:
    if (ad < ads[0]) 
      pr += prob[n][x][0][0];
    else 
      pr += splint(ads, prob[n][x][0],
		   prob[n][x][1], 
		   gridsize, ad);
    
    if (pr >= 0.0 && pr < 1.0);
    else {
  //    printf("1: val=%e\n", pr); // FOR DEBUGGING
   //   printf("%i %i %e\n", nmax-n, x, ad); // FOR DEBUGGING
   //   printf("%e %e\t%e %e\n", ads[gridsize-2], 
	//     prob[n][x][0][gridsize-2],
	//     ads[gridsize-1],
	//     prob[n][x][0][gridsize-1]); // FOR DEBUGGING
      //      debug_splint=1;
      splint(ads, prob[n][x][0],
	     prob[n][x][1],
	     gridsize, ad);
      
    }    
    if (pr < 0.0 || pr >= 1.00000001) {
      if(pr >= 1.00000001) { // Mod for errors in numerical precision
        pr = 1.0;
      }
      else {
        printf("val=%e\n", pr); exit(-1);
      }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc;
    }
    like += log(pr) - data[s].baselike;
  }
  return like;
}

//actually returns LR now
//Uses genetic distance instead of physical locatio
// Uses B-values
double ln_likelihood_with_map_and_bvalues(double *p, double alpha, double sweep_genetic_distance) { // modified to sweep_genetic_distance
  int s, n, x, gridsize=500, first, b; // Added variable b
  double like=0.0, dist, ad, pr;

  if (prob_bvalue==NULL)
    calcprobs_with_bvalues(p, gridsize); // Changed to calcprobs_with_bvalues

  sweep_width=0;
  for (s=0; s<datasize; s++) {
    dist = fabs(data_rec[s].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    // Assuming 1 cM/Mb
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) continue;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0;
    first=1;

    b = GetIndexForBvalue(data_bvalue[s].bvalue); // Find B-value index

  calc:
    if (ad < ads[0]) 
      pr += prob_bvalue[n][x][b][0][0]; // Changed to prob_bvalue and inserted [b]
    else 
      pr += splint(ads, prob_bvalue[n][x][b][0],
		   prob_bvalue[n][x][b][1], 
		   gridsize, ad); // Changed to prob_bvalue and inserted [b] 
    
    if (pr >= 0.0 && pr < 1.0);
    else {
   //   printf("1: val=%e\n", pr); // FOR DEBUGGING
   //   printf("%i %i %e\n", nmax-n, x, ad); // FOR DEBUGGING
   //   printf("bvalue=%lf\n", data_bvalue[s].bvalue); // FOR DEBUGGING
    //  printf("%e %e\t%e %e\n", ads[gridsize-2], 
	 //    prob_bvalue[n][x][b][0][gridsize-2],
	  //   ads[gridsize-1],
	  //   prob_bvalue[n][x][b][0][gridsize-1]); // Changed to prob_bvalue and inserted [b]
      //      debug_splint=1;
      splint(ads, prob_bvalue[n][x][b][0],
	     prob_bvalue[n][x][b][1],
	     gridsize, ad); // Changed to prob_bvalue and inserted [b]
      
    }    
    if (pr < 0.0 || pr >= 1.00000001) {
      if(pr >= 1.00000001) { // Mod for errors in numerical precision
        pr = 1.0;
      }
      else {
        printf("val=%e\n", pr); exit(-1);
      }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc;
    }
    like += log(pr) - data[s].baselike;
  }
  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
double findmax_alpha(double *alpha, double *pi, double sweep) {
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double *vals, *likes, tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  vals = malloc(gridsize*sizeof(double));
  likes = malloc(gridsize*sizeof(double));
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data[i].loc-sweep) < mind)
      mind = fabs(data[i].loc-sweep);
    if (i==0 || fabs(data[i].loc-sweep) > maxd)
      maxd = fabs(data[i].loc-sweep);
    totd += fabs(data[i].loc-sweep);
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
  minalpha = BIG_ENOUGH/(totd/datasize);
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);

  // Ensures that the search is fine enough for a strong sweep
  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0)
	likes[i] = likes[i-1];
      else likes[i] = ln_likelihood(pi, vals[i], sweep);
      //      printf("%i %f %e like=%f sweep_width=%i\n", i, sweep, vals[i], likes[i], sweep_width);
      sw[i] = sweep_width;
      if (i==0 || likes[i] > likes[max])
	max = i;
    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals[max-1];
    if (max==gridsize-1) maxalpha += (vals[max-1]-vals[max-2]);
    else maxalpha = vals[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }
  like = likes[max];
  sweep_width = sw[max];
  *alpha = vals[max];
  free(vals);
  free(likes);
  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
// modifies spectrum using bvalues
double findmax_alpha_with_bvalues(double *alpha, double *pi, double sweep) {
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double *vals, *likes, tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  vals = malloc(gridsize*sizeof(double));
  likes = malloc(gridsize*sizeof(double));
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data[i].loc-sweep) < mind)
      mind = fabs(data[i].loc-sweep);
    if (i==0 || fabs(data[i].loc-sweep) > maxd)
      maxd = fabs(data[i].loc-sweep);
    totd += fabs(data[i].loc-sweep);
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
  minalpha = BIG_ENOUGH/(totd/datasize);
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);

  // Ensures that the search is fine enough for a strong sweep
  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0)
	likes[i] = likes[i-1];
      else likes[i] = ln_likelihood_with_bvalues(pi, vals[i], sweep); // Now calling ln_likelihood_with_bvalues
      //      printf("%i %f %e like=%f sweep_width=%i\n", i, sweep, vals[i], likes[i], sweep_width);
      sw[i] = sweep_width;
      if (i==0 || likes[i] > likes[max])
	max = i;
    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals[max-1];
    if (max==gridsize-1) maxalpha += (vals[max-1]-vals[max-2]);
    else maxalpha = vals[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }
  like = likes[max];
  sweep_width = sw[max];
  *alpha = vals[max];
  free(vals);
  free(likes);
  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
// assumes distances are genetic distances rather than physical
double findmax_alpha_with_map(double *alpha, double *pi, double sweep_genetic_distance) { // Modified to sweep_genetic_distance
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double *vals, *likes, tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  vals = malloc(gridsize*sizeof(double));
  likes = malloc(gridsize*sizeof(double));
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data_rec[i].rate-sweep_genetic_distance) < mind) // modified to data_rec[s].rate-sweep_genetic_distance
      mind = fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    if (i==0 || fabs(data_rec[i].rate-sweep_genetic_distance) > maxd) // modified to data_rec[s].rate-sweep_genetic_distance
      maxd = fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    totd += fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
  minalpha = BIG_ENOUGH/(totd/datasize);
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);

  // Ensures that the search is fine enough for a strong sweep
  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0)
	likes[i] = likes[i-1];
      else likes[i] = ln_likelihood_with_map(pi, vals[i], sweep_genetic_distance);  // Now calling ln_likelihood_with_map and sweep_genetic_distance
      //      printf("%i %f %e like=%f sweep_width=%i\n", i, sweep_genetic_distance, vals[i], likes[i], sweep_width); // Changed to sweep_genetic_distance
      sw[i] = sweep_width;
      if (i==0 || likes[i] > likes[max])
	max = i;
    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals[max-1];
    if (max==gridsize-1) maxalpha += (vals[max-1]-vals[max-2]);
    else maxalpha = vals[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }
  like = likes[max];
  sweep_width = sw[max];
  *alpha = vals[max];
  free(vals);
  free(likes);
  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
// assumes distances are genetic distances rather than physical
// modifies epctrum using B values
double findmax_alpha_with_map_and_bvalues(double *alpha, double *pi, double sweep_genetic_distance) { // Modified to sweep_genetic_distance
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double *vals, *likes, tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  vals = malloc(gridsize*sizeof(double));
  likes = malloc(gridsize*sizeof(double));
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data_rec[i].rate-sweep_genetic_distance) < mind) // modified to data_rec[s].rate-sweep_genetic_distance
      mind = fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    if (i==0 || fabs(data_rec[i].rate-sweep_genetic_distance) > maxd) // modified to data_rec[s].rate-sweep_genetic_distance
      maxd = fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
    totd += fabs(data_rec[i].rate-sweep_genetic_distance); // modified to data_rec[s].rate-sweep_genetic_distance
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
  minalpha = BIG_ENOUGH/(totd/datasize);
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);

  // Ensures that the search is fine enough for a strong sweep
  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0)
	likes[i] = likes[i-1];
      else likes[i] = ln_likelihood_with_map_and_bvalues(pi, vals[i], sweep_genetic_distance);  // Now calling ln_likelihood_with_map_and_bvalues and sweep_genetic_distance
      //      printf("%i %f %e like=%f sweep_width=%i\n", i, sweep_genetic_distance, vals[i], likes[i], sweep_width); // Changed to sweep_genetic_distance
      sw[i] = sweep_width;
      if (i==0 || likes[i] > likes[max])
	max = i;
    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals[max-1];
    if (max==gridsize-1) maxalpha += (vals[max-1]-vals[max-2]);
    else maxalpha = vals[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }
  like = likes[max];
  sweep_width = sw[max];
  *alpha = vals[max];
  free(vals);
  free(likes);
  return like;
}

double findsweeps(char *outfn, double *p, int gridsize, int noisy, int msrep) {
  double alpha, lr, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep;
  int i, rep;
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;

  for (i=0; i<datasize; i++) 
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p));
  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  if (noisy) {
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1) {
    outfile = my_fopen(outfn, "w");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  for (rep=0; rep<gridsize; rep++) {
    if (gridsize==1) { //gridsize=1 is not recommended
      if(flagUserDefinedGrid) { sweep = data_grid[0]; }        
      else { sweep=(smin+smax)/2.0; }
    }
    else {
      if(flagUserDefinedGrid) { sweep = data_grid[rep]; }        
      else { sweep = smin + (smax-smin)*rep/(gridsize-1); }
    }
    lr = findmax_alpha(&alpha, p, sweep);
    if (rep==0 || lr > maxlr) {
      maxlr = lr;
      maxalpha = alpha;
      maxsweep = sweep;
    }
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\n", sweep, lr, alpha);
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\n", sweep, lr, alpha);
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha);
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e", maxlr, maxsweep, maxalpha);
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);
  return maxlr;
}

// Populate the grid of B-values
void PopulateBvalueGrid(void)
{
	int i;
	
	bvalue_grid = malloc(gridsize_bvalue*sizeof(double));
	
	for(i = 0; i < gridsize_bvalue; i++) {
		bvalue_grid[i] = ((double)i) / (gridsize_bvalue - 1);
	}
}

// Find sweeps using a b-value map
double findsweeps_with_bvalues(char *outfn, double *p, int gridsize, int noisy, int msrep) {
  double alpha, lr, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep;
  int i, rep;
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;
  double *p_modified = malloc(xmax*sizeof(double)); // make new allele frequency spectrum (added this)

  PopulateBvalueGrid(); // Populate the grid of B-values (added this)

  // Initialize to modified frequency spectrum (added this)
  for(i = 0; i < xmax; i++) {
    p_modified[i] = 0.0;
  }

  for (i=0; i<datasize; i++) {
    ModifySpectrumWithBvalues(p, data[i].n, data_bvalue[i].bvalue, p_modified); // Modify frequency spectrum locally for each B value (added this)
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p_modified)); // Calculate base likelihood using modified spectrum instead (changed to p_modified)
  }

  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  if (noisy) {
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1) {
    outfile = my_fopen(outfn, "w");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  for (rep=0; rep<gridsize; rep++) {
    if (gridsize==1) { //gridsize=1 is not recommended
      if(flagUserDefinedGrid) { sweep = data_grid[0]; }        
      else { sweep=(smin+smax)/2.0; }
    }
    else {
      if(flagUserDefinedGrid) { sweep = data_grid[rep]; }        
      else { sweep = smin + (smax-smin)*rep/(gridsize-1); }
    }
    lr = findmax_alpha_with_bvalues(&alpha, p, sweep); // Change to findmax_alpha_with_bvalues            [[[SHOULD CHANGE p to p_mofidied??????]]]
    if (rep==0 || lr > maxlr) {
      maxlr = lr;
      maxalpha = alpha;
      maxsweep = sweep;
    }
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\n", sweep, lr, alpha);
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\n", sweep, lr, alpha);
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha);
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e", maxlr, maxsweep, maxalpha);
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);
 
  free(p_modified); // Free temporary Bvalue-corrected frequency spectrum (added this)

  return maxlr;
}

// Find sweeps with recombination map
double findsweeps_with_map(char *outfn, double *p, int gridsize, int noisy, int msrep) {
  double alpha, lr, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep, sweep_genetic_distance=-1.0; // Added sweep_genetic_distance
  int i, rep, lower_sweep_index=0; // Added lower_sweep_index=0
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;

  for (i=0; i<datasize; i++) 
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p));
  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  if (noisy) {
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1) {
    outfile = my_fopen(outfn, "w");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  for (rep=0; rep<gridsize; rep++) {
    if (gridsize==1) { //gridsize=1 is not recommended
       if(flagUserDefinedGrid) { sweep = data_grid[0]; }        
       else { sweep=(smin+smax)/2.0; }  
     
       // Cycle through the dataset to find the pair of markers that sweep lies between (added this whole loop)
       for(i = lower_sweep_index + 1; i < datasize; i++) {
         // Found markers flanking sweep position
         // Obtain recombination distance for sweep
         if(data[i - 1].loc <= sweep && sweep <= data[i].loc) {
           lower_sweep_index = i - 1;
           sweep_genetic_distance = data_rec[i - 1].rate + ((sweep - data[i - 1].loc) / (data[i].loc - data[i - 1].loc)) * (data_rec[i].rate - data_rec[i - 1].rate);
           break;
         }
       } 
    }
    else {
       if(flagUserDefinedGrid) { sweep = data_grid[rep]; }        
       else { sweep = smin + (smax-smin)*rep/(gridsize-1); }
       
       // Set genetic distance of sweep position
       if(sweep == smin) {
         sweep_genetic_distance = data_rec[0].rate;
       }
       else if(sweep == smax) {
         sweep_genetic_distance = data_rec[datasize - 1].rate;
       }
       else {
         // Cycle through the dataset to find the pair of markers that sweep lies between
         for(i = lower_sweep_index + 1; i < datasize; i++) {
           // Found markers flanking sweep position
           // Obtain recombination distance for sweep
           if(data[i - 1].loc <= sweep && sweep <= data[i].loc) {
             lower_sweep_index = i - 1;
             sweep_genetic_distance = data_rec[i - 1].rate + ((sweep - data[i - 1].loc) / (data[i].loc - data[i - 1].loc)) * (data_rec[i].rate - data_rec[i - 1].rate);
             break;
           }
         }
       }
    }
    lr = findmax_alpha_with_map(&alpha, p, sweep_genetic_distance); // Changed to findmax_alpha_with_map_and_bvalues, call using genetic map position for the sweep instead
    if (rep==0 || lr > maxlr) {
      maxlr = lr;
      maxalpha = alpha;
      maxsweep = sweep;
    }
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\n", sweep, lr, alpha*1e6); // multiplied alpha by 1e6
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\n", sweep, lr, alpha*1e6); // Added sweep_genetic_distance and multiplied alpha by 1e6
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha*1e6); // multiplied alpha by 1e6
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e", maxlr, maxsweep, maxalpha*1e6); // multiplied alpha by 1e6
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);
  return maxlr;
}


// Find sweeps with recombination map and B-value map
double findsweeps_with_map_and_bvalues(char *outfn, double *p, int gridsize, int noisy, int msrep) {
  double alpha, lr, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep, sweep_genetic_distance=-1.0; // Added sweep_genetic_distance
  int i, rep, lower_sweep_index=0; // Added lower_sweep_index=0
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;
  double *p_modified = malloc(xmax*sizeof(double)); // make new allele frequency spectrum (added this)

  PopulateBvalueGrid(); // Populate the grid of B-values (added this)

  // Initialize to modified frequency spectrum (added this)
  for(i = 0; i < xmax; i++) {
    p_modified[i] = 0.0;
  }

  for (i=0; i<datasize; i++) {
    ModifySpectrumWithBvalues(p, data[i].n, data_bvalue[i].bvalue, p_modified); // Modify frequency spectrum locally for each B value (added this)
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p_modified)); // Calculate base likelihood using modified spectrum instead (changed to p_modified)
  }
 
  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  if (noisy) {
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1) {
    outfile = my_fopen(outfn, "w");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  for (rep=0; rep<gridsize; rep++) {
    if (gridsize==1) { //gridsize=1 is not recommended
       if(flagUserDefinedGrid) { sweep = data_grid[0]; }        
       else { sweep=(smin+smax)/2.0; }  
     
       // Cycle through the dataset to find the pair of markers that sweep lies between (added this whole loop)
       for(i = lower_sweep_index + 1; i < datasize; i++) {
         // Found markers flanking sweep position
         // Obtain recombination distance for sweep
         if(data[i - 1].loc <= sweep && sweep <= data[i].loc) {
           lower_sweep_index = i - 1;
           sweep_genetic_distance = data_rec[i - 1].rate + ((sweep - data[i - 1].loc) / (data[i].loc - data[i - 1].loc)) * (data_rec[i].rate - data_rec[i - 1].rate);
           break;
         }
       } 
    }
    else {
       if(flagUserDefinedGrid) { sweep = data_grid[rep]; }        
       else { sweep = smin + (smax-smin)*rep/(gridsize-1); }
       
       // Set genetic distance of sweep position
       if(sweep == smin) {
         sweep_genetic_distance = data_rec[0].rate;
       }
       else if(sweep == smax) {
         sweep_genetic_distance = data_rec[datasize - 1].rate;
       }
       else {
         // Cycle through the dataset to find the pair of markers that sweep lies between
         for(i = lower_sweep_index + 1; i < datasize; i++) {
           // Found markers flanking sweep position
           // Obtain recombination distance for sweep
           if(data[i - 1].loc <= sweep && sweep <= data[i].loc) {
             lower_sweep_index = i - 1;
             sweep_genetic_distance = data_rec[i - 1].rate + ((sweep - data[i - 1].loc) / (data[i].loc - data[i - 1].loc)) * (data_rec[i].rate - data_rec[i - 1].rate);
             break;
           }
         }
       }
    }
    lr = findmax_alpha_with_map_and_bvalues(&alpha, p, sweep_genetic_distance); // Changed to findmax_alpha_with_map_and_bvalues, call using genetic map position for the sweep instead, [[[SHOULD CHANGE p to p_mofidied??????]]]
    if (rep==0 || lr > maxlr) {
      maxlr = lr;
      maxalpha = alpha;
      maxsweep = sweep;
    }
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\n", sweep, lr, alpha*1e6); // multiplied alpha by 1e6
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\n", sweep, lr, alpha*1e6); // Added sweep_genetic_distance and multiplied alpha by 1e6
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha*1e6); // multiplied alpha by 1e6
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e", maxlr, maxsweep, maxalpha*1e6); // multiplied alpha by 1e6
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);

  free(p_modified); // Free temporary Bvalue-corrected frequency spectrum (added this)

  return maxlr;
}

// Print the B-value frequency spectrum
void PrintBvalueSpectrum(char *outfn, int bindex, int gridindex)
{
	int i;
	FILE *outfile = fopen(outfn, "w");

	printf("Printing spectrum for bindex = %d and gridindex = %d\n", bindex, gridindex);
	
	if(invar < 2) {
		for(i=1; i<xmax; i++) {
			fprintf(outfile, "%d\t%le\n", i,prob_bvalue[nmax][i][bindex][0][gridindex]);
			printf("%d\t%le\n", i,prob_bvalue[nmax][i][bindex][0][gridindex]);
		}
	}
	else {
		for(i=0; i<xmax; i++) {
			fprintf(outfile, "%d\t%le\n", i,prob_bvalue[nmax][i][bindex][0][gridindex]);
			printf("%d\t%le\n", i,prob_bvalue[nmax][i][bindex][0][gridindex]);
		}
	}
	
	fclose(outfile);
}

void usage() {
  printf("usage:\n");
  printf("\tto get frequency spectrum: ./SweepFinder2 -f CombinedFreqFile SpectFile\n");
  printf("\tto find sweeps:\n");
  printf("\t\t./SweepFinder2 -s G FreqFile OutFile\n");
  printf("\t\t./SweepFinder2 -sg g FreqFile OutFile\n");
  printf("\t\t./SweepFinder2 -su GridFile FreqFile OutFile\n");
//  printf("\tto analyze ms output: ./SweepFinder -m GRIDSIZE msfilename outfilename\n");
  printf("\tto find sweeps using a pre-computed frequency spectrum:\n");
  printf("\t\t./SweepFinder2 -l G FreqFile SpectFile OutFile\n");
  printf("\t\t./SweepFinder2 -lg g FreqFile SpectFile OutFile\n");
  printf("\t\t./SweepFinder2 -lu GridFile FreqFile SpectFile OutFile\n");
  printf("\tto find sweeps using pre-computed frequency spectra given recombination map:\n");
  printf("\t\t./SweepFinder2 -lr G FreqFile SpectFile RecFile OutFile\n");
  printf("\t\t./SweepFinder2 -lrg g FreqFile SpectFile RecFile OutFile\n");
  printf("\t\t./SweepFinder2 -lru GridFile FreqFile SpectFile RecFile OutFile\n");
  printf("\tto find sweeps using pre-computed frequency spectra given B values:\n");
  printf("\t\t./SweepFinder2 -lb G FreqFile SpectFile BValFile N1 N2 T OutFile\n");
  printf("\t\t./SweepFinder2 -lbg g FreqFile SpectFile BValFile N1 N2 T OutFile\n");
  printf("\t\t./SweepFinder2 -lbu GridFile FreqFile SpectFile BValFile N1 N2 T OutFile\n");
  printf("\tto find sweeps using pre-computed frequency spectra given recombination map and B values:\n");
  printf("\t\t./SweepFinder2 -lrb G FreqFile SpectFile RecFile BValFile N1 N2 T OutFile\n");
  printf("\t\t./SweepFinder2 -lrbg g FreqFile SpectFile RecFile BValFile N1 N2 T OutFile\n");
  printf("\t\t./SweepFinder2 -lrbu GridFile FreqFile SpectFile RecFile BValFile N1 N2 T OutFile\n");
  exit(-1);
}

int main(int argc, char *argv[]) {
  double *p;
  double space_between_grid_points=-1.0;
  int gridsize=-1, isnoisy=1;
  char snpfn[1000], outfn[1000], freqfn[1000], recmapfn[1000], bvaluefn[1000], gridfn[1000];
  if (argc < 3) usage();
//  if (strlen(argv[1])!=2 || argv[1][0]!='-')
 //   usage();
  if (strcmp(argv[1], "-f") == 0) {
    if (argc!=4) usage();

    sprintf(snpfn, "%s", argv[2]);
    sprintf(outfn, "%s", argv[3]);
    readsnps(snpfn);
    getfreq(outfn);
  }
  else if ( strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "-sg") == 0 || strcmp(argv[1], "-su") == 0 ) {
    if (argc!=5) usage();

    if(strcmp(argv[1], "-s") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-sg") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-su") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", "tempfreq.txt");
    sprintf(outfn, "%s", argv[4]);
    readsnps(snpfn);

    if(strcmp(argv[1], "-sg") == 0) {
	gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }

    getfreq(freqfn);
    p = loadfreq(freqfn);
    findsweeps(outfn, p, gridsize, isnoisy, -1);
    free(p);
  }
  else if (strcmp(argv[1], "-l") == 0 || strcmp(argv[1], "-lg") == 0  || strcmp(argv[1], "-lu") == 0) {
    if (argc!=6) usage();

    if(strcmp(argv[1], "-l") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-lg") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-lu") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(outfn, "%s", argv[5]);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-lg") == 0) {
	gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    p = loadfreq(freqfn);
    findsweeps(outfn, p, gridsize, isnoisy, -1);
    free(p);
  }
 /* else if (strcmp(argv[1], "-m") == 0) {
    int rep=-1;
    if (argc!=5) usage();
    gridsize = atoi(argv[2]);
    if (gridsize <=0) {
      printf("gridsize should be > 0!\n"); usage();
    }
    sprintf(snpfn, argv[3]);
    sprintf(freqfn, "tempfreq.txt");
    sprintf(outfn, argv[4]);
    rep=0;
    while (readsnps_ms(snpfn)) {
      getfreq(freqfn);
      p = loadfreq(freqfn);
      findsweeps(outfn, p, gridsize, 0, ++rep);
      free(p);
    }
  }*/
  else if(strcmp(argv[1], "-lr") == 0 || strcmp(argv[1], "-lrg") == 0  || strcmp(argv[1], "-lru") == 0) {
    if(argc!=7) {
      usage();
    }
    else {
      printf("You have chosen to get sweeps using pre-computed frequency spectra give recombination map\n");
    }

    if(strcmp(argv[1], "-lr") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-lrg") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-lru") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(recmapfn, "%s", argv[5]);
    sprintf(outfn, "%s", argv[6]);
    readsnps(snpfn);

    if(strcmp(argv[1], "-lrg") == 0) {
	gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    readrecmap(recmapfn);
    p = loadfreq(freqfn);
    findsweeps_with_map(outfn, p, gridsize, isnoisy, -1);
    free(p);
  }
  else if(strcmp(argv[1], "-lb") == 0 || strcmp(argv[1], "-lbg") == 0  || strcmp(argv[1], "-lbu") == 0) {
    if(argc!=10) {
      usage();
    }
    else {
      printf("You have chosen to get sweeps using pre-computed frequency spectra give B values\n");
    }

    if(strcmp(argv[1], "-lb") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-lbg") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-lbu") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(bvaluefn, "%s", argv[5]);
    
    N_anc = atof(argv[6]);
    N_curr = atof(argv[7]);
    split_time = atof(argv[8]);

    if(N_anc <= 0.0) {
	printf("Dipiod ancestral population should be greater than 0!\n");
	usage();
    }
	
    if(N_curr <= 0.0) {
	printf("Diploid current population should be greater than 0!\n");
	usage();
    }

    if(split_time <= 0.0) {
	printf("Divergence time in generations greater than 0!\n");
	usage();
    }

    sprintf(outfn, "%s", argv[9]);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-lbg") == 0) {
	gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    readbvalues(bvaluefn);
    p = loadfreq(freqfn);
    findsweeps_with_bvalues(outfn, p, gridsize, isnoisy, -1);
    free(p);
  }
  else if(strcmp(argv[1], "-lrb") == 0 || strcmp(argv[1], "-lrbg") == 0  || strcmp(argv[1], "-lrbu") == 0) {
    if(argc!=11) {
      usage();
    }
    else {
      printf("You have chosen to get sweeps using pre-computed frequency spectra give recombination map and B values\n");
    }

    if(strcmp(argv[1], "-lrb") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-lrbg") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-lrbu") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(recmapfn, "%s", argv[5]);
    sprintf(bvaluefn, "%s", argv[6]);

    N_anc = atof(argv[7]);
    N_curr = atof(argv[8]);
    split_time = atof(argv[9]);

    if(N_anc <= 0.0) {
	printf("Dipiod ancestral population should be greater than 0!\n");
	usage();
    }
	
    if(N_curr <= 0.0) {
	printf("Diploid current population should be greater than 0!\n");
	usage();
    }

    if(split_time <= 0.0) {
	printf("Divergence time in generations greater than 0!\n");
	usage();
    }


    sprintf(outfn, "%s", argv[10]);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-lrbg") == 0) {
	gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    readrecmap(recmapfn);
    readbvalues(bvaluefn);
    p = loadfreq(freqfn);
    findsweeps_with_map_and_bvalues(outfn, p, gridsize, isnoisy, -1);
    free(p);
  }
  else usage();
  return 0;
}
