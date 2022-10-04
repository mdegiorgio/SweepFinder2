#ifndef __MY_RND_H__
#define __MY_RND_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gammLn(double xx);
void SetSeed3(int seed[3]);
void SetSeed(int seed);
void saveseed();
void loadseed();
double uniform();
int rand_int(int N);
int Binomial(int n, double p);
int geometric(int expected_value);
double expo(double c);
double Rgamma(double a, double lambda);
void RDirichlet(const double * a, const int k, double * b);
double RBeta(const double alpha, const double beta);
double Poisson(double xm);
void seed(char *cmd, int rank, int size);
#endif
