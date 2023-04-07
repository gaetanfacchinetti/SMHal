#ifndef SPECIALINTEGRATIONS
#define SPECIALINTEGRATIONS

#include "mymath.h"

class SpecialIntegrations
{

public:
    SpecialIntegrations(){};

    // Functions for the evaluation of cross sections
    dcomp IntJt(int n, double a, double b, double tmin, double tmax);
    dcomp IntJu(int n, double r, double b, double tmin, double tmax);
    dcomp IntJtComp(int n, int m, double a, double b, double tmin, double tmax);
    dcomp IntJuComp(int n, int m, double r, double b, double tmin, double tmax);

    dcomp IntJttComp(int n, int p, int q, double a, double b, double c, double d, double tmin, double tmax);
    dcomp IntJtuComp(int n, int p, int q, double a, double b, double rc, double d, double tmin, double tmax);
    dcomp IntJuuComp(int n, int p, int q, double ra, double b, double rc, double d, double tmin, double tmax);

    double IntK(int n, double tmin, double tmax);


    // Functions for the evaluation of the s-wave term of thermally averaged cross sections
    dcomp IntHtComp(int n, int m, double a, double b, double t, double pkl);
    dcomp IntHuComp(int n, int m, double a, double b, double t, double pkl);
    dcomp IntHttComp(int n, int p, int q, double a, double b, double c, double d, double t, double pkl);
    dcomp IntHtuComp(int n, int p, int q, double a, double b, double r, double d, double t, double pkl);
    dcomp IntHuuComp(int n, int p, int q, double r1, double b, double r2, double d, double t, double pkl);

    // Other functions
    int Binom(int k, int n);  // Fast computation of the usual binomial coefficient
    int BinomGen(int k, int n);  // General definition of the binomial coefficient
    dcomp FF(int n, int m, int p, dcomp x);
    int GG(int n, int m, int p, int q);
    dcomp IntJtCompTest(int n, int m, double a, double b, double tmin, double tmax);
};

#endif