//
// Created by kamil on 18/11/2020.
//

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../lib/ran_gen.h"
#include "../lib/gamma.h"
#include "../lib/spop2.h"

extern gsl_rng *RAND_GSL;

#define MAXR 1000
// #define MAXX 100
#define MAXX 10000
#define RB   0.05

// Deterministic exponential decay
void survD() {
    double X0 = MAXX;
    double rB = RB;
    double tm = 0;
    for (tm=0; tm<100; tm+=0.01) {
        printf("0,%g,%g\n",tm,X0*exp(-rB*tm));
    }
}

// Discrete-time deterministic mortality model
void survL() {
    double X0 = MAXX;
    double rB = RB;
    double tm = 0;
    for (tm=0; tm<100; tm+=1.0) {
        printf("1,%g,%g\n",tm,X0*pow(1.0-rB,tm));
    }
}

// Stochastic master equation for spontaneous decay
void survM() {
    double dt = 0.01;
    double hB = RB;
    int X = MAXX;
    double tm = 0;
    double tt = 0;
    double to = 0;
    double a = 0;
    int r = 0;
    for (r=0; r<MAXR; r++) {
        X = MAXX;
        tm = 0;
        to = 0;
        while (tm < 100) {
            a = hB * X;
            tt = tm + gsl_ran_exponential(RAND_GSL,1.0/a);
            while (tt > to) {
                printf("2,%g,%d\n",to,X);
                to += dt;
                if (to > 100) goto end;
            }
            X--;
            tm = tt;
        }
        end:;
    }
}

// Discrete-time stochastic mortality model
void survB() {
    double pB = (1.0-RB);
    int X = MAXX;
    int r = 0;
    double tm = 0;
    for (r=0; r<MAXR; r++) {
        X = MAXX;
        for (tm=0; tm<100; tm+=1.0) {
            printf("3,%g,%d\n",tm,X);
            X = gsl_ran_binomial(RAND_GSL,pB,X);
        }
    }
}

int main(void) {
    rng_setup("Survival");

    survD();
    survL();
    survM();
    survB();

    rng_destroy();
    return 0;
}
