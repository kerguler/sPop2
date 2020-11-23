#ifndef RAN_GEN_H
#define RAN_GEN_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *get_RAND_GSL(void);
void rng_setup(char*);
void rng_setup_seed(unsigned int, char*);
void rng_destroy(void);
double rng_exponential(const double);

#endif
