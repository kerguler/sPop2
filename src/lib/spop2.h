#ifndef SPOP_H
#define SPOP_H

#include <float.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "uthash.h"

#define is_greater(s,p) (((s).age > (p).age) || ((s).age == (p).age && (s).devcycle > (p).devcycle) || ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development > (p).development))
#define is_smaller(s,p) (((s).age < (p).age) || ((s).age == (p).age && (s).devcycle < (p).devcycle) || ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development < (p).development))
#define is_equal(s,p) ((s).age == (p).age && (s).devcycle == (p).devcycle && (s).development == (p).development)
#define is_empty(sp,s) (((sp)->stochastic && (s).number.i==0) || (!((sp)->stochastic) && (s).number.d<=DPOP_EPS))

#define get_lchild(num) (((num)<<1)+1)
#define get_rchild(num) (((num)<<1)+2)
#define get_parent(num) ((num) ? ((num)-1)>>1 : 0)

/*
 * Data type for abundance:
 *  unsigned int: stochastic simulations
 *  double:       deterministic simulations
 */
typedef union {
  unsigned int i;
  double d;
} sdnum;

/* ******************************* */

typedef struct qunit_st *qunit;
typedef struct spop2_st *spop2;

struct qunit_st {
    double dev;
    sdnum size;
    UT_hash_handle hh;
};

typedef double (*pfunc)(double, double);

struct spop2_st {
    qunit devc;
    qunit devtable;
    sdnum size;
    sdnum completed;
    unsigned char stochastic;
    unsigned char pdist;
    pfunc cfun;
};

#define MODE_GAMMA_RAW    0
#define MODE_GAMMA_HASH   1
#define MODE_NBINOM_RAW   2
#define MODE_GAMMA_MATRIX 3
#define MODE_BINOM_RAW    4

#define MODE_ACCP_ERLANG  5
#define MODE_ACCP_FIXED   6
#define MODE_ACCP_PASCAL  7
#define MODE_ACCP_GAMMA   8
#define MODE_ACCP_CASWELL 9

#define PRINT_MODE(x) { \
char str[10][20] = {    \
"MODE_GAMMA_RAW",       \
"MODE_GAMMA_HASH",      \
"MODE_NBINOM_RAW",      \
"MODE_GAMMA_MATRIX",    \
"MODE_BINOM_RAW",       \
"MODE_ACCP_ERLANG",     \
"MODE_ACCP_FIXED",      \
"MODE_ACCP_PASCAL",     \
"MODE_ACCP_GAMMA",      \
"MODE_ACCP_CASWELL"     \
};                      \
printf("%s",str[(x)]);  \
}

/* ******************************* */

/*
 * Data structure for a cohort of individuals:
 *  age:         age of the cohort (number of iterations since entry)
 *  devcycle:    number of times development was completed
 *               this is a counter useful for certain life processes
 *  development: number of steps during development
 *               this is usually equal to age
 *  accumulate:  pointer to the accumulative development module
 *  number:      size of the cohort (integer or real)
 */
typedef struct individual_st {
  unsigned int age;
  unsigned int devcycle;
  unsigned int development;
  sdnum number;
} individual_data;

/*
 * Data structure for a population comprising a set of cohorts:
 *  individuals:  an array of pointers for each cohort
 *  ncat:         maximum size of the individuals array
 *  cat:          current size of the individuals array
 *  size:         size of the population (sum of cohort sizes)
 *  dead:         number of individuals dying during an iteration
 *  developed:    number of individuals completing development during an iteration
 *  devtable:     population of individuals completing development during an iteration
 *  gamma_mode:   selection of probability distributions/methods for the method of hazards
 *                MODE_GAMMA_RAW
 *                MODE_GAMMA_HASH
 *                MODE_NBINOM_RAW
 *  stochastic:   0: deterministic
 *                1: stochastic
 *  accumulative: logical indicator for an accumulative development process
 */
typedef struct population_st {
  individual_data *individuals;
  unsigned int ncat;
  unsigned int cat;
  sdnum size;
  sdnum dead;
  sdnum developed;
  void *devtable;
  unsigned char gamma_mode;
  unsigned char stochastic;
  unsigned char accumulative;
} *spop;

/* ******************************* */

gsl_rng *get_RAND_GSL(void);
void rng_setup(char*);
void rng_setup_seed(unsigned int, char*);
void rng_destroy(void);
double rng_exponential(const double);

/* ******************************* */

#define EPS 1e-14

#define n_MAX 400.0
#define n_STEP 1.0
#define n_i_MAX 400
#define row_SIZE 400
#define mean_MAX 201.0
#define mean_STEP 0.01
#define mean_i_MAX 20100
#define gamma_SIZE 8040000

#define MAX_DAYS 1000

#define gamma_matrix_sd 0.375

void set_gamma_mem(uint64_t);
void gamma_dist_destroy(void);
void gamma_dist_check(void);
double gamma_pdf(double, double, double);
char gamma_dist_hash(double, double, double, double *);
double gamma_dist_prob(double, double, double);
double gamma_dist_matrix(double, double);
void prepare_gamma_matrix(void);
double nbinom_prob(unsigned int, double, double);
double nbinom_dist_prob(double, double, unsigned int);

/* ******************************* */

void set_APPROX(double);

double fun_pois_C(double, double);
double fun_fixed_C(double, double);
double fun_daily_C(double, double);
double fun_unif_C(double, double);
double fun_cpois_C(double, double);

qunit qunit_new(double,sdnum);
void qunit_free(qunit *);

char spop2_get_cfun(char, pfunc *);
spop2 spop2_init(unsigned char, unsigned char);
char spop2_empty_devc(qunit *);
char spop2_empty(spop2);
char spop2_destroy(spop2 *);
#define spop2_add(pop,dev,size) {   \
  sdnum tmp;                        \
  if ((pop)->stochastic)            \
    tmp.i = (unsigned int)size;     \
  else                              \
    tmp.d = (double)size;           \
  spop2_sdadd((pop),(dev),tmp);     \
}
char spop2_sdadd(spop2, double, sdnum);
char spop2_sdpopadd(spop2, spop2, char, sdnum *);
void spop2_print(spop2);
void spop2_retrieve(spop2, char, double *, double *, unsigned int *);
char spop2_iterate(spop2, double, double);
char spop2_survive(spop2, double, char, sdnum *);

/* ******************************* */

void set_DPOP_EPS(double);
void set_DPOP_MAX_DAYS(unsigned int);

spop spop_init(unsigned char, unsigned char);
void spop_empty(spop);
void spop_destroy(spop*);
void spop_print(spop);
void spop_print_to_csv(spop);

void swap(spop, individual_data *, individual_data *);

#define spop_add(s,age,devcycle,development,number) {                     \
    sdnum spop_var_tmp;                                                   \
    if ((s)->stochastic)                                                  \
      spop_var_tmp.i = (int)(number);                                     \
    else                                                                  \
      spop_var_tmp.d = (double)(number);                                  \
    spop_sdadd((s),(age),(devcycle),(development),spop_var_tmp,0);        \
  }
void spop_sdadd(spop, unsigned int, unsigned int, unsigned int, sdnum, char);
void spop_popadd(spop, spop, char);

typedef double (*prob_func)(unsigned int, double, double, double);
typedef void (*iter_func)(const individual_data*, double*, double*, double*);
char spop_iterate(spop, double, double, double, iter_func, double, double, double, iter_func, unsigned char);

/* ******************************* */

#define SPOPLIB_MoH 0
#define SPOPLIB_MoA 1

typedef struct spoplib_link_st *spoplib_link;

struct spoplib_link_st {
    unsigned int id;
    char type;
    void *pop;
    UT_hash_handle hh;
};

void spoplib_print(unsigned int);
void spoplib_destroy(unsigned int);
void spoplib_destroy_all(void);
unsigned int spoplib_init(unsigned char, unsigned char);
void spoplib_add(unsigned int, unsigned int, unsigned int, unsigned int, double, double);
void spoplib_iterate(unsigned int, double, double, double, double, double, double);
void spoplib_read(unsigned int, double *, double *, double *);
void spoplib_retrieve(unsigned int, char, double *, double *, unsigned int *);

/* ******************************* */

#endif
