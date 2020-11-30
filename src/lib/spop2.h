#ifndef SPOP_H
#define SPOP_H

#include <float.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

typedef struct stage_st *chain;
typedef struct accumulative_st *accp;

struct stage_st {
    sdnum size;
    chain prev;
    chain next;
};

typedef double (*pfunc)(double, double);

struct accumulative_st {
    chain first;
    chain last;
    unsigned int cat;
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
  accp accumulate;
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

char choose(unsigned int, double *, unsigned int, unsigned int *);
accp accp_init(unsigned char, unsigned char);
char accp_empty(accp);
char accp_destroy(accp *);
#define accp_add(pop,stage,size) {             \
    sdnum accp_var_tmp;                        \
    if ((pop)->stochastic)                     \
      accp_var_tmp.i = (unsigned int)(size);   \
    else                                       \
      accp_var_tmp.d = (double)(size);         \
    accp_sdadd((pop),(stage),accp_var_tmp);    \
  }
char accp_sdadd(accp, unsigned int, sdnum);
void accp_print(accp);
sdnum accp_get_size(accp);
char chain_resize(accp, unsigned int);
char accp_iterate(accp, double, double);
char accp_survive(accp, double, sdnum *);

/* ******************************* */

void set_DPOP_EPS(double);
void set_DPOP_MAX_DAYS(unsigned int);

spop spop_init(unsigned char, unsigned char, unsigned char);
void spop_empty(spop);
void spop_destroy(spop*);
void spop_print(spop);
void spop_print_to_csv(spop);

void swap(spop, individual_data *, individual_data *);

#define spop_add(s,age,devcycle,development,stage,number) {               \
    sdnum spop_var_tmp;                                                   \
    if ((s)->stochastic)                                                  \
      spop_var_tmp.i = (int)(number);                                     \
    else                                                                  \
      spop_var_tmp.d = (double)(number);                                  \
    accp dev = 0;                                                         \
    if ((s)->accumulative) {                                              \
       dev = accp_init((s)->stochastic,(s)->gamma_mode);                  \
       accp_sdadd(dev,(stage),spop_var_tmp);                                  \
    }                                                                     \
    spop_sdadd((s),(age),(devcycle),(development),dev,spop_var_tmp);      \
  }
void spop_sdadd(spop, unsigned int, unsigned int, unsigned int, accp, sdnum);

void spop_popadd(spop, spop);

typedef double (*prob_func)(unsigned int, double, double, double);
typedef void (*iter_func)(const individual_data*, double*, double*, double*);
char spop_iterate(spop, double, double, double, iter_func, double, double, double, iter_func, unsigned char);

/* ******************************* */

#endif
