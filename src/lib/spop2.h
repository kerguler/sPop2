#ifndef SPOP_H
#define SPOP_H

#include <float.h>

// double DPOP_EPS;
// double DPOP_MAX_DAYS;

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

/*
 * Data structure for a cohort of individuals:
 *  age:         age of the cohort (number of iterations since entry)
 *  devcycle:    number of times development was completed
 *               this is a counter useful for certain life processes
 *  development: number of steps during development
 *               this is usually equal to age
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
 *  individuals: an array of pointers for each cohort
 *  ncat:        maximum size of the individuals array
 *  cat:         current size of the individuals array
 *  size:        size of the population (sum of cohort sizes)
 *  dead:        number of individuals dying during an iteration
 *  developed:   number of individuals completing development during an iteration
 *  devtable:    population of individuals completing development during an iteration
 *  gamma_mode:  selection of probability distributions/methods for the method of hazards
 *               MODE_GAMMA_RAW
 *               MODE_GAMMA_HASH
 *               MODE_NBINOM_RAW
 *  stochastic:  0: deterministic
 *               1: stochastic
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
} *spop;

void set_DPOP_EPS(double);
void set_DPOP_MAX_DAYS(unsigned int);

spop spop_init(unsigned char, unsigned char);
void spop_empty(spop);
void spop_destroy(spop*);
void spop_print(spop);
void spop_print_to_csv(spop);

void swap(spop, individual_data *, individual_data *);

#define spop_add(s,age,devcycle,development,number) {   \
    sdnum tmp;                                          \
    if ((s)->stochastic)                                \
      tmp.i = (int)(number);                            \
    else                                                \
      tmp.d = (double)(number);                            \
    spop_sdadd((s),(age),(devcycle),(development),tmp);    \
  }
void spop_sdadd(spop, unsigned int, unsigned int, unsigned int, sdnum);

void spop_popadd(spop, spop);

typedef double (*prob_func)(unsigned int, double, double, double);
typedef void (*iter_func)(const individual_data*, double*, double*, double*);
char spop_iterate(spop, double, double, double, iter_func, double, double, double, iter_func, unsigned char);

#endif
