#ifndef ACCP_H
#define ACCP_H

#include <float.h>
#include "spop2.h"

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

#define ERLANG 0
#define FIXED  1
#define PASCAL 2
#define GAMMA  3

char choose(unsigned int, double *, unsigned int, unsigned int *);
accp accp_init(unsigned char, unsigned char);
char accp_empty(accp);
char accp_destroy(accp *);
#define accp_add(pop,stage,size) {    \
    sdnum tmp;                        \
    if ((pop)->stochastic)            \
      tmp.i = (unsigned int)(size);   \
    else                              \
      tmp.d = (double)(size);         \
    accp_sdadd((pop),(stage),(tmp));  \
  }
char accp_sdadd(accp,unsigned int,sdnum);
void accp_print(accp);
sdnum accp_get_size(accp);
char chain_resize(accp, unsigned int);
char accp_iterate(accp, double, double);

#endif
