#ifndef QNTA_H
#define QNTA_H

#include "spop2/uthash.h"
#include "spop2/spop2.h"

typedef struct qunit_st *qunit;
typedef struct quant_st *quant;

struct qunit_st {
    double dev;
    sdnum size;
    UT_hash_handle hh;
};

struct quant_st {
  qunit devc;
  sdnum size;
  sdnum completed;
  unsigned char stochastic;
  unsigned char pdist;
  pfunc cfun;
};

qunit qunit_new(double,sdnum);
void qunit_free(qunit *);

quant quant_init(unsigned char, unsigned char);
char quant_empty(quant);
char quant_destroy(quant *);
#define quant_add(pop,dev,size) {   \
  sdnum tmp;                        \
  if ((pop)->stochastic)            \
    tmp.i = (unsigned int)size;     \
  else                              \
    tmp.d = (double)size;           \
  quant_sdadd((pop),(dev),tmp);     \
}
char quant_sdadd(quant, double, sdnum);
void quant_print(quant);
char quant_iterate(quant, double, double);

#endif
