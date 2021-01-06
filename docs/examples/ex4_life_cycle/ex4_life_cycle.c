#include "spop2/spop2.h"

void print_out(unsigned int tm, quant *egg, quant *larva, quant *pupa, quant *adult) {
    printf("%d,%g,%g,%g,%g\n",tm,(*egg)->size.d,(*larva)->size.d,(*pupa)->size.d,(*adult)->size.d);
}

int main(void) {
    unsigned char mode = MODE_ACCP_ERLANG;
    unsigned int tm = 0;
    sdnum tmp;
    //
    // set_APPROX(1e-2);
    //
    quant egg = quant_init(0, mode);
    quant larva = quant_init(0, mode);
    quant pupa = quant_init(0, mode);
    quant adult = quant_init(0, mode);
    //
    quant_add(egg, 0, 100);
    print_out(tm, &egg, &larva, &pupa, &adult);
    //
    for (tm=1; tm<100; tm++) {
        quant_iterate(egg,   10, 1);
        quant_iterate(larva, 10, 1);
        quant_iterate(pupa,  10, 1);
        quant_iterate(adult, 10, 1);
        //
        quant_add(egg, 0, 1.0 * adult->size.d);
        quant_sdpopadd(larva, egg, 1, &tmp);
        quant_sdpopadd(pupa, larva, 1, &tmp);
        if (tm==20) {
            //printf("\n");
            //quant_print(pupa);
        }
        quant_survive(pupa,  0.5, 1, &tmp);
        if (tm==20) {
            //quant_print(pupa);
            //printf("\n");
            //quant_print(adult);
        }
        quant_sdpopadd(adult, pupa, 1, &tmp);
        if (tm==20) {
            //quant_print(adult);
            //printf("\n");
        }
        //
        print_out(tm, &egg, &larva, &pupa, &adult);
    }
    //
    quant_destroy(&egg);
    quant_destroy(&larva);
    quant_destroy(&pupa);
    quant_destroy(&adult);
    return 0;
}
