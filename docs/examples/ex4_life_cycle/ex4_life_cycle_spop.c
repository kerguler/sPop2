#include "spop2/spop2.h"

void print_out(unsigned int tm, spop *egg, spop *larva, spop *pupa, spop *adult) {
    printf("%d,%g,%g,%g,%g\n",tm,(*egg)->size.d,(*larva)->size.d,(*pupa)->size.d,(*adult)->size.d);
}

int main(void) {
    unsigned char mode = MODE_ACCP_ERLANG;
    unsigned int tm = 0;
    //
    set_APPROX(1e-2);
    //
    spop egg = spop_init(0, mode);
    spop larva = spop_init(0, mode);
    spop pupa = spop_init(0, mode);
    spop adult = spop_init(0, mode);
    //
    spop_add(egg, 0, 0, 0, 0, 100);
    print_out(tm, &egg, &larva, &pupa, &adult);
    //
    for (tm=1; tm<25; tm++) {
        spop_iterate(egg,  0, 10, 1, 0,  0, 0, 0, 0,  1);
        spop_iterate(larva,  0, 10, 1, 0,  0, 0, 0, 0,  1);
        spop_iterate(pupa,  0, 10, 1, 0,  0, 0, 0, 0,  1);
        spop_iterate(adult,  0, 10, 1, 0,  0, 0, 0, 0,  1);
        //
        // spop_add(egg, 0, 0, 0, 0, 1.0 * adult->size.d);
        spop_popadd(larva, egg, 1);
        spop_popadd(pupa, larva, 1);
        if (tm==20) {
            printf("\n");
            spop_print(pupa);
        }
        spop_iterate(pupa,  0.5, 0, 0, 0,  0, 0, 0, 0,  1);
        if (tm==20) {
            spop_print(pupa);
            printf("\n");
            spop_print(adult);
        }
        spop_popadd(adult, pupa->devtable, 0);
        if (tm==20) {
            spop_print(adult);
            printf("\n");
        }
        //
        print_out(tm, &egg, &larva, &pupa, &adult);
    }
    //
    spop_destroy(&egg);
    spop_destroy(&larva);
    spop_destroy(&pupa);
    spop_destroy(&adult);
    return 0;
}
