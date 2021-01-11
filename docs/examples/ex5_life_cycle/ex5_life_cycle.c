#include "spop2/spop2.h"

void print_out(unsigned int tm, spop2 *egg, spop2 *larva, spop2 *pupa, spop2 *adult) {
    printf("%d,%g,%g,%g,%g\n",tm,(*egg)->size.d,(*larva)->size.d,(*pupa)->size.d,(*adult)->size.d);
}

int main(void) {
    unsigned char mode = MODE_ACCP_ERLANG;
    unsigned int tm = 0;
    //
    spop2 egg = spop2_init(0, mode);
    spop2 larva = spop2_init(0, mode);
    spop2 pupa = spop2_init(0, mode);
    spop2 adult = spop2_init(0, mode);
    //
    spop2_add(egg, 0, 100);
    print_out(tm, &egg, &larva, &pupa, &adult);
    //
    for (tm=1; tm<200; tm++) {
        spop2_iterate(egg,   10, 2, 0, 0);
        spop2_iterate(larva, 10, 2, 0, 0);
        spop2_iterate(pupa,  10, 2, 0, 0);
        spop2_iterate(adult, 10, 2, 0, 0);
        //
        spop2_add(egg, 0, 1.0 * adult->size.d);
        spop2_add(larva, 0, egg->developed.d);
        spop2_add(pupa, 0, larva->developed.d);
        spop2_add(adult, 0, 0.5 * pupa->developed.d);
        //
        print_out(tm, &egg, &larva, &pupa, &adult);
    }
    //
    spop2_destroy(&egg);
    spop2_destroy(&larva);
    spop2_destroy(&pupa);
    spop2_destroy(&adult);
    return 0;
}
