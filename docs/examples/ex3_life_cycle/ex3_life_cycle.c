#include "spop2/spop2.h"

void print_out(unsigned int tm, spop *egg, spop *larva, spop *pupa, spop *adult) {
    printf("%d,%g,%g,%g,%g\n",tm,(*egg)->size.d,(*larva)->size.d,(*pupa)->size.d,(*adult)->size.d);
}

int main(void) {
    unsigned int mode = MODE_GAMMA_HASH;
    unsigned int tm = 0;
    //
    spop egg = spop_init(0, mode);
    spop larva = spop_init(0, mode);
    spop pupa = spop_init(0, mode);
    spop adult = spop_init(0, mode);
    //
    spop_add(egg, 0, 0, 0, 0, 100);
    print_out(tm, &egg, &larva, &pupa, &adult);
    //
    for (tm=1; tm<100; tm++) {
        spop_iterate(egg,  0, 10, 1, 0,  0, 0, 0, 0,  0);
        spop_iterate(larva,  0, 10, 1, 0,  0, 0, 0, 0,  0);
        spop_iterate(pupa,  0, 10, 1, 0,  0, 0, 0, 0,  0);
        spop_iterate(adult,  0, 0, 0, 0,  0, 10, 1, 0,  0);
        //
        spop_add(egg, 0, 0, 0, 0, 1.0 * adult->size.d);
        spop_add(larva, 0, 0, 0, 0, egg->developed.d);
        spop_add(pupa, 0, 0, 0, 0, larva->developed.d);
        spop_add(adult, 0, 0, 0, 0, 0.5 * pupa->developed.d);
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
