#include "spop2/spop2.h"

void print_out(unsigned int tm, double S, spop *E, spop *I, double R) {
    printf("%d,%g,%g,%g,%g\n",tm,S,(*E)->size.d,(*I)->size.d,R);
}

int main(void) {
    unsigned char mode = MODE_GAMMA_HASH;
    unsigned int tm = 0;
    //
    double S = 0;
    spop E = spop_init(0, mode);
    spop I = spop_init(0, mode);
    double R = 0;
    //
    double v = 0;
    //
    S = 99;
    spop_add(E, 0, 0, 0, 1);
    print_out(tm, S, &E, &I, R);
    //
    for (tm=1; tm<100; tm++) {
        spop_iterate(E,  0, 1, 1, 0,  0, 0, 0, 0,  0);
        spop_iterate(I,  0, 7, 1, 0,  0, 0, 0, 0,  0);
        //
        v = 0.4 * I->size.d / (S + E->size.d + I->size.d + R);
        //
        R += I->developed.d;
        spop_add(I, 0, 0, 0, E->developed.d);
        spop_add(E, 0, 0, 0, v * S);
        S -= v * S;
        //
        print_out(tm, S, &E, &I, R);
    }
    //
    spop_destroy(&E);
    spop_destroy(&I);
    return 0;
}
