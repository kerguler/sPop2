#include "spop2/spop2.h"

void run_MoH() {
    unsigned char modes[2] = {
            MODE_GAMMA_HASH,
            MODE_NBINOM_RAW
    };
    unsigned char mode;
    unsigned int tm;
    spop pop;
    //
    for (mode = 0; mode < 2; mode++) {
        pop = spop_init(0, modes[mode]);
        spop_add(pop, 0, 0, 0, 100);
        printf("%d,%d,%g,%g,%g\n", mode, 0, pop->size.d, pop->dead.d, pop->developed.d);
        for (tm = 1; tm < 50; tm++) {
            spop_iterate(pop,  0, 20, 10, 0,  0.1, 0, 0, 0,  0);
            printf("%d,%d,%g,%g,%g\n", mode, tm, pop->size.d, pop->dead.d, pop->developed.d);
        }
        spop_destroy(&pop);
    }
}

void run_MoA() {
    unsigned char modes[3] = {
            MODE_ACCP_ERLANG,
            MODE_ACCP_FIXED,
            MODE_ACCP_PASCAL
    };
    unsigned char mode;
    unsigned int tm;
    spop2 pop2;
    //
    for (mode = 0; mode < 3; mode++) {
        pop2 = spop2_init(0, modes[mode]);
        spop2_add(pop2, 0, 100);
        printf("%d,%d,%g,%g,%g\n", mode, 0, pop2->size.d, pop2->dead.d, pop2->developed.d);
        for (tm = 1; tm < 50; tm++) {
            spop2_iterate(pop2, 20, 10, 0.1, 0);
            printf("%d,%d,%g,%g,%g\n", mode, tm, pop2->size.d, pop2->dead.d, pop2->developed.d);
        }
        spop2_destroy(&pop2);
    }
}

int main(void) {
    run_MoH();
    run_MoA();
    return 0;
}
