#include "spop2/spop2.h"

int main(void) {
    unsigned char modes[5] = {
        MODE_GAMMA_HASH,
        MODE_NBINOM_RAW,
        MODE_ACCP_ERLANG,
        MODE_ACCP_FIXED,
        MODE_ACCP_PASCAL
    };
    unsigned char mode;
    unsigned int tm;
    spop pop;
    //
    for (mode=0; mode<5; mode++) {
        pop = spop_init(0, modes[mode]);
        //
        spop_add(pop, 0, 0, 0, 0, 1000);
        printf("%d,%d,%g\n", mode, 0, pop->size.d);
        //
        for (tm=1; tm<50; tm++) {
            spop_iterate(pop,  0, 20, 10, 0,  0, 0, 0, 0,  0);
            printf("%d,%d,%g\n", mode, tm, pop->size.d);
        }
        //
        spop_destroy(&pop);
    }
    return 0;
}
