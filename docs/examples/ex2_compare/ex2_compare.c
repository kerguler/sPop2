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
    spop2 pop2;
    //
    for (mode=0; mode<5; mode++) {
        switch (modes[mode]) {
            case MODE_GAMMA_HASH:
            case MODE_NBINOM_RAW:
                pop = spop_init(0, modes[mode]);
                spop_add(pop, 0, 0, 0, 1000);
                printf("%d,%d,%g\n", mode, 0, pop->size.d);
                for (tm = 1; tm < 50; tm++) {
                    spop_iterate(pop, 0, 20, 10, 0, 0, 0, 0, 0, 0);
                    printf("%d,%d,%g\n", mode, tm, pop->size.d);
                }
                spop_destroy(&pop);
                break;
            case MODE_ACCP_ERLANG:
            case MODE_ACCP_FIXED:
            case MODE_ACCP_PASCAL:
                pop2 = spop2_init(0, modes[mode]);
                spop2_add(pop2, 0, 1000);
                printf("%d,%d,%g\n", mode, 0, pop2->size.d);
                for (tm = 1; tm < 50; tm++) {
                    spop2_iterate(pop2, 20, 10);
                    printf("%d,%d,%g\n", mode, tm, pop2->size.d);
                }
                spop2_destroy(&pop2);
        }
    }
    return 0;
}
