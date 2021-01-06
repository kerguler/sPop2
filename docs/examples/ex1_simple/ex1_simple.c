#include "spop2/spop2.h"

int main(void) {
    spop pop = spop_init(0, MODE_GAMMA_HASH);
    spop_add(pop, 0, 0, 0, 1000);
    spop_iterate(pop,
                 0, 10, 5, 0,
                 0.25, 0, 0, 0,
                 0);
    spop_print(pop);
    spop_destroy(&pop);
    return 0;
}
