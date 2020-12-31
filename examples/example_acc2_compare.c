/**
 * This is an example file for spop2
 * ----------------------------------------------
 * Compile with:
 * $ gcc -Wall -lm -lspop2 -lgsl -o example_acc2_compare example_acc2_compare.c
 * Execute:
 * ./example_acc2_compare
 */

#include "spop2/spop2.h"

int main(void) {
    spop pop = spop_init(0, MODE_GAMMA_HASH);
    spop_add(pop, 0, 0, 0, 0, 1000);
    spop_iterate(pop,
                 0, 10, 5, 0,
                 0.25, 0, 0, 0,
                 0);
    spop_print(pop);
    spop_destroy(&pop);

    pop = spop_init(0, MODE_ACCP_ERLANG);
    spop_add(pop, 0, 0, 0, 0, 1000);
    spop_iterate(pop,
                 0, 10, 5, 0,
                 0.25, 0, 0, 0,
                 0);
    spop_print(pop);
    spop_destroy(&pop);
    return 0;
}
