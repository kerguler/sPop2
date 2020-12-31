/**
 * This is an example file for spop2
 * ----------------------------------------------
 * Compile with:
 * $ gcc -I/usr/local -Wall -lm -lspop2 -lgsl -o example_acc1_simple example_acc1_simple.c
 * Execute:
 * ./example_acc1_simple
 */

#include "spop2/spop2.h"

int main(void) {
    spop pop = spop_init(0, MODE_ACCP_ERLANG);
    spop_add(pop, 0, 0, 0, 0, 1000);
    spop_iterate(pop,
                 0, 10, 5, 0,
                 0.25, 0, 0, 0,
                 0);
    spop_print(pop);
    spop_destroy(&pop);
    return 0;
}
