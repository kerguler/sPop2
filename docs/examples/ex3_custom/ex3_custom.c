#include "spop2/spop2.h"

void death(const individual_data *ind, double *death_prob, double *death_mean, double *death_sd) {
    (*death_prob) = 0;
    (*death_mean) = 480.0 - (ind->devcycle > 4 ? 240.0 : 48.0 * ind->devcycle);
    (*death_sd) = 0.1 * (*death_mean);
}

int main(void) {
    int i, j;
    spop vec;
    rng_setup("spop.c");
    //
    vec = spop_init(1, MODE_GAMMA_HASH);
    for (j = 0; j < 100; j++) {
        spop_add(vec, 0, 0, 0, 1000);
        for (i = 0; i < 720; i++) {
            spop_iterate(vec,
                         0, 50.0, 10.0, 0,
                         0, 0, 0, death,
                         0);
            spop_popadd(vec, vec->devtable);
            printf("1,%d,%d,%d\n", i, vec->size.i, vec->developed.i);
        }
        spop_empty(vec);
    }
    spop_destroy(&vec);
    //
    vec = spop_init(0, MODE_GAMMA_HASH);
    spop_add(vec, 0, 0, 0, 1000);
    for (i = 0; i < 720; i++) {
        spop_iterate(vec,
                     0, 50.0, 10.0, 0,
                     0, 0, 0, death,
                     0);
        spop_popadd(vec, vec->devtable);
        printf("0,%d,%g,%g\n", i, vec->size.d, vec->developed.d);
    }
    spop_destroy(&vec);
    //
    rng_destroy();
    return 0;
}
