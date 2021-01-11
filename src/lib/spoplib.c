/*
 * WARNING: Do not forget to initiate RAND_GSL from Python for stochastic simulations!
 */

#include "spop2.h"
#include "uthash.h"

unsigned int hash_counter = 0;
spoplib_link pop_link = 0;

void spoplib_destroy(unsigned int id) {
    spop pop;
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (p) {
        switch (p->type) {
            case SPOPLIB_MoH:
                pop = (spop)(p->pop);
                spop_destroy(&pop);
                break;
            case SPOPLIB_MoA:
            default:
                pop2 = (spop2)(p->pop);
                spop2_destroy(&pop2);
                break;
        }
        HASH_DEL(pop_link, p);
        free(p);
    }
}

void spoplib_destroy_all(void) {
    spop pop;
    spop2 pop2;
    spoplib_link p, tmp;
    HASH_ITER(hh, pop_link, p, tmp) {
        switch (p->type) {
            case SPOPLIB_MoH:
                pop = (spop)(p->pop);
                spop_destroy(&pop);
                break;
            case SPOPLIB_MoA:
            default:
                pop2 = (spop2)(p->pop);
                spop2_destroy(&pop2);
                break;
        }
        HASH_DEL(pop_link, p);
        free(p);
    }
    free(pop_link);
    pop_link = 0;
    hash_counter = 0;
}

void spoplib_print(unsigned int id) {
    spop pop;
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (!p) {
        printf("Error: spoplib id %d not found!\n", id);
        return;
    }
    switch (p->type) {
        case SPOPLIB_MoH:
            pop = (spop) (p->pop);
            spop_print(pop);
            break;
        case SPOPLIB_MoA:
        default:
            pop2 = (spop2) (p->pop);
            spop2_print(pop2);
            break;
    }
}

unsigned int spoplib_init(unsigned char stochastic, unsigned char gamma_mode) {
    spop pop;
    spop2 pop2;
    spoplib_link link;
    unsigned int id = ++hash_counter;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(double), p);
    if (p) {
        printf("Error: Multiple counter ids %d!\n", hash_counter);
        return 0;
    } else {
        switch (gamma_mode) {
            case MODE_GAMMA_RAW:
            case MODE_GAMMA_HASH:
            case MODE_NBINOM_RAW:
            case MODE_GAMMA_MATRIX:
            case MODE_BINOM_RAW:
                pop = spop_init(stochastic, gamma_mode);
                link = (spoplib_link) malloc(sizeof(struct spoplib_link_st));
                link->id = id;
                link->type = SPOPLIB_MoH;
                link->pop = (void *) pop;
                HASH_ADD(hh, pop_link, id, sizeof(unsigned int), link);
                break;
            case MODE_ACCP_ERLANG:
            case MODE_ACCP_FIXED:
            case MODE_ACCP_PASCAL:
            case MODE_ACCP_GAMMA:
            case MODE_ACCP_CASWELL:
                pop2 = spop2_init(stochastic, gamma_mode);
                link = (spoplib_link) malloc(sizeof(struct spoplib_link_st));
                link->id = id;
                link->type = SPOPLIB_MoA;
                link->pop = (void *) pop2;
                HASH_ADD(hh, pop_link, id, sizeof(unsigned int), link);
                break;
            default:
                printf("Error: Unsupported probability distribution %d!\n", gamma_mode);
                return 0;
                break;
        }
    }
    return id;
}

void spoplib_add(unsigned int id, unsigned int age, unsigned int devcycle, unsigned int development, double stage, double number) {
    spop pop;
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (!p) {
        printf("Error: spoplib id %d not found!\n", id);
        return;
    }
    switch (p->type) {
        case SPOPLIB_MoH:
            pop = (spop)(p->pop);
            spop_add(pop,age,devcycle,development,number);
            break;
        case SPOPLIB_MoA:
            pop2 = (spop2)(p->pop);
            spop2_add(pop2,stage,number);
            break;
        default:
            printf("Error: Unsupported spoplib type %d!\n", p->type);
            return;
            break;
    }
}

void spoplib_iterate(unsigned int id,
                     double dev_prob,
                     double dev_mean,
                     double dev_sd,
                     double death_prob,
                     double death_mean,
                     double death_sd) {
    spop pop;
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (!p) {
        printf("Error: spoplib id %d not found!\n", id);
        return;
    }
    switch (p->type) {
        case SPOPLIB_MoH:
            pop = (spop)(p->pop);
            spop_iterate(pop,
                         dev_prob,
                         dev_mean,
                         dev_sd,
                         0,
                         death_prob,
                         death_mean,
                         death_sd,
                         0,
                         0);
            break;
        case SPOPLIB_MoA:
            pop2 = (spop2)(p->pop);
            spop2_iterate(pop2,
                          dev_mean,
                          dev_sd,
                          death_prob,
                          0);
            break;
        default:
            printf("Error: Unsupported spoplib type %d!\n", p->type);
            return;
            break;
    }
}

void spoplib_read(unsigned int id, double *size, double *developed, double *dead) {
    spop pop;
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (!p) {
        printf("Error: spoplib id %d not found!\n", id);
        return;
    }
    switch (p->type) {
        case SPOPLIB_MoH:
            pop = (spop)(p->pop);
            if (pop->stochastic) {
                *size = (double)(pop->size.i);
                *developed = (double)(pop->developed.i);
                *dead = (double)(pop->dead.i);
            } else {
                *size = pop->size.d;
                *developed = pop->developed.d;
                *dead = pop->dead.d;
            }
            break;
        case SPOPLIB_MoA:
            pop2 = (spop2)(p->pop);
            if (pop2->stochastic) {
                *size = (double)(pop2->size.i);
                *developed = (double)(pop2->developed.i);
                *dead = (double)(pop2->dead.i);
            } else {
                *size = pop2->size.d;
                *developed = pop2->developed.d;
                *dead = pop2->dead.d;
            }
            break;
        default:
            printf("Error: Unsupported spoplib type %d!\n", p->type);
            return;
            break;
    }
}

void spoplib_retrieve(unsigned int id, char devtable, double *dev, double *size, unsigned int *limit) {
    spop2 pop2;
    spoplib_link p;
    HASH_FIND(hh, pop_link, &id, sizeof(unsigned int), p);
    if (!p) {
        printf("Error: spoplib id %d not found!\n", id);
        return;
    }
    if (!(p->type == SPOPLIB_MoA)) return;
    pop2 = (spop2) (p->pop);
    spop2_retrieve(pop2, devtable, dev, size, limit);
}
