/*
 * WARNING: Do not forget to initiate RAND_GSL from Python for stochastic simulations!
 */

unsigned int MAXPOPS = 1;
spop *pop_list = 0;
spop2 *pop2_list = 0;

unsigned int spoplib_init(unsigned char stochastic, unsigned char gamma_mode) {
    if (pop_list)
        spoplib_destroy_all();
    pop_list = (spop *)calloc(MAXPOPS,sizeof(spop));
    pop_list[0] = spop_init(stochastic,gamma_mode);
    return 0;
}

void spoplib_add(unsigned int id, unsigned int age, unsigned int devcycle, unsigned int development, double stage, double number) {
    if (pop_list[id])
        spop_add(pop_list[id],age,devcycle,development,stage,number);
}

void spoplib_iterate(unsigned int id,
                     double dev_prob,
                     double dev_mean,
                     double dev_sd,
                     double death_prob,
                     double death_mean,
                     double death_sd) {
    if (pop_list[id])
        spop_iterate(pop_list[id],
                     dev_prob,
                     dev_mean,
                     dev_sd,
                     0,
                     death_prob,
                     death_mean,
                     death_sd,
                     0,
                     0);
}

void spoplib_read(unsigned int id, double *size, double *developed, double *dead) {
    if (pop_list[id]) {
        if (pop_list[id]->stochastic) {
            *size = (double)(pop_list[id]->size.i);
            *developed = (double)(pop_list[id]->developed.i);
            *dead = (double)(pop_list[id]->dead.i);
        } else {
            *size = pop_list[id]->size.d;
            *developed = pop_list[id]->developed.d;
            *dead = pop_list[id]->dead.d;
        }
    }
}

void spoplib_retrieve(unsigned int id, char devtable, double *dev, double *size, unsigned int *limit) {
    if (pop_list[id] && pop_list[id]->accumulative) {
        unsigned int count = 0;
        for (count = 0; count < pop_list[id]->cat; count++) {
            if (pop_list[id]->individuals[count].accumulate) {
                quant_retrieve(pop_list[id]->individuals[count].accumulate, devtable, dev, size, limit);
            }
        }
    }
}

void spoplib_print(unsigned int id) {
    if (pop_list[id])
        spop_print(pop_list[id]);
}

void spoplib_destroy(unsigned int id) {
    if (pop_list[id])
        spop_destroy(&(pop_list[id]));
    pop_list[id] = 0;
}

void spoplib_destroy_all() {
    unsigned int i;
    for (i=0; i<MAXPOPS; i++)
        spoplib_destroy(i);
    free(pop_list);
}
