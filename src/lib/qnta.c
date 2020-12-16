#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "spop2.h"

/*
#include <signal.h>
void (*signal(int signum, void (*sighandler)(int)))(int);
void clean_exit_on_sig(int sig_num) {
    printf ("\n Signal %d received\n",sig_num);
    exit(1);
}
*/

/*
#define free(pointer) {printmem(0,0,sizeof(*(pointer)),__FILE__,__LINE__); free((pointer));}
void *printmem(int, void *, size_t, char *, int);

void *printmem(int type, void *pointer, size_t size, char *filen, int linen) {
    static size_t total = 0;
    total += (type?1:-1)*size;
    printf("%d: %s: %d: %d: %d\n",type,filen,linen,(int)(size),(int)(total));
    return pointer;
}
*/

double QSIZE_MAX = 10000.0;
#define QSIZE_MAX_EPS   1e-3
#define QSIZE_MAX_NONE  0.0

double QSIZE_ROUND_EPS = QSIZE_MAX_NONE;
double QSIZE_EPS = 0.0;
void set_APPROX(char on) {
    if (on) {
        QSIZE_ROUND_EPS = QSIZE_MAX_EPS;
        QSIZE_EPS = 1e-8;
    } else {
        QSIZE_ROUND_EPS = QSIZE_MAX_NONE;
        QSIZE_EPS = 0.0;
    }
}

double QSIZE_ROUND(double val) {
    return round(val / QSIZE_ROUND_EPS) * QSIZE_ROUND_EPS;
}

gsl_rng *RANDOM = 0;

double fun_pois_C(double x, double par) {
    return gsl_cdf_poisson_P((unsigned int)x, 1.0/par);
}

double fun_fixed_C(double x, double par) {
    return (unsigned int)x >= par;
}

double fun_Caswell_C(double x, double par) {
    return ((unsigned int)x == 0) ? 1.0-par : 1.0;
}

double fun_unif_C(double x, double par) {
    return 1.0-pow(par, x+1.0);
}

double fun_cpois_C(double x, double par) {
    return gsl_sf_gamma_inc_Q(x + 1.0, 1.0 / par);
}

int aorder (const void *a, const void *b) {
    return ( *(double *)a > *(double *)b ) ? 1 : ( ( *(double *)a < *(double *)b ) ? -1 : 0 );
}

char fun_rdist (unsigned int pop, double *p, unsigned int size, unsigned int *ret) {
    if (!RANDOM)
        RANDOM = get_RAND_GSL();
    double max=0.0;
    unsigned int i=0, j=0;
    //
    for (i=0; i<size; i++) ret[i] = 0;
    max = p[size-1];
    //
    double *vec = (double *)calloc(pop,sizeof(double));
    if (!vec) return 1;
    for (i=0; i<pop; i++)
        vec[i] = max * gsl_rng_uniform(RANDOM);
    qsort(vec, pop, sizeof(double), aorder);
    //
    for (i=0, j=0; i<pop; ) {
        if (vec[i] <= p[j]) {
            ret[j] += 1;
            i++;
            continue;
        }
        j++;
    }
    //
    free(vec);
    return 0;
}

quant quant_init(unsigned char stochastic, unsigned char pdist) {
    /*
     * signal(SIGSEGV, clean_exit_on_sig);
     */
    set_APPROX(0);
    //
    quant pop = (quant) calloc(1, sizeof(struct quant_st));
    if (!pop) return 0;
    pop->devc = 0;
    pop->stochastic = stochastic;
    pop->pdist = pdist;
    switch (pop->pdist) {
        case MODE_ACCP_ERLANG:
            pop->cfun = fun_pois_C;
            break;
        case MODE_ACCP_FIXED:
            pop->stochastic = 0;
            pop->cfun = fun_fixed_C;
            break;
        case MODE_ACCP_PASCAL:
            pop->cfun = fun_unif_C;
            break;
        case MODE_ACCP_GAMMA:
            pop->cfun = fun_cpois_C;
            break;
        case MODE_ACCP_CASWELL:
            pop->cfun = fun_Caswell_C;
            break;
        default:
            printf("I don't know what to do with this: %d\n", pop->pdist);
            return 0;
            break;
    }
    if (stochastic) {
        RANDOM = get_RAND_GSL();
        pop->size.i = 0;
        pop->completed.i = 0;
    } else {
        pop->size.d = 0.0;
        pop->completed.d = 0.0;
    }
    return pop;
}

char quant_empty(quant s) {
    qunit p, tmp;
    HASH_ITER(hh, s->devc, p, tmp) {
        HASH_DEL(s->devc, p);
        free(p);
    }
    //
    if (s->stochastic) {
        s->size.i = 0;
        s->completed.i = 0;
    } else {
        s->size.d = 0.0;
        s->completed.d = 0.0;
    }
    //
    return 0;
}

char quant_destroy(quant *s) {
    char ret = 0;
    if (*s) {
        ret += quant_empty(*s);
        free((*s));
        (*s) = 0;
    }
    return ret;
}

qunit qunit_new(double dev, sdnum size) {
    qunit tmp = (qunit)calloc(1,sizeof(struct qunit_st));
    if (!tmp) return 0;
    tmp->dev = dev;
    tmp->size = size;
    return tmp;
}

void qunit_free(qunit *tmp) {
    if (!(*tmp)) return;
    free((*tmp));
    (*tmp) = 0;
}

char quant_sdadd(quant pop, double dev, sdnum size) {
    char ret = 0;
    //
    if (QSIZE_ROUND_EPS) dev = QSIZE_ROUND(dev);
    //
    qunit qnt;
    HASH_FIND(hh, pop->devc, &dev, sizeof(double), qnt);
    if (qnt) {
        if (pop->stochastic)
            qnt->size.i += size.i;
        else
            qnt->size.d += size.d;
    } else {
        qnt = qunit_new(dev, size);
        if (!qnt) return 0;
        HASH_ADD(hh, pop->devc, dev, sizeof(double), qnt);
    }
    //
    if (pop->stochastic)
        pop->size.i += qnt->size.i;
    else
        pop->size.d += qnt->size.d;
    //
    return ret;
}

char quant_sdpopadd(quant pop, quant add) {
    if (pop->stochastic != add->stochastic) {
        printf("Incompatible population types!\n");
        return 1;
    }
    qunit p = 0, tmp = 0, pp = 0;
    HASH_ITER(hh, add->devc, p, tmp) {
        HASH_FIND(hh, pop->devc, &(p->dev), sizeof(double), pp);
        if (pp) {
            if (pop->stochastic)
                pp->size.i += p->size.i;
            else
                pp->size.d += p->size.d;
        } else {
            qunit ppp = qunit_new(p->dev,p->size);
            HASH_ADD(hh, pop->devc, dev, sizeof(double), ppp);
        }
        //
        if (pop->stochastic)
            pop->size.i += p->size.i;
        else
            pop->size.d += p->size.d;
    }
    return 0;
}

void quant_print(quant s) {
    if (!s) return;
    printf("# s%d n%g c%g\n#dev,#size\n",
           s->stochastic,
           s->stochastic ? s->size.i : s->size.d,
           s->stochastic ? s->completed.i : s->completed.d);
    qunit p, tmp;
    HASH_ITER(hh, s->devc, p, tmp) {
        printf("%g,%g\n",
               p->dev,
               s->stochastic ? p->size.i : p->size.d);
    }
}

char quant_iterate_stochastic(quant pop,
                              double gamma_k,
                              double gamma_theta) {
    pop->completed.i = 0;
    unsigned int item = 0;
    double accd = 0.0;
    qunit p = 0, tmp = 0;
    //
    qunit devc = 0, qnt = 0;
    qunit pp = 0;
    unsigned int dev = 0, i = 0;
    sdnum itm;
    //
    unsigned int counter = 0;
    HASH_ITER(hh, pop->devc, p, tmp) {
        dev = floor(p->dev * gamma_k);
        item = gsl_ran_binomial(RANDOM,
                                1.0 - (pop->cfun)(gamma_k - dev - 1, gamma_theta),
                                p->size.i);
        p->size.i -= item;
        pop->size.i -= item;
        pop->completed.i += item;
        //
        if (p->size.i > 0) {
            unsigned int rsize = gamma_k - dev;
            unsigned int *vec = (unsigned int *) calloc(rsize, sizeof(unsigned int));
            if (!vec) return 1;
            double *prob = (double *) calloc(rsize, sizeof(double));
            if (!prob) {
                free(vec);
                return 1;
            }
            for (i = 0; i < rsize; i++)
                prob[i] = (pop->cfun)(i, gamma_theta);
            if (fun_rdist(p->size.i, prob, rsize, vec)) {
                free(prob);
                free(vec);
                return 1;
            }
            for (i = 0; i < rsize; i++) {
                if (!vec[i]) continue;
                accd = p->dev + ((double) i / gamma_k);
                //
                if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
                //
                itm.i = vec[i];
                if (itm.i) {
                    qnt = qunit_new(accd, itm);
                    if (!qnt) return 1;
                    HASH_FIND(hh, devc, &accd, sizeof(double), pp);
                    if (pp)
                        pp->size.i += itm.i;
                    else {
                        HASH_ADD(hh, devc, dev, sizeof(double), qnt);
                        counter++;
                    }
                }
            }
        }
        //
        HASH_DEL(pop->devc, p);
        free(p);
    };
    //
    if (counter > QSIZE_MAX) {
        set_APPROX(1);
        printf("Warning! Hash size = %d\n",counter);
    } else {
        set_APPROX(0);
    }
    pop->devc = devc;
    return 0;
}

char quant_iterate_deterministic(quant pop,
                                 double gamma_k,
                                 double gamma_theta) {
    pop->completed.d = 0.0;
    unsigned int acc = 0;
    double accd = 0.0;
    qunit p = 0, tmp = 0;
    //
    pop->size.d = 0.0;
    sdnum item;
    unsigned int dev = 0;
    qunit devc = 0, qnt = 0;
    qunit pp = 0;
    //
    unsigned int counter = 0;
    HASH_ITER(hh, pop->devc, p, tmp) {
        if (p->size.d > QSIZE_EPS) {
            dev = floor(p->dev * gamma_k);
            pop->completed.d += p->size.d * (1.0 - (pop->cfun)(gamma_k - dev - 1, gamma_theta));
            //
            for (acc = dev; acc < gamma_k; acc++) {
                item.d = p->size.d * ((pop->cfun)(acc - dev, gamma_theta) -
                                      (acc == dev ? 0.0 : (pop->cfun)(acc - dev - 1, gamma_theta)));
                accd = p->dev + ((double) (acc - dev) / gamma_k);
                //
                if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
                //
                if (item.d) {
                    qnt = qunit_new(accd, item);
                    if (!qnt) return 1;
                    HASH_FIND(hh, devc, &accd, sizeof(double), pp);
                    if (pp)
                        pp->size.d += item.d;
                    else {
                        HASH_ADD(hh, devc, dev, sizeof(double), qnt);
                        counter++;
                    }
                    pop->size.d += item.d;
                }
            }
        }
        //
        HASH_DEL(pop->devc, p);
        free(p);
    };
    pop->devc = devc;
    //
    if (counter > QSIZE_MAX) {
        set_APPROX(1);
        printf("Warning! Hash size = %d\n",counter);
    } else {
        set_APPROX(0);
    }
    //
    return 0;
}

char quant_iterate(quant pop,
                   double dev_mean,       // mean development time
                   double dev_sd) {       // sd development time
    double gamma_k = 0.0,
            gamma_theta = 0.0;
    switch (pop->pdist) {
        case MODE_ACCP_ERLANG:
            gamma_theta = dev_sd * dev_sd / dev_mean;
            gamma_k = dev_mean / gamma_theta;
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                /*
                double m = gamma_k * gamma_theta;
                double s = sqrt(gamma_theta * m);
                printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n", m, s);
                */
            }
            break;
        case MODE_ACCP_GAMMA:
            gamma_theta = dev_sd * dev_sd / dev_mean;
            gamma_k = dev_mean / gamma_theta;
            // printf("k=%g, theta=%g\n",gamma_k,gamma_theta);
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                /*
                double m = gamma_k * gamma_theta;
                double s = sqrt(gamma_theta * m);
                printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n", m, s);
                */
            }
            break;
        case MODE_ACCP_FIXED:
            gamma_k = round(dev_mean);
            gamma_theta = 1.0;
            break;
        case MODE_ACCP_CASWELL:
            // gamma_theta = dev_mean / (dev_mean + (dev_sd * dev_sd));
            gamma_theta = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta >= 1.0 || gamma_theta == 0.0) {
                printf("The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 0;
            }
            // gamma_k = (dev_mean * dev_mean) / (dev_mean + (dev_sd * dev_sd));
            gamma_k = dev_mean * gamma_theta / (1.0 - gamma_theta);
            // printf("k=%g, theta=%g\n",gamma_k,gamma_theta);
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                /*
                double m = gamma_k * (1.0 - gamma_theta) / gamma_theta;
                double s = sqrt(dev_mean / gamma_theta);
                printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n", m, s);
                */
            }
            break;
        case MODE_ACCP_PASCAL:
            gamma_theta = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta >= 1.0 || gamma_theta == 0.0) {
                printf("The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 0;
            }
            gamma_k = dev_mean * gamma_theta / (1.0 - gamma_theta);
            // printf("k=%.18f, theta=%.18f\n",gamma_k,gamma_theta);
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                /*
                double m = gamma_k * (1.0 - gamma_theta) / gamma_theta;
                double s = sqrt(dev_mean / gamma_theta);
                printf("For the MODE_ACCP_PASCAL distribution, the shape parameter will be adjusted to yield\nk = %.18f, theta = %.18f, Mean = %.18f, St.dev. = %.18f\n", gamma_k, gamma_theta, m, s);
                */
            }
            // printf("p=%g, r=%g\n",gamma_theta,gamma_k);
            break;
        default:
            printf("I don't know what to do with this: %d\n", pop->pdist);
            return 0;
            break;
    }
    //
    if (pop->stochastic) {
        pop->completed.i = 0;
    } else {
        pop->completed.d = 0.0;
    }
    //
    if (!pop->devc) return 0;
    //
    if (pop->stochastic)
        return quant_iterate_stochastic(pop, gamma_k, gamma_theta);
    else
        return quant_iterate_deterministic(pop, gamma_k, gamma_theta);
    //
    return 0;
}

char quant_survive(quant pop, double prob, sdnum *ret) {
    ret->i = 0;
    ret->d = 0.0;
    pop->completed.i = 0;
    pop->completed.d = 0.0;
    sdnum item;
    qunit p, tmp;
    if (pop->stochastic) {
        HASH_ITER(hh, pop->devc, p, tmp) {
            if (!(p->size.i)) continue;
            item.i = gsl_ran_binomial(RANDOM,
                                      prob,
                                      p->size.i);
            p->size.i -= item.i;
            pop->size.i -= item.i;
            ret->i += item.i;
            if (!(p->size.i)) {
                HASH_DEL(pop->devc, p);
                free(p);
            }
        }
    } else {
        HASH_ITER(hh, pop->devc, p, tmp) {
            if (!(p->size.d)) continue;
            item.d = p->size.d * prob;
            p->size.d -= item.d;
            pop->size.d -= item.d;
            ret->d += item.d;
            if (!(p->size.d)) {
                HASH_DEL(pop->devc, p);
                free(p);
            }
        }
    }
    return 0;
}
