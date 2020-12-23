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

#define ONE ((double)(1.0))

double QSIZE_MAX = 10000.0;
#define QSIZE_MAX_EPS   1e-3
#define QSIZE_MAX_NONE  0.0

double QSIZE_ROUND_EPS = QSIZE_MAX_NONE;
double QSIZE_EPS = 0.0;
void set_APPROX(char on) {
    if (on) {
        QSIZE_ROUND_EPS = QSIZE_MAX_EPS;
        // QSIZE_EPS = 1e-8;
    } else {
        QSIZE_ROUND_EPS = QSIZE_MAX_NONE;
        // QSIZE_EPS = 0.0;
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

char quant_get_cfun(char mode, pfunc *cfun) {
    switch (mode) {
        case MODE_ACCP_ERLANG:
            *cfun = fun_pois_C;
            break;
        case MODE_ACCP_FIXED:
            *cfun = fun_fixed_C;
            break;
        case MODE_ACCP_PASCAL:
            *cfun = fun_unif_C;
            break;
        case MODE_ACCP_GAMMA:
            *cfun = fun_cpois_C;
            break;
        case MODE_ACCP_CASWELL:
            *cfun = fun_Caswell_C;
            break;
        default:
            printf("Error: I don't know what to do with this: %d\n", mode);
            return 0;
            break;
    }
    return 1;
}

quant quant_init(unsigned char stochastic, unsigned char pdist) {
    // signal(SIGSEGV, clean_exit_on_sig);
    //
    quant pop = (quant) calloc(1, sizeof(struct quant_st));
    if (!pop) return 0;
    pop->devc = 0;
    pop->stochastic = stochastic;
    pop->pdist = pdist;
    if (!quant_get_cfun(pop->pdist,&(pop->cfun))) return 0;
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
        printf("Error: Incompatible population types!\n");
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

void quant_retrieve(quant s, double *dev, double *size, unsigned int *limit) {
    if (!s) return;
    unsigned int i = 0;
    qunit p, tmp;
    HASH_ITER(hh, s->devc, p, tmp) {
        dev[i] = p->dev;
        size[i] = s->stochastic ? (double)(p->size.i) : p->size.d;
        printf("Retrieve %d %g %g\n",i,dev[i],size[i]);
        i++;
    }
    limit[0] = i;
}

char quant_iterate_hazards(quant pop,
                           double gamma_k,
                           double gamma_theta,
                           pfunc cfun) {
    if (pop->stochastic) {
        pop->completed.i = 0;
        pop->size.i = 0;
    } else {
        pop->completed.d = 0.0;
        pop->size.d = 0.0;
    }
    //
    unsigned int dev = 0;
    double accd = 0.0;
    double haz = 0.0, h0 = 0.0, h1 = 0.0;
    qunit p = 0, tmp = 0;
    //
    sdnum item;
    qunit devc = 0, qnt = 0;
    qunit pp = 0;
    //
    unsigned int counter = 0;
    HASH_ITER(hh, pop->devc, p, tmp) {
        for (dev = 0;
             pop->stochastic ? p->size.i > 0 : p->size.d > QSIZE_EPS;
             dev++) {
            //
            accd = p->dev + ((double)dev / gamma_k);
            if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
            if (accd >= ONE) {
                if (pop->stochastic) {
                    pop->completed.i += p->size.i;
                    p->size.i = 0;
                } else {
                    pop->completed.d += p->size.d;
                    p->size.d = 0.0;
                }
                break;
            }
            //
            h0 = !dev ? 0.0 : (cfun)(dev - 1, gamma_theta);
            h1 = (cfun)(dev, gamma_theta);
            haz = h0 == ONE ? 1.0 : (h1 - h0) / (ONE - h0);
            //
            if (pop->stochastic) {
                item.i = gsl_ran_binomial(RANDOM,
                                          haz,
                                          p->size.i);
                if (!item.i) continue;
                pop->size.i += item.i;
                p->size.i -= item.i;
            } else {
                item.d = p->size.d * haz;
                if (!item.d) continue;
                pop->size.d += item.d;
                p->size.d -= item.d;
            }
            //
            qnt = qunit_new(accd, item);
            if (!qnt) return 1;
            HASH_FIND(hh, devc, &accd, sizeof(double), pp);
            if (pp) {
                if (pop->stochastic)
                    pp->size.i += item.i;
                else
                    pp->size.d += item.d;
            } else {
                HASH_ADD(hh, devc, dev, sizeof(double), qnt);
                counter++;
            }
        }
        //
        HASH_DEL(pop->devc, p);
        free(p);
    }
    pop->devc = devc;
    //
    if (counter > QSIZE_MAX)
        printf("Warning: Hash size = %d\n", counter);
    //
    return 0;
}

char quant_iterate_twin_hazards(quant pop,
                               double *gamma_k,
                               double *gamma_theta,
                               double *gamma_p,
                               pfunc cfun) {
    if (pop->stochastic) {
        pop->completed.i = 0;
        pop->size.i = 0;
    } else {
        pop->completed.d = 0.0;
        pop->size.d = 0.0;
    }
    //
    unsigned int dev = 0, acc = 0;
    double accd = 0.0;
    double haz[2] = {0.0, 0.0}, h0 = 0.0, h1 = 0.0;
    qunit p = 0, tmp = 0;
    //
    sdnum item;
    qunit devc = 0, qnt = 0;
    qunit pp = 0;
    //
    unsigned int counter = 0;
    HASH_ITER(hh, pop->devc, p, tmp) {
        for (dev = 0;
             pop->stochastic ? p->size.i > 0 : p->size.d > QSIZE_EPS;
             dev++) {
            //
            accd = p->dev + ((double)dev/gamma_k[0]);
            if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
            if (accd >= ONE) {
                if (pop->stochastic) {
                    pop->completed.i += p->size.i;
                    p->size.i = 0;
                } else {
                    pop->completed.d += p->size.d;
                    p->size.d = 0.0;
                }
                break;
            }
            //
            h0 = !dev ? 0.0 : (cfun)(dev - 1, gamma_theta[0]);
            h1 = (cfun)(dev, gamma_theta[0]);
            haz[0] = gamma_p[0] * (h0 == ONE ? 1.0 : (h1 - h0) / (ONE - h0));
            //
            for (acc = 0;
                 pop->stochastic ? p->size.i > 0 : p->size.d > QSIZE_EPS;
                 acc++) {
                //
                accd += p->dev + ((double)acc/gamma_k[1]);
                if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
                if (accd >= ONE) {
                    haz[1] = gamma_p[1] * (ONE - (!acc ? 0.0 : (cfun)(acc - 1, gamma_theta[1])));
                    //
                    if (pop->stochastic) {
                        item.i = gsl_ran_binomial(RANDOM,
                                                  (haz[0]+haz[1]),
                                                  p->size.i);
                        if (!item.i) continue;
                        pop->completed.i += item.i;
                        p->size.i -= item.i;
                    } else {
                        item.d = p->size.d * (haz[0]+haz[1]);
                        if (!item.d) continue;
                        pop->completed.d += item.d;
                        p->size.d -= item.d;
                    }
                    //
                    break;
                }
                //
                h0 = !acc ? 0.0 : (cfun)(acc - 1, gamma_theta[1]);
                h1 = (cfun)(acc, gamma_theta[1]);
                haz[1] = gamma_p[1] * (h0 == ONE ? 1.0 : (h1 - h0) / (ONE - h0));
                //
                if (pop->stochastic) {
                    item.i = gsl_ran_binomial(RANDOM,
                                              (haz[0]+haz[1]),
                                              p->size.i);
                    if (!item.i) continue;
                    pop->size.i += item.i;
                    p->size.i -= item.i;
                } else {
                    item.d = p->size.d * (haz[0]+haz[1]);
                    if (!item.d) continue;
                    pop->size.d += item.d;
                    p->size.d -= item.d;
                }
                //
                qnt = qunit_new(accd, item);
                if (!qnt) return 1;
                HASH_FIND(hh, devc, &accd, sizeof(double), pp);
                if (pp) {
                    if (pop->stochastic)
                        pp->size.i += item.i;
                    else
                        pp->size.d += item.d;
                } else {
                    HASH_ADD(hh, devc, dev, sizeof(double), qnt);
                    counter++;
                }
            }
        }
        //
        HASH_DEL(pop->devc, p);
        free(p);
    }
    /*
     * MUST HAVE TWO ALTERNATIVE DEVELOPMENT PATHS!
     * pop->devc0;
     * pop->devc1;
     */
    pop->devc = devc;
    //
    if (counter > QSIZE_MAX)
        printf("Warning: Hash size = %d\n", counter);
    //
    return 0;
}

char quant_iterate(quant pop,
                   double dev_mean,       // mean development time
                   double dev_sd) {       // sd development time
    double gamma_k[2] = {0.0, 0.0},
            gamma_theta[2] = {0.0, 0.0},
            gamma_p[2] = {1.0, 0.0};
    unsigned int gamma_i = 1;
    char pdist = pop->pdist;
    pfunc cfun = pop->cfun;
    if (dev_sd == 0) {
        pdist = MODE_ACCP_FIXED;
        if (!quant_get_cfun(pdist, &cfun)) return 1;
    }
    switch (pdist) {
        case MODE_ACCP_ERLANG:
            gamma_theta[0] = dev_sd * dev_sd / dev_mean;
            gamma_k[0] = dev_mean / gamma_theta[0];
            if (gamma_k[0] != round(gamma_k[0])) {
                gamma_k[0] = round(gamma_k[0]);
            }
            break;
        case MODE_ACCP_GAMMA:
            gamma_theta[0] = dev_sd * dev_sd / dev_mean;
            gamma_k[0] = dev_mean / gamma_theta[0];
            if (gamma_k[0] != round(gamma_k[0])) {
                static double k;
                gamma_i = 2;
                k = gamma_k[0];
                gamma_k[1] = ceil(gamma_k[0]);
                gamma_k[0] = floor(gamma_k[0]);
                gamma_p[0] = (gamma_k[1]-k) / (gamma_k[1]-gamma_k[0]);
                gamma_p[1] = 1.0 - gamma_p[0];
                gamma_theta[0] = dev_mean / gamma_k[0];
                gamma_theta[1] = dev_mean / gamma_k[1];
                printf("Gamma: %g %g, %g %g, %g %g\n",
                       gamma_k[0],
                       gamma_k[1],
                       gamma_p[0],
                       gamma_p[1],
                       gamma_theta[0],
                       gamma_theta[1]
                );
            }
            break;
        case MODE_ACCP_FIXED:
            gamma_k[0] = round(dev_mean);
            gamma_theta[0] = 1.0;
            break;
        case MODE_ACCP_CASWELL:
            // gamma_theta[0] = dev_mean / (dev_mean + (dev_sd * dev_sd));
            gamma_theta[0] = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta[0] >= 1.0 || gamma_theta[0] == 0.0) {
                printf("Error: The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 1;
            }
            // gamma_k[0] = (dev_mean * dev_mean) / (dev_mean + (dev_sd * dev_sd));
            gamma_k[0] = dev_mean * gamma_theta[0] / (1.0 - gamma_theta[0]);
            if (gamma_k[0] != round(gamma_k[0])) {
                gamma_k[0] = round(gamma_k[0]);
            }
            break;
        case MODE_ACCP_PASCAL:
            gamma_theta[0] = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta[0] >= 1.0 || gamma_theta[0] == 0.0) {
                printf("Error: The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 1;
            }
            gamma_k[0] = dev_mean * gamma_theta[0] / (1.0 - gamma_theta[0]);
            if (gamma_k[0] != round(gamma_k[0])) {
                gamma_k[0] = round(gamma_k[0]);
            }
            break;
        default:
            printf("Error: I don't know what to do with this: %d\n", pdist);
            return 1;
            break;
    }
    //
    if (gamma_k[0] == 0) { // This should be the minimum one!
        printf("Error: 0 mean is not acceptable!\n");
        return 1;
    }
    //
    if (pop->stochastic) {
        pop->completed.i = 0;
    } else {
        pop->completed.d = 0.0;
    }
    //
    if (!pop->devc) return 1;
    //
    if (gamma_i == 1)
        return quant_iterate_hazards(pop, gamma_k[0], gamma_theta[0], cfun);
    else
        return quant_iterate_twin_hazards(pop, gamma_k, gamma_theta, gamma_p, cfun);
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
