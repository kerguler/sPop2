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
#define realloc(pointer,size) printmem(1,(void *)realloc((pointer),(size)),(size),__FILE__,__LINE__)
#define calloc(number,size) printmem(1,(void *)calloc((number),(size)),(number)*(size),__FILE__,__LINE__)
#define malloc(size) printmem(1,(void *)malloc((size)),(size),__FILE__,__LINE__)
#define free(pointer) {printmem(0,0,sizeof(*(pointer)),__FILE__,__LINE__); free((pointer));}
void *printmem(int, void *, size_t, char *, int);

void *printmem(int type, void *pointer, size_t size, char *filen, int linen) {
    static size_t total = 0;
    total += (type?1:-1)*size;
    printf("%d: %s: %d: %d: %d\n",type,filen,linen,(int)(size),(int)(total));
    return pointer;
}
*/
int qunits = 0;
void set_qunits(int more) {
    qunits += more;
    printf("qunits = %d\n",qunits);
}

#define ZERO ((double)(0.0))
#define ONE ((double)(1.0))
/*
 * WARNING: If this is not high enough,
 * not all the developed will be transferred to the next stage!
 *
 * QUESTION: How high is high enough?
 */
#define ACCTHR ((double)(2.0))

double QSIZE_MAX = 10000.0;
double QSIZE_ROUND_EPS = 0.0;
double QSIZE_EPS = 0.0;
void set_APPROX(double eps) {
    QSIZE_ROUND_EPS = eps;
}

double QSIZE_ROUND(double val) {
    return round(val / QSIZE_ROUND_EPS) * QSIZE_ROUND_EPS;
}

gsl_rng *RANDOM = 0;

double fun_pois_C(double x, double par) {
    return gsl_cdf_poisson_P((unsigned int)x, ONE/par);
}

double fun_fixed_C(double x, double par) {
    return (unsigned int)x >= par;
}

double fun_Caswell_C(double x, double par) {
    return ((unsigned int)x == 0) ? ONE - par : ONE;
}

double fun_unif_C(double x, double par) {
    return ONE - pow(par, x + ONE);
}

double fun_cpois_C(double x, double par) {
    return gsl_sf_gamma_inc_Q(x + ONE, ONE / par);
}

char spop2_get_cfun(char mode, pfunc *cfun) {
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

spop2 spop2_init(unsigned char stochastic, unsigned char pdist) {
    // signal(SIGSEGV, clean_exit_on_sig);
    //
    spop2 pop = (spop2) calloc(1, sizeof(struct spop2_st));
    if (!pop) return 0;
    pop->devc = 0;
    pop->devtable = 0;
    pop->stochastic = stochastic;
    pop->pdist = pdist;
    if (!spop2_get_cfun(pop->pdist,&(pop->cfun))) return 0;
    if (stochastic) {
        RANDOM = get_RAND_GSL();
        pop->size.i = 0;
        pop->dead.i = 0;
        pop->developed.i = 0;
    } else {
        pop->size.d = 0.0;
        pop->dead.d = 0.0;
        pop->developed.d = 0.0;
    }
    return pop;
}

char spop2_empty_devc(qunit *dev) {
    if (!(*dev)) return 0;
    qunit p, tmp;
    HASH_ITER(hh, (*dev), p, tmp) {
        HASH_DEL((*dev), p);
        qunit_free(&p);
    }
    qunit_free(dev);
    return 0;
}

char spop2_empty(spop2 s) {
    char ret = spop2_empty_devc(&(s->devc));
    if (ret) return 1;
    ret += spop2_empty_devc(&(s->devtable));
    if (ret) return 1;
    //
    if (s->stochastic) {
        s->size.i = 0;
        s->dead.i = 0;
        s->developed.i = 0;
    } else {
        s->size.d = 0.0;
        s->dead.d = 0.0;
        s->developed.d = 0.0;
    }
    //
    return ret;
}

char spop2_destroy(spop2 *s) {
    char ret = 0;
    if (*s) {
        ret += spop2_empty(*s);
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

char spop2_sdadd(spop2 pop, double dev, sdnum size) {
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
        pop->size.i += size.i;
    else
        pop->size.d += size.d;
    //
    return ret;
}

char spop2_sdpopadd(spop2 pop, spop2 add, char devtable, sdnum *size) {
    if (pop->stochastic != add->stochastic) {
        printf("Error: Incompatible population types!\n");
        return 1;
    }
    if (pop->stochastic)
        size->i = 0;
    else
        size->d = 0.0;
    //
    qunit addc = devtable ? add->devtable : add->devc;
    qunit p = 0, tmp = 0, pp = 0;
    HASH_ITER(hh, addc, p, tmp) {
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
        if (pop->stochastic) {
            pop->size.i += p->size.i;
            size->i += p->size.i;
        } else {
            pop->size.d += p->size.d;
            size->d += p->size.d;
        }
    }
    return 0;
}

unsigned int spop2_nqunit(spop2 s) {
    unsigned int num = 0;
    if (!s) return num;
    qunit p, tmp;
    HASH_ITER(hh, s->devc, p, tmp) {
        ++num;
    }
    return num;
}

void spop2_print(spop2 s) {
    if (!s) return;
    printf("# m");
    PRINT_MODE(s->pdist);
    printf(" s%d n%g d%g c%g\n#dev,#size\n",
           s->stochastic,
           s->stochastic ? s->size.i : s->size.d,
           s->stochastic ? s->dead.i : s->dead.d,
           s->stochastic ? s->developed.i : s->developed.d);
    qunit p, tmp;
    HASH_ITER(hh, s->devc, p, tmp) {
        printf("%g,%g\n",
               p->dev,
               s->stochastic ? p->size.i : p->size.d);
    }
    if (s->devtable) {
        printf("# completed\n#dev,#size\n");
        HASH_ITER(hh, s->devtable, p, tmp) {
            printf("%g,%g\n",
                   p->dev,
                   s->stochastic ? p->size.i : p->size.d);
        }
    }
}

void spop2_retrieve(spop2 s, char devtable, double *dev, double *size, unsigned int *limit) {
    if (!s) return;
    unsigned int i = 0;
    qunit devc = devtable ? s->devtable : s->devc;
    qunit p, tmp;
    HASH_ITER(hh, devc, p, tmp) {
        dev[i] = p->dev;
        size[i] = s->stochastic ? (double)(p->size.i) : p->size.d;
        i++;
    }
    limit[0] = i;
}

char spop2_development(spop2 pop,
                       double gamma_k,
                       double gamma_theta,
                       pfunc cfun) {
    unsigned int dev = 0;
    double accd = 0.0;
    double haz = 0.0, h0 = 0.0, h1 = 0.0;
    qunit p = 0, tmp = 0;
    //
    sdnum item;
    qunit devc = 0, qnt = 0;
    qunit pp = 0;
    //
    if (pop->stochastic)
        pop->size.i = 0;
    else
        pop->size.d = 0.0;
    //
    if (!(HASH_COUNT(pop->devc))) return 0;
    unsigned int counter = 0;
    HASH_ITER(hh, pop->devc, p, tmp) {
        for (dev = 0;
             pop->stochastic ? p->size.i > 0 : p->size.d > QSIZE_EPS;
             dev++) {
            //
            accd = p->dev + ((double)dev / gamma_k);
            if (QSIZE_ROUND_EPS) accd = QSIZE_ROUND(accd);
            if (accd >= ACCTHR) {
                if (pop->stochastic) {
                    item.i = p->size.i;
                    pop->developed.i += p->size.i;
                    p->size.i = 0;
                } else {
                    item.d = p->size.d;
                    pop->developed.d += p->size.d;
                    p->size.d = 0.0;
                }
                accd = ACCTHR - ONE;
                HASH_FIND(hh, pop->devtable, &accd, sizeof(double), pp);
                if (pp) {
                    if (pop->stochastic)
                        pp->size.i += item.i;
                    else
                        pp->size.d += item.d;
                } else {
                    qnt = qunit_new(accd, item);
                    if (!qnt) return 1;
                    HASH_ADD(hh, pop->devtable, dev, sizeof(double), qnt);
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
                p->size.i -= item.i;
            } else {
                item.d = p->size.d * haz;
                if (!item.d) continue;
                p->size.d -= item.d;
            }
            //
            if (accd >= ONE) {
                if (pop->stochastic) {
                    pop->developed.i += item.i;
                } else {
                    pop->developed.d += item.d;
                }
                accd -= ONE;
                HASH_FIND(hh, pop->devtable, &accd, sizeof(double), pp);
                if (pp) {
                    if (pop->stochastic)
                        pp->size.i += item.i;
                    else
                        pp->size.d += item.d;
                } else {
                    qnt = qunit_new(accd, item);
                    if (!qnt) return 1;
                    HASH_ADD(hh, pop->devtable, dev, sizeof(double), qnt);
                }
            } else {
                if (pop->stochastic) {
                    pop->size.i += item.i;
                } else {
                    pop->size.d += item.d;
                }
                //
                HASH_FIND(hh, devc, &accd, sizeof(double), pp);
                if (pp) {
                    if (pop->stochastic)
                        pp->size.i += item.i;
                    else
                        pp->size.d += item.d;
                } else {
                    qnt = qunit_new(accd, item);
                    if (!qnt) return 1;
                    HASH_ADD(hh, devc, dev, sizeof(double), qnt);
                    counter++;
                }
            }
        }
        //
        HASH_DEL(pop->devc, p);
        qunit_free(&p);
    }
    //
    qunit_free(&(pop->devc));
    pop->devc = devc;
    //
    if (counter > QSIZE_MAX)
        printf("Warning: Hash size = %d\n", counter);
    //
    return 0;
}

char spop2_mortality(spop2 pop,
                     double prob,
                     char devtable) {
    sdnum item;
    qunit p, tmp;
    qunit devc = devtable ? pop->devtable : pop->devc;
    //
    if (!(HASH_COUNT(devc))) return 0;
    if (pop->stochastic) {
        HASH_ITER(hh, devc, p, tmp) {
            if (p->size.i) {
                item.i = gsl_ran_binomial(RANDOM,
                                          prob,
                                          p->size.i);
                p->size.i -= item.i;
                pop->size.i -= item.i;
                if (!devtable)
                    pop->dead.i += item.i;
            }
            if (p->size.i == 0) {
                HASH_DEL(devc, p);
                qunit_free(&p);
            }
        }
    } else {
        HASH_ITER(hh, devc, p, tmp) {
            if (p->size.d) {
                item.d = p->size.d * prob;
                p->size.d -= item.d;
                pop->size.d -= item.d;
                if (!devtable)
                    pop->dead.d += item.d;
            }
            if (p->size.d == ZERO) {
                HASH_DEL(devc, p);
                qunit_free(&p);
            }
        }
    }
    if (devtable)
        pop->devtable = devc;
    else
        pop->devc = devc;
    return 0;
}

char spop2_iterate(spop2 pop,
                   double dev_mean,       // mean development time
                   double dev_sd,         // sd development time
                   double death,          // probability of death per iteration
                   char devtable) {
    double gamma_k = 0.0,
           gamma_theta = 0.0;
    char pdist = pop->pdist;
    pfunc cfun = pop->cfun;
    //
    if (pop->stochastic) {
        pop->developed.i = 0;
        pop->dead.i = 0;
    } else {
        pop->developed.d = 0.0;
        pop->dead.d = 0.0;
    }
    spop2_empty_devc(&(pop->devtable));
    if (!pop->devc) return 1;
    //
    if (dev_sd == 0) {
        pdist = MODE_ACCP_FIXED;
        if (!spop2_get_cfun(pdist, &cfun)) return 1;
    }
    switch (pdist) {
        case MODE_ACCP_ERLANG:
            gamma_theta = dev_sd * dev_sd / dev_mean;
            gamma_k = dev_mean / gamma_theta;
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                gamma_theta = dev_mean / gamma_k;
            }
            break;
        case MODE_ACCP_GAMMA:
            gamma_theta = dev_sd * dev_sd / dev_mean;
            gamma_k = dev_mean / gamma_theta;
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                gamma_theta = dev_mean / gamma_k;
            }
            break;
        case MODE_ACCP_FIXED:
            gamma_k = round(dev_mean);
            gamma_theta = 1.0;
            break;
        case MODE_ACCP_CASWELL:
            gamma_theta = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta >= 1.0 || gamma_theta == 0.0) {
                printf("Error: The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 1;
            }
            gamma_k = dev_mean * gamma_theta / (1.0 - gamma_theta);
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                gamma_theta = dev_mean / gamma_k;
            }
            break;
        case MODE_ACCP_PASCAL:
            gamma_theta = dev_mean / (dev_sd * dev_sd);
            if (gamma_theta >= 1.0 || gamma_theta == 0.0) {
                printf("Error: The negative binomial cannot yield mean=%g and sd=%g\n", dev_mean, dev_sd);
                return 1;
            }
            gamma_k = dev_mean * gamma_theta / (1.0 - gamma_theta);
            if (gamma_k != round(gamma_k)) {
                gamma_k = round(gamma_k);
                gamma_theta = gamma_k / (dev_mean + gamma_k);
            }
            break;
        default:
            printf("Error: I don't know what to do with this: %d\n", pdist);
            return 1;
            break;
    }
    //
    if (gamma_k == 0) {
        printf("Error: k=0 is not acceptable!\n");
        return 1;
    }
    //
    if (devtable) {
        spop2_development(pop, gamma_k, gamma_theta, cfun);
        if (death)
            spop2_mortality(pop, death, devtable);
    } else {
        if (death)
            spop2_mortality(pop, death, devtable);
        spop2_development(pop, gamma_k, gamma_theta, cfun);
    }
    //
    return 0;
}

