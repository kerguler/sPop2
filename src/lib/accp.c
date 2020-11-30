#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "spop2.h"

gsl_rng *RANDOM;

double fun_pois_C(double x, double par) {
  return gsl_cdf_poisson_P((unsigned int)x, 1.0/par);
}

double fun_fixed_C(double x, double par) {
  return (unsigned int)x > 0;
}

double fun_unif_C(double x, double par) {
  return 1.0-pow(par, x+1.0);
}

double fun_cpois_C(double x, double par) {
  return gsl_sf_gamma_inc_Q(x+1.0, 1.0/par);
}

int ascending (const void *a, const void *b) {
   return ( *(double *)a > *(double *)b ) ? 1 : ( ( *(double *)a < *(double *)b ) ? -1 : 0 );
}

char choose(unsigned int pop, double *p, unsigned int size, unsigned int *ret) {
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
  qsort(vec, pop, sizeof(double), ascending);
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

accp accp_init(unsigned char stochastic, unsigned char pdist) {
  accp pop = (accp)calloc(1,sizeof(struct accumulative_st));
  if (!pop) return 0;
  pop->first = 0;
  pop->last = 0;
  pop->cat = 0;
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
  default:
    printf("I don't know what to do with this: %d\n",pop->pdist);
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

char accp_empty(accp s) {
  char ret = chain_resize(s,0);
  if (s->stochastic) {
    s->size.i = 0;
    s->completed.i = 0;
  } else {
    s->size.d = 0.0;
    s->completed.d = 0.0;
  }
  return ret;
}

char accp_destroy(accp *s) {
  char ret = 0;
  if (*s) {
    ret += accp_empty(*s);
    free((*s));
    (*s) = 0;
  }
  return ret;
}

sdnum accp_get_size(accp s) {
  if (s->stochastic)
    s->size.i = 0;
  else
    s->size.d = 0.0;
  if (!(s->first))
    return s->size;
  //
  chain tmp;
  for (tmp = s->first; tmp; tmp = tmp->next) {
    if (s->stochastic)
      s->size.i += tmp->size.i;
    else
      s->size.d += tmp->size.d;
  }
  return s->size;
}

char accp_sdadd(accp pop, unsigned int stage, sdnum size) {
    char ret = 0;
    if (stage >= pop->cat) ret += chain_resize(pop,stage+1);
    if (ret) return ret;
    //
    unsigned int count = 0;
    chain tmp;
    for (count=0, tmp=pop->first; count<stage && tmp; count++, tmp=tmp->next);
    if (pop->stochastic) {
      tmp->size.i += size.i;
      pop->size.i += size.i;
    } else {
      tmp->size.d += size.d;
      pop->size.d += size.d;
    }
    //
    return ret;
}

void accp_print(accp s) {
  if (!s) return;
  unsigned int count;
  chain tmp;
  printf(">------------------>\n Population type: ");
  if (s->stochastic) {
    printf("Stochastic\n Population size: %d\n Completed: %d\n Dev.stage\tNumber\n ---------\t------\n",
           s->size.i,
           s->completed.i);
    for (tmp=s->first, count=0; tmp; tmp=tmp->next, count++)
      printf(" %d          \t%d\n",
             count,
             tmp->size.i);
  } else {
    printf("Deterministic\n Population size: %g\n Completed: %g\n Dev.stage\tNumber\n ---------\t------\n",
           s->size.d,
           s->completed.d);
    for (tmp=s->first, count=0; tmp; tmp=tmp->next, count++)
      printf(" %d          \t%g\n",
             count,
             tmp->size.d);
  }
  printf(" ---------\t------\n Number of stages to complete: %d\n>------------------>\n",s->cat);
}

char chain_resize(accp pop, unsigned int cat) {
  // Do nothing
  if (cat == pop->cat) return 0;
  // Add more links
  if (cat > pop->cat) {
    chain tmp;
    while (cat > pop->cat) {
      tmp = (chain)calloc(1,sizeof(struct stage_st));
      if (!tmp) return 1;
      tmp->prev = pop->last;
      tmp->next = 0;
      if (!(pop->first)) {
        pop->first = pop->last = tmp;
      } else {
        pop->last->next = tmp;
        pop->last = tmp;
      }
      pop->cat += 1;
    }
    return 0;
  }
  // Remove links and store the ones completed
  {
    chain tmp;
    for (tmp=pop->last; cat < pop->cat && tmp; tmp=pop->last) {
      if (pop->stochastic) {
        pop->size.i -= tmp->size.i;
        pop->completed.i += tmp->size.i;
      } else {
        pop->size.d -= tmp->size.d;
        pop->completed.d += tmp->size.d;
      }
      pop->last = tmp->prev;
      if (pop->last) pop->last->next = 0;
      pop->cat -= 1;
      tmp->prev = 0;
      tmp->next = 0;
      free(tmp);
      tmp = 0;
    }
    return 0;
  }
}

char accp_iterate_deterministic(accp pop,
                                double gamma_k,
                                double gamma_theta) {
    unsigned int count = 0,
                 devc = 0;
    double item = 0.0;
    chain tmp;
    pop->completed.d = 0.0;
    for (tmp=pop->first, count=0; tmp; tmp=tmp->next, count++) {
      item = tmp->size.d * (1.0-(pop->cfun)(pop->cat-1-count,gamma_theta));
      pop->completed.d += item;
    }
    //
    double *devp = (double *)calloc(pop->cat,sizeof(double));
    if (!devp) return 1;
    for (tmp=pop->first, count=0; tmp; tmp=tmp->next, count++)
      for (devc=count; devc < pop->cat; devc++)
        devp[devc] += tmp->size.d * ((pop->cfun)(devc-count,gamma_theta)-(devc==count ? 0.0 : (pop->cfun)(devc-count-1,gamma_theta)));
    //
    pop->size.d = 0.0;
    for (tmp=pop->first, count=0; tmp; tmp=tmp->next, count++) {
      tmp->size.d = devp[count];
      pop->size.d += devp[count];
    }
    //
    free(devp);
    return 0;
}

char accp_iterate_stochastic(accp pop,
                             double gamma_k,
                             double gamma_theta) {
    unsigned int item = 0,
                 count = 0;
    chain tmp;
    pop->completed.i = 0;
    for (tmp=pop->first, count=0; tmp; tmp=tmp->next, count++) {
      item = gsl_ran_binomial(RANDOM,
                              1.0-(pop->cfun)(pop->cat-1-count,gamma_theta),
                              tmp->size.i);
      tmp->size.i -= item;
      pop->completed.i += item;
    }
    //
    unsigned int devc = 0,
                 i = 0;
    unsigned int *devp = (unsigned int *)calloc(pop->cat,sizeof(unsigned int));
    if (!devp) return 1;
    for (tmp=pop->first, count=0; tmp; tmp=tmp->next, count++) {
      if (!(tmp->size.i) || count >= pop->cat) continue;
      unsigned int rsize = pop->cat - count;
      unsigned int *ret = (unsigned int *)calloc(rsize,sizeof(unsigned int));
      if (!ret) { free(devp); return 1; }
      double *p = (double *)calloc(rsize,sizeof(double));
      if (!p) { free(ret); free(devp); return 1; }
      //
      for (i=0; i<rsize; i++)
        p[i] = (pop->cfun)(i,gamma_theta);
      if (choose(tmp->size.i,p,rsize,ret)) { free(p); free(ret); free(devp); return 1; }
      for (i=0,devc=count; devc < pop->cat; i++,devc++)
        devp[devc] += ret[i];
      //
      free(p);
      free(ret);
    }
    //
    pop->size.i = 0;
    for (tmp=pop->first, devc=0; tmp; tmp=tmp->next, devc++) {
      tmp->size.i = devp[devc];
      pop->size.i += devp[devc];
    }
    //
    free(devp);
    return 0;
}

char accp_iterate(accp pop,
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
      double m = gamma_k * gamma_theta;
      double s = sqrt(gamma_theta * m);
      printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n",m,s);
    }
  break;
  case MODE_ACCP_GAMMA:
    gamma_theta = dev_sd * dev_sd / dev_mean;
    gamma_k = dev_mean / gamma_theta;
    // printf("k=%g, theta=%g\n",gamma_k,gamma_theta);
    if (gamma_k != round(gamma_k)) {
      gamma_k = round(gamma_k);
      double m = gamma_k * gamma_theta;
      double s = sqrt(gamma_theta * m);
      printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n",m,s);
    }
  break;
  case MODE_ACCP_FIXED:
    gamma_k = round(dev_mean);
    gamma_theta = 1.0 / gamma_k; // This is not used.
  break;
  case MODE_ACCP_PASCAL:
    gamma_theta = dev_mean / (dev_sd * dev_sd);
    if (gamma_theta >= 1.0 || gamma_theta == 0.0) {
      printf("The negative binomial cannot yield mean=%g and sd=%g\n",dev_mean,dev_sd);
      return 0;
    }
    gamma_k = dev_mean * gamma_theta / (1.0 - gamma_theta);
    //printf("k=%g, theta=%g\n",gamma_k,gamma_theta);
    if (gamma_k != round(gamma_k)) {
      gamma_k = round(gamma_k);
      double m = gamma_k * (1.0-gamma_theta) / gamma_theta;
      double s = sqrt(dev_mean / gamma_theta);
      printf("For the MODE_ACCP_ERLANG distribution, the shape parameter will be adjusted to yield\nMean = %g, St.dev. = %g\n",m,s);
    }
    // printf("p=%g, r=%g\n",gamma_theta,gamma_k);
  break;
  default:
    printf("I don't know what to do with this: %d\n",pop->pdist);
    return 0;
  break;
  }
  //
  if (pop->stochastic) {
    pop->completed.i = 0;
  } else {
    pop->completed.d = 0.0;
  }
  if (gamma_k != pop->cat)
    chain_resize(pop,gamma_k);
  //
  if (pop->stochastic)
    return accp_iterate_stochastic(pop, gamma_k, gamma_theta);
  else
    return accp_iterate_deterministic(pop, gamma_k, gamma_theta);
  //
  return 0;
}

char accp_survive(accp pop,
                  double prob,
                  sdnum *ret) {
    ret->i = 0;
    ret->d = 0.0;
    pop->completed.i = 0;
    pop->completed.d = 0.0;
    sdnum item;
    chain tmp;
    for (tmp=pop->first; tmp; tmp=tmp->next) {
        if (pop->stochastic) {
            if (!(tmp->size.i)) continue;
            item.i = gsl_ran_binomial(RANDOM,
                                      prob,
                                      tmp->size.i);
            tmp->size.i -= item.i;
            pop->size.i -= item.i;
            ret->i += item.i;
        } else {
            if (!(tmp->size.d)) continue;
            item.d = tmp->size.d * prob;
            tmp->size.d -= item.d;
            pop->size.d -= item.d;
            ret->d += item.d;
        }
    }
    return 0;
}
