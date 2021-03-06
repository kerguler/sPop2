#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "spop2.h"

// Memory monitoring
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
// ---

double DPOP_EPS      = DBL_MIN;
double DPOP_MAX_DAYS = 1000000;

gsl_rng *RAND_GSL;

void set_DPOP_EPS(double eps) {
    DPOP_EPS = eps;
}

void set_DPOP_MAX_DAYS(unsigned int days) {
    DPOP_MAX_DAYS = days;
}

/*
 * Initialisation routine for an age-structured sPop2 population
 *  stochastic:   logical indicator for a stochastic or a deterministic setup
 *                0: deterministic
 *                1: stochastic
 *  gamma_mode:   selection of probability distributions/methods for the method of hazards
 *                MODE_GAMMA_RAW
 *                MODE_GAMMA_HASH
 *                MODE_NBINOM_RAW
 *                MODE_ACCP_ERLANG (MODE_GAMMA_HASH)
 *                MODE_ACCP_FIXED  (MODE_GAMMA_HASH)
 *                MODE_ACCP_PASCAL (MODE_NBINOM_RAW)
 *                MODE_ACCP_GAMMA  (MODE_GAMMA_HASH) *under construction
 */
spop spop_init(unsigned char stochastic, unsigned char gamma_mode) {
  spop pop = (spop)malloc(sizeof(struct population_st));
  pop->individuals = 0;
  pop->ncat = 0;
  pop->cat = 0;
  if (stochastic) {
    RAND_GSL = get_RAND_GSL();
    pop->size.i = 0;
    pop->dead.i = 0;
    pop->developed.i = 0;
  } else {
    pop->size.d = 0.0;
    pop->dead.d = 0.0;
    pop->developed.d = 0.0;
  }
  pop->devtable = 0;
  pop->gamma_mode = gamma_mode;
  pop->stochastic = stochastic;
  switch (pop->gamma_mode) {
      case MODE_ACCP_ERLANG:
      case MODE_ACCP_FIXED:
      case MODE_ACCP_PASCAL:
      case MODE_ACCP_GAMMA:
      case MODE_ACCP_CASWELL:
          printf("Error: Please use the spop2 class for the method of accumulation\n");
          return 0;
          break;
      default:
          break;
  }
  return pop;
}

void spop_empty(spop s) {
  if (s->individuals) {
    free(s->individuals);
    s->individuals = 0;
  }
  s->ncat = 0;
  s->cat = 0;
  if (s->stochastic) {
    s->size.i = 0;
    s->dead.i = 0;
    s->developed.i = 0;
  } else {
    s->size.d = 0;
    s->dead.d = 0;
    s->developed.d = 0;
  }
  if (s->devtable)
    spop_destroy((spop *)(&(s->devtable)));
}

void spop_destroy(spop *s) {
  if (*s) {
    spop_empty(*s);
    free((*s));
    (*s) = 0;
  }
}

void spop_retrieve(spop s, char devtable, double *dev, double *size, unsigned int *limit) {
    if (!s) return;
    unsigned int count;
    individual_data *table = devtable ? ((spop) (s->devtable))->individuals : s->individuals;
    for (count = 0; count < s->cat; count++) {
        dev[count] = table[count].development;
        size[count] = s->stochastic ? (double)(table[count].number.i) : table[count].number.d;
    }
    *limit = count;
}

void spop_print(spop s) {
    if (!s) return;
    unsigned int count;
    printf("/------------------>\nGamma mode: ");
    PRINT_MODE(s->gamma_mode);
    printf("\nPopulation type: ");
    if (s->stochastic) {
        printf("Stochastic\nPopulation size: %d\nDead: %d\nDeveloped: %d\nAge\t#Cycle\tDev.\tNumber\n",
               s->size.i,
               s->dead.i,
               s->developed.i);
        for (count = 0; count < s->cat; count++) {
            printf("%d\t%d\t%d\t%d\n",
                   s->individuals[count].age,
                   s->individuals[count].devcycle,
                   s->individuals[count].development,
                   s->individuals[count].number.i);
        }
    } else {
        printf("Deterministic\nPopulation size: %g\nDead: %g\nDeveloped: %g\nAge\t#Cycle\tDev.\tNumber\n",
               s->size.d,
               s->dead.d,
               s->developed.d);
        for (count = 0; count < s->cat; count++) {
            printf("%d\t%d\t%d\t%g\n",
                   s->individuals[count].age,
                   s->individuals[count].devcycle,
                   s->individuals[count].development,
                   s->individuals[count].number.d);
        }
    }
    printf("Number of categories: %d\nOccupied categories: %d\n\\------------------>\n", s->ncat, s->cat);
}

void spop_print_to_csv(spop s) {
  if (!s) return;
  unsigned int count;
  if (s->stochastic) {
    printf("age,dev_cycle,dev,count\n");
    for (count=0; count < s->cat; count++)
      printf("%d,%d,%d,%d\n",
             s->individuals[count].age,
             s->individuals[count].devcycle,
             s->individuals[count].development,
             s->individuals[count].number.i);
  } else {
    printf("age,dev_cycle,dev,count\n");
    for (count=0; count < s->cat; count++)
      printf("%d,%d,%d,%g\n",
             s->individuals[count].age,
             s->individuals[count].devcycle,
             s->individuals[count].development,
             s->individuals[count].number.d);
  }
}

void swap(spop s, individual_data *from, individual_data *to) {
  individual_data tmp;
  tmp.age = (*to).age; tmp.devcycle = (*to).devcycle; tmp.development = (*to).development;
  (*to).age = (*from).age; (*to).devcycle = (*from).devcycle; (*to).development = (*from).development;
  (*from).age = tmp.age; (*from).devcycle = tmp.devcycle; (*from).development = tmp.development;
  if (s->stochastic) {
    tmp.number.i = (*to).number.i;
    (*to).number.i = (*from).number.i;
    (*from).number.i = tmp.number.i;
  } else {
    tmp.number.d = (*to).number.d;
    (*to).number.d = (*from).number.d;
    (*from).number.d = tmp.number.d;
  }
}

int cmpfunc(const void *a, const void *b) {
  if (is_smaller(*(individual_data*)a,*(individual_data*)b)) return -1;
  if (is_greater(*(individual_data*)a,*(individual_data*)b)) return 1;
  return 0;
}

int search(individual_data *arr, int l, int r, individual_data *x) {
  if (r >= l) {
    int mid = l + (r - l)/2;
    if (is_equal(arr[mid],(*x))) return mid;
    if (is_greater(arr[mid],(*x))) return search(arr, l, mid-1, x);
    return search(arr, mid+1, r, x);
  }
  return -1;
}

void spop_sdadd(spop s,
                unsigned int age,
                unsigned int devcycle,
                unsigned int development,
                sdnum number) {
    if (s->stochastic) {
        if (number.i <= 0) return;
    } else {
        if (number.d <= DPOP_EPS) return;
    }
    //
    individual_data tmp;
    tmp.age = age;
    tmp.devcycle = devcycle;
    tmp.development = development;
    // WARNING: This should be enough for the "search" algorithm
    //
    int cat = -1;
    if (s->individuals && s->ncat && s->cat)
        cat = search(s->individuals, 0, s->cat - 1, &tmp);
    if (cat == -1) { // Not found
        if (s->cat >= s->ncat) { // Expand buffer
            s->ncat = 2 * s->ncat + 1;
            if (!(s->individuals)) {
                s->individuals = (individual_data *) malloc(s->ncat * sizeof(individual_data));
            } else {
                s->individuals = (individual_data *) realloc(s->individuals, s->ncat * sizeof(individual_data));
            }
            unsigned int i;
            for (i = s->cat; i < s->ncat; i++) {
                s->individuals[i].age = 0;
                s->individuals[i].devcycle = 0;
                s->individuals[i].development = 0;
                if (s->stochastic)
                    s->individuals[i].number.i = 0;
                else
                    s->individuals[i].number.d = 0.0;
            }
        }
        cat = s->cat;
        s->cat += 1;
    }
    //
    s->individuals[cat].age = tmp.age;
    s->individuals[cat].devcycle = tmp.devcycle;
    s->individuals[cat].development = tmp.development;
    //
    if (s->stochastic) {
        s->individuals[cat].number.i += number.i;
        s->size.i += number.i;
    } else {
        s->individuals[cat].number.d += number.d;
        s->size.d += number.d;
    }
    //
    qsort(s->individuals, s->cat, sizeof(individual_data), cmpfunc);
}

void spop_popadd(spop s, spop d) {
  unsigned int i;
  if (!d || !(d->individuals)) return;
  for (i=0; i<d->cat; i++)
    spop_sdadd(s,
               d->individuals[i].age,
               d->individuals[i].devcycle,
               d->individuals[i].development,
               d->individuals[i].number);
}

char calc_spop(spop s, individual_data *tmpn, double prob, sdnum *k) {
  if (s->stochastic) {
    (*k).i = gsl_ran_binomial(RAND_GSL,
                              prob,
                              tmpn->number.i);
    tmpn->number.i -= (*k).i;
    s->size.i -= (*k).i;
    if (tmpn->number.i == 0)
      return 1;
  } else {
    (*k).d = tmpn->number.d * prob;
    tmpn->number.d -= (*k).d;
    s->size.d -= (*k).d;
    if (tmpn->number.d <= DPOP_EPS || s->size.d <= DPOP_EPS) {
      tmpn->number.d = 0;
      if (s->size.d <= DPOP_EPS)
        s->size.d = 0;
      return 1;
    }
  }
  return 0;
}

double calc_prob_nbinom_raw(unsigned int age, double prob, double mean, double sd) {
  return nbinom_dist_prob(mean,sd,age);
}

double calc_prob_gamma_raw(unsigned int age, double prob, double mean, double sd) {
  return gamma_dist_prob(mean,sd,age);
}

double calc_prob_gamma_hash(unsigned int age, double prob, double mean, double sd) {
  if (!gamma_dist_hash(mean,sd,age,&prob)) {
    printf("ERROR: Gamma distribution failed!\n");
    return NAN;
  }
  return prob;
}

double calc_prob_fixed(unsigned int age, double prob, double mean, double sd) {
  return age >= mean - 1.0;
}

double calc_prob_raw(unsigned int age, double prob, double mean, double sd) {
  return prob;
}

prob_func assign_prob_func(spop s, double prob, double mean, double sd) {
  if (mean == 0 && sd == 0 && prob >= 0 && prob <= 1) {
    // printf("Assigning calc_prob_raw\n");
    return calc_prob_raw;
  }
  // For a wider range and an alternative implementation:
  // http://people.math.sc.edu/Burkardt/c_src/prob/prob.html
  if (mean >= 0 && sd > 0) {
    switch (s->gamma_mode) {
    case MODE_GAMMA_RAW:
      // printf("Assigning calc_prob_gamma_raw\n");
      return calc_prob_gamma_raw;
      break;
    case MODE_ACCP_GAMMA:
    case MODE_ACCP_FIXED:
    case MODE_ACCP_ERLANG:
    case MODE_GAMMA_HASH:
      // printf("Assigning calc_prob_gamma_hash\n");
      return calc_prob_gamma_hash;
      break;
    case MODE_ACCP_PASCAL:
    case MODE_ACCP_CASWELL:
    case MODE_NBINOM_RAW:
      // printf("Assigning calc_prob_gamma_hash\n");
      return calc_prob_nbinom_raw;
      break;
    default:
      printf("ERROR: Wrong distribution option selected: %d\n",s->gamma_mode);
      return 0;
      break;
    }
  }
  //
  if (mean > 0 && sd == 0) {
    // printf("Assigning calc_prob_fixed\n");
    return calc_prob_fixed;
  }
  //
  printf("ERROR: Wrong options for probability function: prob = %g, mean = %g,  sd = %g\n",prob,mean,sd);
  return 0;
}

void spop_removeitem(spop s, unsigned int *i, char *remove) {
  unsigned int j;
  switch ((*remove)) {
  case 0:
    break;
  default:
  case -1:
    printf("ERROR: Error in population removal!\n");
    spop_print(s);
    s->dead.d = NAN;
    s->developed.d = NAN;
    (*i) = NAN;
    return;
    break;
  case 1:
    for (j=(*i); j<s->cat-1; j++)
      swap(s, &(s->individuals[j]), &(s->individuals[j+1]));
    s->cat -= 1;
    (*i)--;
    //
    if (s->ncat > 1 && s->cat < (s->ncat >> 1)) { // Shrink buffer
        unsigned int tmp = (s->ncat >> 1) - 1;
        s->ncat = tmp;
        s->individuals = (individual_data *)realloc(s->individuals, s->ncat * sizeof(individual_data));
    }
    //
    if ((s->stochastic && s->size.i == 0) || (!(s->stochastic) && s->size.d <= DPOP_EPS))
      break;
    (*remove) = 0;
    break;
  }
}

char spop_iterate(spop  s,
                  double dev_prob,       // fixed development probability
                  double dev_mean,       // mean development time
                  double dev_sd,         // sd development time
                  iter_func dev_fun,     //
                  double death_prob,     // fixed death probability
                  double death_mean,     // mean time of death
                  double death_sd,       // sd time of death
                  iter_func death_fun,   //
                  unsigned char pause) { // pause aging and development for the iteration
    if (s->stochastic) {
        s->dead.i = 0;
        s->developed.i = 0;
    } else {
        s->dead.d = 0.0;
        s->developed.d = 0.0;
    }
    if (s->devtable) {
        spop_empty((spop) (s->devtable));
    } else {
        s->devtable = (void *) spop_init(s->stochastic, s->gamma_mode);
    }
    //
    prob_func calc_prob_death;
    prob_func calc_prob_dev;
    //
    char remove = 0;
    sdnum k;
    double prob;
    individual_data *tmpn;
    unsigned int i;
    for (i = 0; i < s->cat; i++) {
        tmpn = &(s->individuals[i]);
        //
        if (death_fun)
            death_fun(tmpn, &death_prob, &death_mean, &death_sd);
        calc_prob_death = assign_prob_func(s, death_prob, death_mean, death_sd);
        if (calc_prob_death) {
            // Death
            prob = calc_prob_death(tmpn->age, death_prob, death_mean, death_sd);
            remove = calc_spop(s, tmpn, prob, &k);
            if (s->stochastic)
                s->dead.i += k.i;
            else
                s->dead.d += k.d;
        }
        //
        if (remove == 0) {
            if (dev_fun)
                dev_fun(tmpn, &dev_prob, &dev_mean, &dev_sd);
            calc_prob_dev = assign_prob_func(s, dev_prob, dev_mean, dev_sd);
            if (calc_prob_dev) {
                // Develop
                prob = calc_prob_dev(tmpn->development, dev_prob, dev_mean, dev_sd);
                remove = calc_spop(s, tmpn, prob, &k);
                if (s->stochastic)
                    s->developed.i += k.i;
                else
                    s->developed.d += k.d;
            }
            //
            if (!pause) {
                tmpn->age += 1;
                tmpn->development += 1;
            }
            //
            if (tmpn->age >= DPOP_MAX_DAYS || tmpn->development >= DPOP_MAX_DAYS) {
                printf("ERROR: Development time is too high!\n");
                spop_print(s);
                s->developed.d = NAN;
                return 1;
            }
            //
            if ((s->stochastic && k.i > 0) || (!(s->stochastic) && k.d >= DPOP_EPS)) {
                spop_sdadd((spop) (s->devtable),
                           tmpn->age,
                           pause ? tmpn->devcycle : tmpn->devcycle + 1,
                           pause ? tmpn->development : 0,
                           k);
            }
        }
        //
        spop_removeitem(s, &i, &remove);
    }
    //
    if (s->gamma_mode == MODE_GAMMA_HASH)
        gamma_dist_check();
    //
    return 0;
}

