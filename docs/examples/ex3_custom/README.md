## Age-structured population dynamics

This is the hypothetical age-structured population dynamics model with a custom survival schedule as described in <a href="https://doi.org/10.12688/f1000research.15824.3" target="_blank">Erguler (2018)</a>.

**C code snippet**

*ex3_custom.c*

Custom survival function
```c
void death(const individual_data *ind, double *death_prob, double *death_mean, double *death_sd) {
    (*death_prob) = 0;
    (*death_mean) = 480.0 - (ind->devcycle > 4 ? 240.0 : 48.0 * ind->devcycle);
    (*death_sd) = 0.1 * (*death_mean);
}
```
Initiate the random number generator.
```c
rng_setup("spop.c");
```

Initiate the spop data structure for a stochastic population using the gamma distribution as the basis of survival and development.
```c
vec = spop_init(1,                  // This is a stochastic population 
                MODE_GAMMA_HASH);
```
Add 1000 individuals to the population (age=0, devcycle=0, development=0) and iterate for 720 time units (hours).
```c
spop_add(vec, 0, 0, 0, 1000);
for (i = 0; i < 720; i++) {
    spop_iterate(vec,
                 0, 50.0, 10.0, 0,
                 0, 0, 0, death,
                 0);
    /*
     * At each iteration, reintroduce to the population all individuals which completed their development
     * for the next round of development.
     */
     spop_popadd(vec, vec->devtable);
}
```
Clear population contents.
```c
spop_empty(vec);
```
Free the space occupied by the population data structure
```c
spop_destroy(&vec);
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex3_custom ex3_custom.c
```
