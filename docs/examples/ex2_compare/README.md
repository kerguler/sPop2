## Comparing sPop population dynamics

Age-dependent and accumulative development can be modelled using the *spop* and *spop2* data structures, respectively. 

* Method of hazards
    * *MODE_GAMMA_HASH*: Gamma distribution
    * *MODE_NBINOM_RAW*: Negative binomial distribution
* Method of accumulation
    * *MODE_ACCP_FIXED*: Fixed duration
    * *MODE_ACCP_ERLANG*: Erlang distribution
    * *MODE_ACCP_PASCAL*: Pascal distribution

**C code snippet**

*ex2_compare.c*

The following code creates an age-structured *spop* model and iterates 100 individuals for 50 time steps. Only 90% of the individuals survive each day, and they develop according to a probability distribution (given by *mode*) with a mean of 20 and a standard deviation of 10 steps.
```c
pop = spop_init(0, mode);
spop_add(pop, 0, 0, 0, 100);
for (tm=1; tm<50; tm++) {
    spop_iterate(pop,  
                 0, 20, 10, 0,  
                 0.1, 0, 0, 0,  
                 0);
    printf("%d,%g\n", tm, pop->size.d);
}
spop_destroy(&pop);
```
This is the equivalent population with the method of hazards.
```c
pop2 = spop2_init(0, mode);
spop2_add(pop2, 0, 100);
for (tm = 1; tm < 50; tm++) {
    spop2_iterate(pop2, 
                  20,     // Mean development time
                  10,     // St. dev. of development time
                  0.1,    // Daily mortality
                  0);     // (Logical) Mortality applied to devtable instead
    printf("%d,%g\n", tm, pop2->size.d);
}
spop2_destroy(&pop2);
```

**Compile and run**

```bash
$ gcc -o ex2_compare ex2_compare.c -Wall -lm -lspop2 -lgsl
```

**Note**

Please note that the shift of 1 time unit between *MODE_NBINOM_RAW* and *MODE_ACCP_PASCAL* is fixed in sPop2.
