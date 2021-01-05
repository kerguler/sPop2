## Comparing sPop population dynamics

Different types of structured populations can be created with the sPop2 library.

**C code snippet**

*ex2_compare.c*

Compare the following probability distributions of process duration.
```c
unsigned int modes[5] = {
    MODE_GAMMA_HASH,  // Method of hazards - Gamma distribution
    MODE_NBINOM_RAW,  // Method of hazards - Negative binomial distribution
    MODE_ACCP_FIXED,  // Method of accumulation - Fixed duration
    MODE_ACCP_ERLANG, // Method of accumulation - Erlang distribution
    MODE_ACCP_PASCAL  // Method of accumulation - Pascal distribution
};
```
For each **mode**, a population is created and iterated for 50 steps with an expectation of 20 steps of development time (with a standard deviation of 10 steps).
```c
for (mode=0; mode<5; mode++) {
    pop = spop_init(0, modes[mode]);
    //
    spop_add(pop, 0, 0, 0, 0, 1000);
    printf("%d,%d,%g\n", mode, 0, pop->size.d);
    //
    for (tm=1; tm<50; tm++) {
        spop_iterate(pop,  0, 20, 10, 0,  0, 0, 0, 0,  0);
        printf("%d,%d,%g\n", mode, tm, pop->size.d);
    }
    //
    spop_destroy(&pop);
}
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex2_compare ex2_compare.c
```

**Note**

Please note that the shift of 1 time unit between *MODE_NBINOM_RAW* and *MODE_ACCP_PASCAL* is fixed in sPop2.
