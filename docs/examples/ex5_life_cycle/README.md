## Insect life cycle

This is an example of a typical insect life cycle modelled with the spop2 data structure (model of accumulation).

**C code snippet**

*ex5_life_cycle.c*

Declare a distinct pseudo-stage-structured population for each of the four life stages.
```c
unsigned char mode = MODE_ACCP_ERLANG;
spop2 egg = spop2_init(0, mode);
spop2 larva = spop2_init(0, mode);
spop2 pupa = spop2_init(0, mode);
spop2 adult = spop2_init(0, mode);
```
Introduce 100 individuals (no history of prior development) to the population of eggs.
```c
spop2_add(egg, 0, 100);
```
At each time step, the following processes should be executed in the given order.
1. Apply the **survival** and **development** events
2. Perform the **stage transformations** and 
```c
for (tm=1; tm<100; tm++) {
    // The following lines will apply the survival (first) and development (next) 
    // processes to the respective life stages.
    spop2_iterate(egg,  10, 2, 0, 0); // 10+-2 steps of development, no death
    spop2_iterate(larva,10, 2, 0, 0); // 10+-2 steps of development, no death
    spop2_iterate(pupa, 10, 2, 0, 0); // 10+-2 steps of development, no death
    spop2_iterate(adult,10, 2, 0, 0); // 10+-2 steps of development (as lifetime)
    //
    // The following lines perform appropriate stage transformations.
    spop2_add(egg,  0, 1.0 * adult->size.d);     // Adult females lay 1 egg per day
    spop2_add(larva,0, egg->developed.d);        // Eggs completing development will hatch into larvae
    spop2_add(pupa, 0, larva->developed.d);      // Larvae completing development will hatch into pupae
    spop2_add(adult,0, 0.5 * pupa->developed.d); // Half the pupae will become adult females
}
```
Release the memory used by the populations.
```c
spop2_destroy(&egg);
spop2_destroy(&larva);
spop2_destroy(&pupa);
spop2_destroy(&adult);
```

**Compile and run**

```bash
$ gcc -o ex5_life_cycle ex5_life_cycle.c -Wall -lm -lspop2 -lgsl
```
