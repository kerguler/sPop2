## Insect life cycle

This is an example of a typical insect life cycle modelled with sPop.

**C code snippet**

*ex5_life_cycle.c*

Declare a distinct population for each of the four life stages.
```c
unsigned char mode = MODE_GAMMA_HASH;
spop egg = spop_init(0, mode);
spop larva = spop_init(0, mode);
spop pupa = spop_init(0, mode);
spop adult = spop_init(0, mode);
```
Introduce 100 individuals (no history of prior development) to the population of eggs.
```c
spop_add(egg, 0, 0, 0, 0, 100);
```
At each time step, the following processes should be executed in the given order.
1. Apply the **survival** and then the **development** events
2. Perform the **stage transformations** and 
```c
for (tm=1; tm<100; tm++) {
    // The following lines will apply the survival (first) and development (next) 
    // processes to the respective life stages.
    spop_iterate(egg,  0, 10, 1, 0,  0, 0, 0, 0,  0); // 20+-1 steps of development, no death
    spop_iterate(larva,  0, 10, 1, 0,  0, 0, 0, 0,  0); // 20+-1 steps of development, no death
    spop_iterate(pupa,  0, 10, 1, 0,  0, 0, 0, 0,  0); // 20+-1 steps of development, no death
    spop_iterate(adult,  0, 0, 0, 0,  0, 10, 1, 0,  0); // no development, 10+-1 steps of lifetime
    //
    // The following lines perform appropriate stage transformations.
    spop_add(egg, 0, 0, 0, 0, 1.0 * adult->size.d); // Adult females lay 1 egg per day
    spop_add(larva, 0, 0, 0, 0, egg->developed.d); // Eggs completing development will hatch into larvae
    spop_add(pupa, 0, 0, 0, 0, larva->developed.d); // Larvae completing development will hatch into pupae
    spop_add(adult, 0, 0, 0, 0, 0.5 * pupa->developed.d); // Half the pupae will become adult females
}
```
Release the memory used by the populations.
```c
spop_destroy(&egg);
spop_destroy(&larva);
spop_destroy(&pupa);
spop_destroy(&adult);
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex5_life_cycle ex5_life_cycle.c
```
