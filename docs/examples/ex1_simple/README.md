## Creating an sPop population

The sPop2 library can be used to create a structured population.

**C code snippet**

*ex1_simple.c*

```c
#include "spop2/spop2.h"
```
Declare an spop object
```c
spop pop = spop_init(0,                  // Deterministic
                     MODE_ACCP_ERLANG);  // Method of accumulation with Erlang distribution
```
Introduce individuals to the population
```c
spop_add(pop,    // The spop object 
         0,      // Age                   (method of hazards)
         0,      // Development cycle     (method of hazards)
         0,      // Development indicator (method of hazards)
         0,      // Acc. dev. object      (method of accumulation)
         1000);  // Population size
```
Take one time step (survive and develop)
```c
spop_iterate(pop,     // The spop object
             0,       // Development (fixed daily probability)
             10,      // Development (mean number of steps)
             5,       // Development (st.dev. of the number of steps) 
             0,       // Development (function to determine daily probability)
             0.25,    // Death (fixed daily probability)
             0,       // Death (mean number of steps)
             0,       // Death (st.dev. of the number of steps) 
             0,       // Death (function to determine daily probability)
             0);      // Pause aging and development
```
Inspect the population
```c
spop_print(pop);
```
Cleanup and exit
```c
spop_destroy(&pop);
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex1_simple ex1_simple.c
```
