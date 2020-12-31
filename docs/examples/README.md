## Using the sPop2 library

Using the sPop2 library to develop a population dynamics model.

**C code snippet**

*example_acc1_simple.c*

```c
#include "spop2/spop2.h"
spop pop = spop_init(0, MODE_ACCP_ERLANG);
spop_add(pop, 0, 0, 0, 0, 1000);
spop_iterate(pop,
             0, 10, 5, 0,
             0.25, 0, 0, 0,
             0);
spop_print(pop);
spop_destroy(&pop);
```
