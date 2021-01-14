# sPop2: a dynamically-structured matrix population model

<p align="center">
<img width="623" height=211" src="docs/figures/logo_sPop2.jpg"/>
</p>

This is the standalone C library of the dynamically-structured matrix population model sPop2.
This version implements both hazard-based and accumulative processes. While the method of hazards is implemented in the *spop* data structure, the method of accumulation is implemented in *spop2*.

## Installation

**Configuration**

Construct configuration and Makefiles

```bash
./configure
```

**Compiling / linking**

Compile and install the library to the OS defaults (`/usr/local/lib` and `/usr/local/include`)

```bash
make
sudo make install
```

## Using the library

Include the created library in your project.

**C code snippet**

Please see <a href="docs/examples/">docs/examples</a> for further documentation and usage examples.
The following is an excerpt from <a href="docs/examples/ex1_simple">ex1_simple</a>.

```c
#include "spop2/spop2.h"
// ...
spop pop = spop_init(0, MODE_GAMMA_HASH);
spop_add(pop, 0, 0, 0, 1000);
spop_iterate(pop,
             0, 10, 5, 0,
             0.25, 0, 0, 0,
             0);
spop_print(pop);
spop_destroy(&pop);
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex1_simple ex1_simple.c
```

## Removing Clutter

You have to become root to cleanup, because `sudo make install` generates `src/.libs/*`.

```bash
sudo make clean
```

## Uninstalling Library

```bash
sudo make uninstall
```

## References

* Erguler K. and Mendel J. A dynamically-structured matrix population model based on renewal processes for accumulative development under variable environmental conditions. BioRxiv. 2021
* Erguler K. sPop: Age-structured discrete-time population dynamics model in C, Python, and R [version 3; peer review: 2 approved]. F1000Research 2020, 7:1220 (https://doi.org/10.12688/f1000research.15824.3)
