# sPop2: an age-structured population dynamics model

This is the standalone C library of the age structured population dynamics model sPop2.
This version implements both hazard-based and accumulative processes in development and survival.

## Installation

**Configuration**

Constructing configuration and Makefiles

```bash
./configure
```

**Compiling / linking**

Compile and install the library to the OS defaults (`/usr/local/lib` and `/usr/local/include`)

```bash
sudo make install
```

## Using the library

Including the created library with your project.

*NOTE: you have to install the library before headers can be found, and linking can be done*

**C code snippet**

*This example can be found in `examples/`.*

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

## Removing Clutter

You have to become root to cleanup, because `sudo make install` generates `src/.libs/*`.

```bash
sudo make clean
```

## Uninstalling Library

```bash
sudo make uninstall
```
