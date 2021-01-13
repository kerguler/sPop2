## Infectious disease dynamics (SEIR)

This is an example of a typical susceptible-exposed-infectious-recovered disease dynamics.

**C code snippet**

*ex4_SEIR.c*

**Visualisation in R**

```R
d<-read.csv("out.csv",header=F)
plot(d[,1],d[,2],t="l",ylim=c(0,100))
lines(d[,1],d[,3])
lines(d[,1],d[,4])
lines(d[,1],d[,5])
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex4_SEIR ex4_SEIR.c
```
