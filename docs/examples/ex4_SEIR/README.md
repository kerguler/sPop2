## Infectious disease dynamics (SEIR)

This is an example of a typical susceptible-exposed-infectious-recovered disease dynamics.

**C code snippet**

*ex4_SEIR.c*

**Visualisation in R**

```R
d<-read.csv("out.csv",header=F)
plot(d[,2],d[,3],t="l",ylim=c(0,100))
lines(d[,2],d[,4])
lines(d[,2],d[,5])
lines(d[,2],d[,6])

d1<-read.csv("out.csv",header=F)
plot(d1[,2],d1[,3],t="l",xlim=c(0,100),ylim=c(0,100))
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex4_SEIR ex4_SEIR.c
```
