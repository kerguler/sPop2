## Infectious disease dynamics (SEIR)

This is an example of a typical susceptible-exposed-infectious-recovered disease dynamics.

**C code snippet**

The following parameters are taken from Anderson, May, Anderson (1992) Table 1.3 on page 10.
```C
double N = 100.0;                     // Total population size (constant)
double sigma_m = 7.5, sigma_s = 7.5;  // Duration of the latent period (mean, std. dev.)
double gamma_m = 6.5, gamma_s = 6.5;  // Duration of the infectious period (mean, std. dev.)
double mu = 0.01;                     // Death and reproduction rate
double beta = 1.0;                    // Transmission coefficient
```


*ex4_SEIR.c*

**Visualisation in R**

```R
d<-read.csv("out.csv",header=F)
plot(c(0,300),c(0,30),t="n",xlab="Time (days)",ylab="Number of infectives",frame=FALSE)
xr <- d[,1]==0; lines(d[xr,3],d[xr,6],col="blue",lwd=4)
xr <- (d[,1]==1) & (d[,2]==1); lines(d[xr,3],d[xr,6],col="red",lwd=3)
xr <- (d[,1]==1) & (d[,2]==0.5); lines(d[xr,3],d[xr,6],col="black",lwd=2)
legend("topright",c("ODE","sPop - Exponential","sPop - Gamma"),col=c("blue","red","black"),lwd=c(4,3,2))
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex4_SEIR ex4_SEIR.c
```
