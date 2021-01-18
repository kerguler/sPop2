## Infectious disease dynamics (SEIR)

This is an example of a typical susceptible-exposed-infectious-recovered disease dynamics.

**C code snippet**

The following parameters are taken from Anderson, May, Anderson (1992) Table 1.3 on page 10, and represent the transmission dynamics of measles.
```c
double N = 100.0;                     // Total population size (constant)
double sigma_m = 7.5, sigma_s = 7.5;  // Duration of the latent period (mean, std. dev.)
double gamma_m = 6.5, gamma_s = 6.5;  // Duration of the infectious period (mean, std. dev.)
double mu = 0.01;                     // Death and reproduction rate
double beta = 1.0;                    // Transmission coefficient
```

*ex4_SEIR.c*

Here, we compare the ODE and sPop representations of the disease dynamics. The following function is required for the ODE model.
```c
int func (double t, const double y[], double f[], void *params) {
    (void) (t); /* avoid unused parameter warning */
    /* S */ f[0] = mu * N                   - (beta / N) * y[2] * y[0] - mu * y[0];
    /* E */ f[1] = - (1.0 / sigma_m) * y[1] + (beta / N) * y[2] * y[0] - mu * y[1];
    /* I */ f[2] =   (1.0 / sigma_m) * y[1] - (1.0 / gamma_m) * y[2]   - mu * y[2];
    /* R */ f[3] =                            (1.0 / gamma_m) * y[2]   - mu * y[3];
    return GSL_SUCCESS;
}
```
The following routine is used to numerically integrate the above equations system using the explicit embedded Runge-Kutta Prince-Dormand method implemented in GSL.
```c
gsl_odeiv2_system sys = {func, 0, 4, 0};
gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
double t = 0.0, t1 = 300.0;
double y[4] = {N-1.0, 0.0, 1.0, 0.0};
int i;
for (i = 1; i <= 1000; i++) {
    double ti = i * t1 / 1000.0;
    int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
    if (status != GSL_SUCCESS) {
        printf("error, return value=%d\n", status);
        break;
    }
}
gsl_odeiv2_driver_free(d);
```
The following routine is used to represent the same system using sPop. Here, **sd** is reserved as a scaling factor. Setting **sd=1** indicates an exponentially-distributed stage duration, which leads to a dynamics analogous to the one simulated with the above ODEs. On the other hand, **tau** is used to scale time.
```c
unsigned char mode = MODE_GAMMA_HASH;
unsigned int tm = 0;
double tau = 10.0;

double tau_sigma_m = sigma_m * tau;
double tau_sigma_s = sigma_s * sd * tau;
double tau_gamma_m = gamma_m * tau;
double tau_gamma_s = gamma_s * sd * tau;
double tau_mu = 1.0-pow(1.0-mu, 1.0/tau);
double tau_beta = beta / tau;

double S = 0;
spop E = spop_init(0, mode);
spop I = spop_init(0, mode);
double R = 0;

double v = 0;

S = N-1.0;
spop_add(I, 0, 0, 0, 1.0);
for (tm=1; tm<300*tau; tm++) {
    v = tau_beta * I->size.d * S / N;

    spop_iterate(E,  0, tau_sigma_m, tau_sigma_s, 0,  tau_mu, 0, 0, 0,  0);
    spop_iterate(I,  0, tau_gamma_m, tau_gamma_s, 0,  tau_mu, 0, 0, 0,  0);

    R = I->developed.d + (1.0 - tau_mu) * R;
    spop_add(I, 0, 0, 0, E->developed.d);
    spop_add(E, 0, 0, 0, v);
    S = tau_mu * N + (1.0 - tau_mu) * S - v;
}
spop_destroy(&E);
spop_destroy(&I);
```

**Compile and run**

```bash
$ gcc -Wall -lm -lspop2 -lgsl -o ex4_SEIR ex4_SEIR.c
```

**Visualisation in R**

```r
d<-read.csv("out.csv",header=F)
plot(c(0,300),c(0,30),t="n",xlab="Time (days)",ylab="Number of infectives",frame=FALSE)
xr <- d[,1]==0; lines(d[xr,3],d[xr,6],col="blue",lwd=4)
xr <- (d[,1]==1) & (d[,2]==1); lines(d[xr,3],d[xr,6],col="red",lwd=3)
xr <- (d[,1]==1) & (d[,2]==0.5); lines(d[xr,3],d[xr,6],col="black",lwd=2)
legend("topright",c("ODE","sPop - Exponential","sPop - Gamma"),col=c("blue","red","black"),lwd=c(4,3,2))
```
