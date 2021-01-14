#include "spop2/spop2.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

// For measles (Anderson, May, Anderson - 1992 - Table 1.3)
double N = 100.0;
double sigma_m = 7.5, sigma_s = 1.5;
double gamma_m = 6.5, gamma_s = 0.5;
double mu = 0.01;
double beta = 1.0;

void print_out(double tm, double S, spop *E, spop *I, double R) {
    printf("1,%g,%g,%g,%g,%g\n",tm,S,(*E)->size.d,(*I)->size.d,R);
}

int func (double t, const double y[], double f[], void *params) {
    (void) (t); /* avoid unused parameter warning */
    /* S */ f[0] = mu * N - mu * y[0] - (beta / N) * y[2] * y[0];
    /* E */ f[1] = (beta / N) * y[2] * y[0] - (1.0 / sigma_m) * y[1] - mu * y[1];
    /* I */ f[2] = (1.0 / sigma_m) * y[1] - (1.0 / gamma_m) * y[2] - mu * y[2];
    /* R */ f[3] = (1.0 / gamma_m) * y[2] - mu * y[3];
    return GSL_SUCCESS;
}

void sim_ode(void) {
    gsl_odeiv2_system sys = {func, 0, 4, 0};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 300.0;
    double y[4] = {N-1.0, 0.0, 1.0, 0.0};
    printf("0,%g,%g,%g,%g,%g\n", t, y[0], y[1], y[2], y[3]);

    for (i = 1; i <= 1000; i++) {
        double ti = i * t1 / 1000.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        printf("0,%g,%g,%g,%g,%g\n", t, y[0], y[1], y[2], y[3]);
    }
    gsl_odeiv2_driver_free(d);
}

void sim_spop(double tau) {
    unsigned char mode = MODE_GAMMA_HASH;
    unsigned int tm = 0;
    //
    double tau_sigma_m = 1.0-pow(1.0-(1.0/sigma_m), 1.0/tau); // sigma_m * tau;
    double tau_sigma_s = sigma_s * tau;
    double tau_gamma_m = 1.0-pow(1.0-(1.0/gamma_m), 1.0/tau); // gamma_m * tau;
    double tau_gamma_s = gamma_s * tau;
    double tau_mu = 1.0-pow(1.0-mu, 1.0/tau);
    double tau_beta = beta / tau;
    //
    double S = 0;
    spop E = spop_init(0, mode);
    spop I = spop_init(0, mode);
    double R = 0;
    //
    double v = 0;
    //
    S = N-1.0;
    spop_add(E, 0, 0, 0, 1.0);
    print_out((double)(tm)/tau, S, &E, &I, R);
    //
    for (tm=1; tm<300*tau; tm++) {
        v = tau_beta * I->size.d * S / N;
        //
        S -= tau_mu * S;
        //spop_iterate(E,  0, tau_sigma_m, tau_sigma_s, 0,  tau_mu, 0, 0, 0,  0);
        //spop_iterate(I,  0, tau_gamma_m, tau_gamma_s, 0,  tau_mu, 0, 0, 0,  0);
        spop_iterate(E,  tau_sigma_m, 0, 0, 0,  tau_mu, 0, 0, 0,  0);
        spop_iterate(I,  tau_gamma_m, 0, 0, 0,  tau_mu, 0, 0, 0,  0);
        R -= tau_mu * R;
        //
        R += I->developed.d;
        spop_add(I, 0, 0, 0, E->developed.d);
        spop_add(E, 0, 0, 0, v);
        S += tau_mu * N - v;
        //
        print_out((double)(tm)/tau, S, &E, &I, R);
    }
    //
    spop_destroy(&E);
    spop_destroy(&I);
}

int main(void) {
    sim_ode();
    sim_spop(5.0);
    return 0;
}
