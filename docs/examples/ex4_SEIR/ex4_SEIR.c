#include "spop2/spop2.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

void print_out(unsigned int tm, double S, spop *E, spop *I, double R) {
    printf("%d,%g,%g,%g,%g\n",tm,S,(*E)->size.d,(*I)->size.d,R);
}

int func (double t, const double y[], double f[], void *params) {
    (void) (t); /* avoid unused parameter warning */
    // For measles (Anderson, May, Anderson - 1992 - Table 1.3)
    double N = 100.0;
    double sigma_m = 7.5;
    double gamma_m = 6.5;
    double mu = 0.01;
    double beta = 1.0;
    //
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
    double y[4] = {99.0, 0.0, 1.0, 0.0};
    printf("%g,%g,%g,%g,%g\n", t, y[0], y[1], y[2], y[3]);

    for (i = 1; i <= 1000; i++) {
        double ti = i * t1 / 1000.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        printf("%g,%g,%g,%g,%g\n", t, y[0], y[1], y[2], y[3]);
    }
    gsl_odeiv2_driver_free(d);
}

void sim_spop(void) {
    unsigned char mode = MODE_GAMMA_HASH;
    unsigned int tm = 0;
    //
    double N = 100;
    double S = 0;
    spop E = spop_init(0, mode);
    spop I = spop_init(0, mode);
    double R = 0;
    //
    // For measles (Anderson, May, Anderson - 1992 - Table 1.3)
    double sigma_m = 7.5, sigma_s = 1.5;
    double gamma_m = 6.5, gamma_s = 0.5;
    double mu = 0.01;
    double beta = 1.0;
    //
    double v = 0;
    //
    S = 99;
    spop_add(I, 0, 0, 0, 1);
    print_out(tm, S, &E, &I, R);
    //
    for (tm=1; tm<300; tm++) {
        v = beta * I->size.d * S / N;
        //
        S -= mu * S;
        spop_iterate(E,  0, sigma_m, sigma_s, 0,  mu, 0, 0, 0,  0);
        spop_iterate(I,  0, gamma_m, gamma_s, 0,  mu, 0, 0, 0,  0);
        R -= mu * R;
        //
        R += I->developed.d;
        spop_add(I, 0, 0, 0, E->developed.d);
        spop_add(E, 0, 0, 0, v);
        S += mu * N - v;
        //
        print_out(tm, S, &E, &I, R);
    }
    //
    spop_destroy(&E);
    spop_destroy(&I);
}

int main(void) {
    sim_ode();
    //sim_spop();
    return 0;
}
