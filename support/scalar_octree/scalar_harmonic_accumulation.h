#ifndef CBALLS_SCALAR_HARMONIC_ACCUMULATION_H
#define CBALLS_SCALAR_HARMONIC_ACCUMULATION_H
/* Private numerical implementation: include after the engine types.
 * Arithmetic and iteration order are preserved; see engineering_contracts.rst. */
#ifdef THREEPCFCONVERGENCE
/* Keep the product/subtraction/weight order, but avoid eight temporary matrices
 * and the separate clear, scale, and add passes for every pivot and mode. */
local inline void ggg_accumulate_outer_products(
    int bins, real xi, const real *c, const real *s,
    const real *diag_c, const real *diag_s, const real *diag_sc,
    real **cc, real **ss, real **sc, real **cs)
{
    for (int i = 1; i <= bins; i++) {
        real *restrict out_cc = cc[i];
        real *restrict out_ss = ss[i];
        real *restrict out_sc = sc[i];
        real *restrict out_cs = cs[i];
        for (int j = 1; j <= bins; j++) {
            real cc_value = c[i]*c[j];
            real ss_value = s[i]*s[j];
            real sc_value = s[i]*c[j];
            real cs_value = c[i]*s[j];
            if (i == j) {
                cc_value -= diag_c[i];
                ss_value -= diag_s[i];
                sc_value -= diag_sc[i];
                cs_value -= diag_sc[i];
            }
            out_cc[j] += cc_value*xi;
            out_ss[j] += ss_value*xi;
            out_sc[j] += sc_value*xi;
            out_cs[j] += cs_value*xi;
        }
    }
}
#endif

#if defined(NMultipoles) && defined(THREEPCFCONVERGENCE)
/* One angular recurrence feeds both the signal and window, including the
 * self-pair terms. Edge coupling needs window modes through 2*mChebyshev. */
local inline void ggg_accumulate_modes(
    struct cmdline_data *cmd, gdhistptr_sincos_omp_ggg hist,
    gdhistptr_sincos_omp_ggg_N histN, int n,
    real cosphi, real sinphi, real xi, real xiN, real field2, real weight2)
{
    const int signal_orders = cmd->mChebyshev + 1;
    const int orders = histN != NULL
        ? GGG_WINDOW_ORDERS(cmd) : signal_orders;
    real t_previous = 1.0, t_before = 0.0;
    real u_previous = 0.0, u_before = 0.0;

    for (int m = 1; m <= orders; m++) {
        real t, u;
        if (m == 1) {
            t = 1.0;
            u = 0.0;
        } else if (m == 2) {
            t = cosphi;
            u = 1.0;
        } else {
            t = 2.0*cosphi*t_previous - t_before;
            u = 2.0*cosphi*u_previous - u_before;
        }
        const real sine = u*sinphi;
        if (histN != NULL) {
            const real nc = xiN*t;
            const real ns = xiN*sine;
            histN->histXithreadcos[m][n] += nc;
            histN->histXithreadsin[m][n] += ns;
            histN->histXithreaddiagcos[m][n] += weight2*t*t;
            histN->histXithreaddiagsin[m][n] += weight2*sine*sine;
            histN->histXithreaddiagsincos[m][n] += weight2*sine*t;
        }
        if (m <= signal_orders) {
            const real kc = xi*t;
            const real ks = (xi*u)*sinphi;
            hist->ChebsT[m] = t;
            hist->ChebsU[m] = u;
            hist->histXithreadcos[m][n] += kc;
            hist->histXithreadsin[m][n] += ks;
            hist->histXithreaddiagcos[m][n] += field2*t*t;
            hist->histXithreaddiagsin[m][n] += field2*sine*sine;
            hist->histXithreaddiagsincos[m][n] += field2*sine*t;
        }
        t_before = t_previous;
        t_previous = t;
        u_before = u_previous;
        u_previous = u;
    }
}
#endif

#endif
