#ifndef CBALLS_SCALAR_MOMENTS_H
#define CBALLS_SCALAR_MOMENTS_H

static inline bool cballs_raw_legacy_multipoles(const struct cmdline_data *cmd)
{
    return cballs_opt_no_normalize_histzeta(cmd)
        && (!strcmp(cmd->searchMethod, "kdtree-omp")
            || !strcmp(cmd->searchMethod, "balltree-omp")
            || !strcmp(cmd->searchMethod, "balltree-mpi"));
}

/* Subtract each neighbor's second moment before the pivot outer product.
 * For an aggregate, field2 is a sum of squares, not square(sum). */
static inline void cballs_accumulate_raw_moments(struct cmdline_data *cmd,
        bodyptr pivot, gdhistptr_sincos_omp hist, int bin,
        real field, real field2, real cosine, real sine)
{
#ifdef TPCF
    const real pivot_field = Weight(pivot)*Kappa(pivot);
    real cm = 1.0, sm = 0.0;
    for (int m = 1; m <= cmd->mChebyshev+1; m++) {
        hist->histXithreadcos[m][bin] += field*cm;
        hist->histXithreadsin[m][bin] += field*sm;
        hist->histZetaMthreadcos[m][bin][bin] -= pivot_field*field2*cm*cm;
        hist->histZetaMthreadsin[m][bin][bin] -= pivot_field*field2*sm*sm;
        hist->histZetaMthreadsincos[m][bin][bin] -= pivot_field*field2*sm*cm;
        hist->histZetaMthreadcossin[m][bin][bin] -= pivot_field*field2*cm*sm;
        const real next = cm*cosine-sm*sine;
        sm = sm*cosine+cm*sine;
        cm = next;
    }
#else
    (void)cmd; (void)pivot; (void)hist; (void)bin;
    (void)field; (void)field2; (void)cosine; (void)sine;
#endif
}
#endif
