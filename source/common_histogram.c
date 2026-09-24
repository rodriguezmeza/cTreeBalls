/* ==============================================================================
 MODULE: cballsutils.c			    [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "tree_contracts.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

#include <errno.h>
#include <sys/wait.h>

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdhistptr_sincos_omp hist;
} sincos_hist_init_context;

local int search_init_sincos_omp_unguarded(void *argument)
{
    sincos_hist_init_context *context = argument;
    struct cmdline_data *cmd = context->cmd;
    gdhistptr_sincos_omp hist = context->hist;
    int n;
    int m;

#ifdef TPCF
    if (!cballs_opt_only_2pcf(cmd)) {
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    } /* scalar 3PCF requested */
#endif
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
//B 2pcf
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
//E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
#ifdef TPCF
    if (!cballs_opt_only_2pcf(cmd)) {
        hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);

        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,
                      1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

        hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    } /* scalar 3PCF requested */
#endif

   for (n = 1; n <= cmd->sizeHistN; n++) {
       hist->histNthread[n] = 0.0;
       hist->histNNSubthread[n] = 0.0;
       hist->histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
       hist->histNNSubXi2pcfthreadp[n] = 0.0;
       hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
       hist->histXi2pcfthread[n] = 0.0;
       hist->histXi2pcfthreadsub[n] = 0.0;
   }

#ifdef TPCF
    if (!cballs_opt_only_2pcf(cmd)) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        }
    } /* scalar 3PCF requested */
#endif

    return SUCCESS;
}

global int search_init_sincos_omp(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  gdhistptr_sincos_omp hist)
{
    sincos_hist_init_context context;
    ErrorMsg allocation_error;

    memset(hist, 0, sizeof(*hist));
    context.cmd = cmd;
    context.gd = gd;
    context.hist = hist;
    if (cballs_allocation_guard(search_init_sincos_omp_unguarded,
                                &context, allocation_error,
                                sizeof(allocation_error)) == FAILURE) {
        search_free_sincos_omp(cmd, gd, hist);
        return FAILURE;
    }
    return SUCCESS;
}

global int search_free_sincos_omp(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
#define FREE_DVECTOR_IF_SET(p,nl,nh) \
    do { if ((p) != NULL) { free_dvector((p),(nl),(nh)); (p) = NULL; } } while (0)
#define FREE_DMATRIX_IF_SET(p,nrl,nrh,ncl,nch) \
    do { if ((p) != NULL) { free_dmatrix((p),(nrl),(nrh),(ncl),(nch)); (p) = NULL; } } while (0)
#define FREE_DMATRIX3D_IF_SET(p,nrl,nrh,ncl,nch,ndl,ndh) \
    do { if ((p) != NULL) { free_dmatrix3D((p),(nrl),(nrh),(ncl),(nch),(ndl),(ndh)); (p) = NULL; } } while (0)
#ifdef TPCF
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX_IF_SET(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadcossin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadsincos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif
    FREE_DVECTOR_IF_SET(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histXi2pcfthread,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNNSubthread,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNthread,1,cmd->sizeHistN);
#ifdef TPCF
        FREE_DVECTOR_IF_SET(hist->ChebsU,1,cmd->mChebyshev+1);
        FREE_DVECTOR_IF_SET(hist->ChebsT,1,cmd->mChebyshev+1);
#endif

#undef FREE_DVECTOR_IF_SET
#undef FREE_DMATRIX_IF_SET
#undef FREE_DMATRIX3D_IF_SET

    return SUCCESS;
}

global int computeBodyProperties_sincos(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp hist)
{
    int n;
    int m;
    real xi = 0.0;
    real xi_2p = 0.0;

    //B Normalization of histograms
    if (Type(p) == BODY) {
#ifdef NOSTANDARNORMHIST
        xi = Kappa(p);
        xi_2p = Kappa(p);
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
            xi_2p = KappaRmin(p);
            xi = NbRmin(p)*xi_2p;
    }
#endif
#else // ! NOSTANDARNORMHIST
        xi = Kappa(p)/nbody;
#ifdef BALLS4SCANLEV
        if (cballs_method_needs_balls4_scan(gd->searchMethod_int))
            xi_2p = (Weight(p)/MAX((real)Nb(p), 1.0))*Kappa(p);
        else
            xi_2p = Weight(p)*Kappa(p);
#else
        xi_2p = Weight(p)*Kappa(p);
#endif
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
#ifdef BALLS4SCANLEV
            if (cballs_method_needs_balls4_scan(gd->searchMethod_int))
                xi_2p = KappaRmin(p);
#endif
            xi = NbRmin(p)*xi_2p/nbody;
    }
#endif
#endif // ! NOSTANDARNORMHIST
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
    //E Normalization of histograms
    const bool raw_multipoles = cballs_raw_legacy_multipoles(cmd);
    if (raw_multipoles) {
        xi = Weight(p)*Kappa(p);
        xi_2p = xi;
#ifdef SMOOTHPIVOT
        if (cballs_opt_smooth_pivot(cmd)
            && !strncmp(cmd->searchMethod, "kdtree-", 7)) {
            xi = KappaRmin(p)/MAX((real)NbRmin(p), 1.0);
            xi_2p = KappaRmin(p);
        }
#endif
    }

#ifdef TPCF
    if (!cballs_opt_only_2pcf(cmd)) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            //B Normalization of histograms
#ifdef NOSTANDARNORMHIST
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= 1.0;
                hist->histXithreadsin[m][n] /= 1.0;
            }
#else
            for (n=1; n<=cmd->sizeHistN; n++) {
                if (!raw_multipoles) {
                    hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                    hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                }
            }
#endif
            //E
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos,
                      hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin,
                      hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos,
                      hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            OUTVP_ext(hist->xiOUTVPcossin,
                      hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMtmpcossin,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,
                      hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            MULMS_ext(hist->histZetaMtmpcossin,
                      hist->xiOUTVPcossin,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],
                     hist->histZetaMthreadcos[m],
                     hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],
                     hist->histZetaMthreadsin[m],
                     hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],
                     hist->histZetaMthreadsincos[m],
                     hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            ADDM_ext(hist->histZetaMthreadcossin[m],
                     hist->histZetaMthreadcossin[m],
                     hist->histZetaMtmpcossin,cmd->sizeHistN);
        }
    } /* scalar 3PCF requested */
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
    }

    return SUCCESS;
}

global int search_init_gd_hist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;
    int n1, n2, l;

#ifdef TPCF
        if (gd->common_scalar_3pcf) for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaM[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
    }

#ifdef TPCF
        if (gd->common_scalar_3pcf) for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaGmRe[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaGmIm[m], cmd->sizeHistN);
        }
#endif

    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

global int search_init_gd_hist_sincos(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

#ifdef TPCF
        if (gd->common_scalar_3pcf) for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
#ifdef TPCF
            if (gd->common_scalar_3pcf) for (m = 1; m <= cmd->mChebyshev+1; m++) {
                // HERE MUST BE gd->histXicos and gd->histXisin
            }
#endif
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

/* HistN producers supply ordered, self-excluded auto-pair counts. Publish
   unordered DD and (optionally) the existing N^2 density estimator:
       xi = 2 DD V / (N^2 shell_volume) - 1.
   Keep catalog sizes as INTEGER and convert before all count products. */
global int search_normalize_count_histograms(struct cmdline_data *cmd,
                                             struct global_data *gd,
                                             INTEGER nbody,
                                             real *histNN, real *histCF)
{
    const double count = (double)nbody;
    const int compute_cf = cballs_opt_and_cf(cmd);
    double volume = 1.0;

    if (histNN == NULL || cmd->sizeHistN < 1 || (compute_cf && histCF == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "count normalization: invalid histogram storage");
        return FAILURE;
    }
    if (compute_cf) {
        if (nbody <= 0 || !isfinite(count*count)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "count normalization: body count must be positive and representable");
            return FAILURE;
        }
        for (int k = 0; k < NDIM; k++) {
            if (!isfinite(gd->Box[k]) || gd->Box[k] <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "count normalization: box lengths must be finite and positive");
                return FAILURE;
            }
            volume *= (double)gd->Box[k];
        }
        if (!isfinite(volume) || volume <= 0.0
            || !isfinite(cmd->rminHist) || cmd->rminHist < 0.0
            || !isfinite(cmd->rangeN) || cmd->rangeN <= cmd->rminHist
            || !isfinite(gd->deltaR) || gd->deltaR <= 0.0
            || (cmd->useLogHist && cmd->rminHist == 0.0 && cmd->logHistBinsPD <= 0)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "count normalization: invalid volume or radial domain");
            return FAILURE;
        }
    }
    for (int n = 1; n <= cmd->sizeHistN; n++) {
        if (!isfinite(histNN[n]) || histNN[n] < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "count normalization: invalid pair count in bin %d", n);
            return FAILURE;
        }
        if (compute_cf) {
            double r0, r1, shell;
            if (cmd->useLogHist) {
                if (cmd->rminHist == 0.0) {
                    /* The zero-cutoff legacy bin index truncates toward zero:
                       bin 1 also accepts the interval immediately below the
                       nominal grid. Normalize the actual accepted interval. */
                    const int lower = n == 1 ? -1 : n-1;
                    r0 = cmd->rangeN * pow(10.0,
                        ((double)lower-cmd->sizeHistN)/cmd->logHistBinsPD);
                    r1 = cmd->rangeN * pow(10.0,
                        ((double)n-cmd->sizeHistN)/cmd->logHistBinsPD);
                } else {
                    r0 = cmd->rminHist * pow(10.0, (double)(n-1)*gd->deltaR);
                    r1 = cmd->rminHist * pow(10.0, (double)n*gd->deltaR);
                }
            } else {
                r0 = cmd->rminHist + (double)(n-1)*gd->deltaR;
                r1 = cmd->rminHist + (double)n*gd->deltaR;
            }
#if NDIM == 3
            shell = (4.0*PI/3.0)*(r1-r0)*(r1*r1+r1*r0+r0*r0);
#else
            shell = PI*(r1-r0)*(r1+r0);
#endif
            if (!isfinite(shell) || shell <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "count normalization: invalid shell volume in bin %d", n);
                return FAILURE;
            }
            /* histNN is still ordered here, hence no additional factor two. */
            const double xi = ((double)histNN[n]/count)/count*(volume/shell)-1.0;
            if (!isfinite(xi)) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "count normalization: non-finite correlation in bin %d", n);
                return FAILURE;
            }
            histCF[n] = xi;
        }
        histNN[n] *= 0.5;
    }
    return SUCCESS;
}

global int search_compute_HistN(struct cmdline_data *cmd,
                                struct global_data *gd, INTEGER nbody)
{
    return search_normalize_count_histograms(cmd, gd, nbody, gd->histNN, gd->histCF);
}
