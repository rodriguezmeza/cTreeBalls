/* Shared numerical kernel for native and saved-file scalar windows. */
#ifndef CBALLS_SCALAR_WINDOW_SOLVER_H
#define CBALLS_SCALAR_WINDOW_SOLVER_H
#include <complex.h>
#include <float.h>
#include <math.h>
#include <string.h>

/* C99 complex representation is two adjacent doubles. Construct without
 * arithmetic so infinite imaginary components cannot contaminate the real part.
 * Some supported C99 system headers do not provide the C11 CMPLX macro. */
static inline double complex cballs_scalar_complex(double re, double im)
{
    double complex value;
    const double parts[2] = {re, im};
    memcpy(&value, parts, sizeof(value));
    return value;
}

/* Scaled complex LU with partial pivoting. Preserve the acceptance threshold;
 * expose rejected solves rather than publishing them as measured zeros. */
static int cballs_scalar_edge_solve(double complex *a, double complex *rhs, int n,
                                double *pivot_ratio)
{
    double smallest_pivot = DBL_MAX, largest_pivot = 0.0;
    *pivot_ratio = NAN;
    const double tolerance = 128.0 * DBL_EPSILON * n;
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double largest = cabs(a[(size_t)col * n + col]);
        for (int row = col + 1; row < n; row++) {
            const double value = cabs(a[(size_t)row * n + col]);
            if (!isfinite(value)) return CBALLS_WINDOW_NONFINITE;
            if (value > largest) {
                pivot = row;
                largest = value;
            }
        }
        if (!isfinite(largest)) return CBALLS_WINDOW_NONFINITE;
        if (largest <= tolerance) {
            *pivot_ratio = 0.0;
            return CBALLS_WINDOW_SINGULAR;
        }
        smallest_pivot = fmin(smallest_pivot, largest);
        largest_pivot = fmax(largest_pivot, largest);
        if (pivot != col) {
            for (int j = col; j < n; j++) {
                const double complex swap = a[(size_t)col * n + j];
                a[(size_t)col * n + j] = a[(size_t)pivot * n + j];
                a[(size_t)pivot * n + j] = swap;
            }
            const double complex swap = rhs[col];
            rhs[col] = rhs[pivot];
            rhs[pivot] = swap;
        }
        for (int row = col + 1; row < n; row++) {
            const double complex factor = a[(size_t)row * n + col]
                                        / a[(size_t)col * n + col];
            for (int j = col + 1; j < n; j++)
                a[(size_t)row * n + j] -= factor * a[(size_t)col * n + j];
            rhs[row] -= factor * rhs[col];
        }
    }
    for (int row = n - 1; row >= 0; row--) {
        for (int j = row + 1; j < n; j++)
            rhs[row] -= a[(size_t)row * n + j] * rhs[j];
        rhs[row] /= a[(size_t)row * n + row];
        if (!isfinite(creal(rhs[row])) || !isfinite(cimag(rhs[row])))
            return CBALLS_WINDOW_NONFINITE;
    }
    *pivot_ratio = smallest_pivot / largest_pivot;
    return CBALLS_WINDOW_VALID;
}

/* Nonnegative modes are contiguous. Negative modes follow by conjugation.
 * Caller owns n*n matrix and n-element rhs scratch, n=2*mmax+1.
 * A rejected bin never contains a published solution. */
static inline int cballs_scalar_edge_deconvolve(
    const double complex *signal, const double complex *window, int mmax,
    double complex *matrix, double complex *rhs, double *pivot_ratio)
{
    const int n = 2*mmax+1;
    *pivot_ratio = NAN;
    for (int m=0; m<n; m++)
        if (!isfinite(creal(window[m])) || !isfinite(cimag(window[m])))
            return CBALLS_WINDOW_NONFINITE;
    for (int m=0; m<=mmax; m++)
        if (!isfinite(creal(signal[m])) || !isfinite(cimag(signal[m])))
            return CBALLS_WINDOW_NONFINITE;
    const double w0 = creal(window[0]);
    if (!(w0 > 0)) return CBALLS_WINDOW_EMPTY;
    for (int row=0; row<n; row++) {
        const int ell = row-mmax;
        rhs[row] = (ell < 0 ? conj(signal[-ell]) : signal[ell])/w0;
        for (int col=0; col<n; col++) {
            const int difference = row-col;
            matrix[(size_t)row*n+col] = (difference < 0
                ? conj(window[-difference]) : window[difference])/w0;
        }
    }
    return cballs_scalar_edge_solve(matrix, rhs, n, pivot_ratio);
}
#endif
