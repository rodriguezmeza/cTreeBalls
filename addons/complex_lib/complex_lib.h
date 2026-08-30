#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {double r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

//#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(double re, double im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b);
double Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(double x, fcomplex a);

#define II Complex(0.0,1.0)

#define CLRMC_ext(p,dim)                                                \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 1; _i <= dim; _i++)                                       \
        for (_j = 1; _j <= dim; _j++)                                   \
            (p)[_i][_j] = Complex(0.0,0.0);                             \
}


#define CLRMC_ext_ext(p,dim,dim2)                                       \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 1; _i <= dim; _i++)                                       \
        for (_j = 1; _j <= dim2; _j++)                                  \
            (p)[_i][_j] = Complex(0.0,0.0);                             \
}

#define OUTVPC_ext(p,v,u,dim)                                            \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 1; _i <= dim; _i++)                                       \
        for (_j = 1; _j <= dim; _j++)                                   \
            (p)[_i][_j] = Cmul((v)[_i], (u)[_j]);                            \
}

#define MULMCRS_ext(p,q,s,dim)                                            \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 1; _i <= dim; _i++)                                       \
        for (_j = 1; _j <= dim; _j++)                                   \
            (p)[_i][_j] = RCmul((s), (q)[_i][_j]);                     \
}

#define ADDMC_ext(p,q,r,dim)                                             \
{                                                                       \
    int _i, _j;                                                         \
    for (_i = 1; _i <= dim; _i++)                                       \
        for (_j = 1; _j <= dim; _j++)                                   \
            (p)[_i][_j] = Cadd((q)[_i][_j], (r)[_i][_j]);                    \
}

fcomplex *cmplx_vector(long nl, long nh);
fcomplex **cmplx_matrix(long nrl, long nrh, long ncl, long nch);
fcomplex ***cmplx_matrix3D(long nrl, long nrh, long ncl, long nch, long nc2l, long nc2h);

void free_cmplx_vector(fcomplex *v, long nl, long nh);
void free_cmplx_matrix(fcomplex **m, long nrl, long nrh, long ncl, long nch);
void free_cmplx_matrix3D(fcomplex ***m, long nrl, long nrh, long ncl, long nch, long nc2l, long nc2h);

#endif /* _NR_COMPLEX_H_ */
