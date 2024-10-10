
#include "stdinc.h"
#include "mathfns.h"
#include "constant.h"
#include "mathutil.h"
#include "vectmath.h"

void RotationVecAWRtoVecB(real C[], real A[], real B[], real theta)
{
    real x1, x2;
    real magA, magW;

#if defined(THREEDIM)
    vector W;
    CROSSVP(W,B,A);
    ABSV(magA, A);
    ABSV(magW, W);
    x1 = rcos(theta)/magA;
    x2 = rsin(theta)/magW;
    C[0] = magA * ( x1*A[0] + x2*W[0] );
    C[1] = magA * ( x1*A[1] + x2*W[1] );
    C[2] = magA * ( x1*A[2] + x2*W[2] );
//    printf("%g %g %g %g\n",theta,C[0],C[1],C[2]);
#else
// Check that these lines applied to 2D case!!!!
    real W;
    CROSSVP(W,B,A);
    ABSV(magA, A);
    magW = ABS(W);
    x1 = rcos(theta)/magA;
    x2 = rsin(theta)/magW;
    C[0] = magA * ( x1*A[0] + x2*W );
    C[1] = magA * ( x1*A[1] + x2*W );
//    printf("%g %g %g %g\n",theta,C[0],C[1]);
#endif

}

// https://en.wikipedia.org/wiki/Euler_angles
// Transpuesta de Z1Y2Z3
// Trans(R) = Inv(R)
void Rotation3D(real vec[], real alpha, real beta, real gamma)
{
    real a11, a12, a13;
    real a21, a22, a23;
    real a31, a32, a33;
    real vecp[3];
    real alp, bet, gam;
    int i, ndim;

    ndim=3;

    alp = radian(alpha);
    bet = radian(beta);
    gam = radian(gamma);

    a11 = rcos(gam)*rcos(bet)*rcos(alp) - rsin(gam)*rsin(alp);
    a12 = rcos(gam)*rcos(bet)*rsin(alp) + rsin(gam)*rcos(alp);
    a13 = -rcos(gam)*rsin(bet);

    a21 = -rsin(gam)*rcos(bet)*rcos(alp) - rcos(gam)*rsin(alp);
    a22 = -rsin(gam)*rcos(bet)*rsin(alp) + rcos(gam)*rcos(alp);
    a23 = rsin(gam)*rsin(bet);

    a31 = rsin(bet)*rcos(alp);
    a32 = rsin(bet)*rsin(alp);
    a33 = rcos(bet);

    vecp[0] = a11*vec[0] + a12*vec[1] + a13*vec[2];
    vecp[1] = a21*vec[0] + a22*vec[1] + a23*vec[2];
    vecp[2] = a31*vec[0] + a32*vec[1] + a33*vec[2];

    for (i = 0; i < ndim; i++)
        vec[i] = vecp[i];
}

/*
// Small rotations: Angles in radians
// Neglect: O(eps^3)
void Rotation3D(real vec[], real alpha, real beta, real gamma)
{
    real a11, a12, a13;
    real a21, a22, a23;
    real a31, a32, a33;
    real vecp[3];
    real alp, bet, gam;
    int i, ndim;

    ndim=3;

    alp = radian(alpha);
    bet = radian(beta);
    gam = radian(gamma);

    a11 = rcos(gam)*rcos(bet)*rcos(alp) - rsin(gam)*rsin(alp);
    a12 = rcos(gam)*rcos(bet)*rsin(alp) + rsin(gam)*rcos(alp);
    a13 = -rcos(gam)*rsin(bet);

    a21 = -rsin(gam)*rcos(bet)*rcos(alp) - rcos(gam)*rsin(alp);
    a22 = -rsin(gam)*rcos(bet)*rsin(alp) + rcos(gam)*rcos(alp);
    a23 = rsin(gam)*rsin(bet);

    a31 = rsin(bet)*rcos(alp);
    a32 = rsin(bet)*rsin(alp);
    a33 = rcos(bet);

    vecp[0] = a11*vec[0] + a12*vec[1] + a13*vec[2];
    vecp[1] = a21*vec[0] + a22*vec[1] + a23*vec[2];
    vecp[2] = a31*vec[0] + a32*vec[1] + a33*vec[2];

    for (i = 0; i < ndim; i++)
        vec[i] = vecp[i];
}
*/
// Small rotations: Angles in radians
#ifdef SINGLEP
void dRotation3D(float vec[], real alpha, real beta, real gamma, float *vecp)
#else
void dRotation3D(real vec[], real alpha, real beta, real gamma, real *vecp)
#endif
{
    real a11, a12, a13;
    real a21, a22, a23;
    real a31, a32, a33;
//    real vecp[3];
//    int i, ndim;

//    ndim=3;

    a11 = 1.0 - gamma*alpha;
    a12 = alpha + gamma;
    a13 = -beta;

//    a21 = -gamma - alpha;
    a21 = -a12;
//    a22 = -gamma*alpha + 1.0;
    a22 = a11;
    a23 = gamma*beta;

    a31 = beta;
    a32 = beta*alpha;
    a33 = 1.0;

    vecp[0] = a11*vec[0] + a12*vec[1] + a13*vec[2];
    vecp[1] = a21*vec[0] + a22*vec[1] + a23*vec[2];
    vecp[2] = a31*vec[0] + a32*vec[1] + a33*vec[2];

//    for (i = 0; i < ndim; i++)
//        vec[i] = vecp[i];
}

// Small rotations (linear order): Angle in radians
void d1Rotation3D(real vec[], real alpha)
{
//    real a11;
    real a12, a13;
    real a21;
//    , a22, a23;
    real a31;
//    , a32, a33;
    real vecp[3];
    int i, ndim;

    ndim=3;

//    a11 = 1.0 - alpha*alpha;
//    a11 = 1.0;
    a12 = 2.0*alpha;
    a13 = -alpha;

    a21 = -a12;
//    a22 = a11;
//    a23 = alpha*alpha;
//    a23 = 0.0;

    a31 = alpha;
//    a32 = a23;
//    a33 = 1.0;

//    vecp[0] = a11*vec[0] + a12*vec[1] + a13*vec[2];
    vecp[0] = vec[0] + a12*vec[1] + a13*vec[2];
//    vecp[1] = a21*vec[0] + a22*vec[1] + a23*vec[2];
//    vecp[1] = a21*vec[0] + a22*vec[1];
    vecp[1] = a21*vec[0] + vec[1];
//    vecp[2] = a31*vec[0] + a32*vec[1] + a33*vec[2];
//    vecp[2] = a31*vec[0] + a33*vec[2];
    vecp[2] = a31*vec[0] + vec[2];

    for (i = 0; i < ndim; i++)
        vec[i] = vecp[i];
}

// Cross vector product: vec2-vec1 vs vec3-vec1
bool crossVecProdSign(real vec1[], real vec2[], real vec3[])
{
    vector vec21, vec31, vec;
    real s;
    SUBV(vec21,vec2,vec1);
    SUBV(vec31,vec3,vec1);
    CROSSVP(vec,vec21,vec31);
    DOTVP(s,vec,vec1);
    if (s<0.0)
        return TRUE;
    else
        return FALSE;
}

real radian(real degree)
{
	return(degree * PI_D / ((real) 180.0) );
}

real degree(real radian)
{
	return(radian * ((real) 180.0) / PI_D );
}

real angle(real XI, real YI, real XF, real YF)
{
	real DX,DY, ang;

    DX = XF - XI;
    DY = YF - YI;
    if (DY == ZERO) {
		ang = ZERO;
        if ( DX < ZERO ) ang = ang + PI_D;
        return(ang);
    }
    if (DX == ZERO) {
		ang = PI_D/((real) 2.0);
        if ( DY < ZERO ) ang = ang + PI_D;
        return(ang);
	}
    ang = ratan( rabs(DY/DX) );
    if (DX < ZERO && DY > ZERO) {
        ang=PI_D-ang;
        return(ang);
	}
    if (DX < ZERO && DY < ZERO) {
		ang=PI_D+ang;
        return(ang);
	}
	if (DX > ZERO && DY < ZERO) {
        ang= ((real) 2.0) * PI_D - ang;
        return(ang);
	}
	return(ang);
}

real angle_dxdy(real DX, real DY)
{
    real ang;

//    DX = XF - XI;
//    DY = YF - YI;
    if (DY == ZERO) {
        ang = ZERO;
        if ( DX < ZERO ) ang = ang + PI_D;
        return(ang);
    }
    if (DX == ZERO) {
        ang = PI_D/((real) 2.0);
        if ( DY < ZERO ) ang = ang + PI_D;
        return(ang);
    }
    ang = ratan( rabs(DY/DX) );
    if (DX < ZERO && DY > ZERO) {
        ang=PI_D-ang;
        return(ang);
    }
    if (DX < ZERO && DY < ZERO) {
        ang=PI_D+ang;
        return(ang);
    }
    if (DX > ZERO && DY < ZERO) {
        ang= ((real) 2.0) * PI_D - ang;
        return(ang);
    }
    return(ang);
}


/*
double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}


double dmin(double x,double y)
{
  if(x<y)
    return x;
  else
    return y;
}


int imax(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}


int imin(int x,int y)
{
  if(x<y)
    return x;
  else
    return y;
}
*/

/*
void int_piksrt(int n, int arr[])	// Ordena un arreglo de enteros
{									// de menor a mayor.
	int i, j;
	int a;

	for (j=1; j<n; j++) {
		a = arr[j];
		i=j-1;
		while (i>=0 && arr[i] > a) {
			arr[i+1]=arr[i];
			i--;
		}
		arr[i+1]=a;
	}
}
*/

// Sorting using the Heapsort Method (pr_anvorpol.c) (Ascending order)

void HeapSort (real *a, int *seq, int n)
{
  real q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] < a[seq[j + 1]]) ++ j;
        if (q < a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting an array of integers using the Heapsort Method (pr_anvorpol.c) (Ascending order)

void HeapSortInt (int *a, int *seq, int n)
{
  int q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] < a[seq[j + 1]]) ++ j;
        if (q < a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting using the Heapsort Method (Descending order)

void HeapSortDescend (real *a, int *seq, int n)
{
  real q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] > a[seq[j + 1]]) ++ j;
        if (q > a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Sorting integer array using the Heapsort Method (Descending order)

void HeapSortIntDescend (int *a, int *seq, int n)
{
  int q;
  int i, ir, ixt, j, k;

  for (j = 0; j < n; j ++) seq[j] = j;
  if (n > 1) {
    k = n / 2;
    ir = n - 1;
    while (1) {
      if (k > 0) {
        -- k;
        ixt = seq[k];
        q = a[ixt];
      } else {
        ixt = seq[ir];
        q = a[ixt];
        seq[ir] = seq[0];
        -- ir;
        if (ir == 0) {
          seq[0] = ixt;
          break;
        }
      }
      i = k;
      j = 2 * k + 1;
      while (j <= ir) {
        if (j < ir && a[seq[j]] > a[seq[j + 1]]) ++ j;
        if (q > a[seq[j]]) {
          seq[i] = seq[j];
          i = j;
          j = 2 * j + 1;
        } else j = ir + 1;
      }
      seq[i] = ixt;
    }
  }
}

// Fast Fourier Transform (pr_anspcor.c)
// Use: FftComplex (a, n);
//  Cmplx *work;
//  AllocMem (work, 2 * (nValCorr - 1), Cmplx);
//  CSet (work[n], corrSum[j][k * nValCorr + n] * damp, 0.);
//  work[n] = work[2 * (nValCorr - 1) - n];
void FftComplex (Cmplx *a, int size)
{
  Cmplx t, w, wo;
  real theta;
  int i, j, k, n;

  k = 0;
  for (i = 0; i < size; i ++) {
    if (i < k) {
      t = a[i];
      a[i] = a[k];
      a[k] = t;
    }
    n = size / 2;
    while (n >= 1 && k >= n) {
      k -= n;
      n /= 2;
    }
    k += n;
  }
  for (n = 1; n < size; n *= 2) {
    theta = M_PI / n;
    CSet (wo, cos (theta) - 1., sin (theta));
    CSet (w, 1., 0.);
    for (k = 0; k < n; k ++) {
      for (i = k; i < size; i += 2 * n) {
        j = i + n;
        CMul (t, w, a[j]);
        CSub (a[j], a[i], t);
        CAdd (a[i], a[i], t);
      }
      CMul (t, w, wo);
      CAdd (w, w, t);
    }
  }
}

real Integrate (real *f, int nf)
{
  real s;
  int i;

  s = 0.5 * (f[0] + f[nf - 1]);
  for (i = 1; i < nf - 1; i ++) s += f[i];
  return (s);
}


