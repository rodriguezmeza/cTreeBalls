
#include "stdinc.h"
#include "mathfns.h"
#include "constant.h"
#include "mathutil.h"
#include "vectmath.h"
#include "numrec.h"

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

// consider using atan2: tan-1(y/x) = atan2(y,x)
//  atan2() takes two arguments: x-coordinate and y-coordinate,
//  and calculate the angle in radians for the quadrant.
//  is defined in <math.h> header file.
//  Caution while using atan2():
//      The value of second argument passed should not be 0.
//      If the second argument passed is 0, the program will not run correctly.
// example:
//  y = 2.53; x = -10.2; result = atan2(y, x); result = result * 180.0/PI;
// will give:
//      Tangent inverse for(x = -10.2, y = 2.53) is 166.1 degrees.
// here is a routine to compute theta and phi from a vector:
/*
 void vec2ang(const double *vec, double *theta, double *phi)
   {
   *theta = atan2(sqrt(vec[0]*vec[0]+vec[1]*vec[1]),vec[2]);
   *phi = atan2 (vec[1],vec[0]);
   if (*phi<0.) *phi += twopi;
   }
 */
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

bool cross_product(real DX, real DY)
{
    bool answer;


    return(answer);
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


void covarianceMatrix(double *x, double *y, int n, double **cvm)
{
    int j;
    real xm;
    real ym;
    real s;
    real covxx;
    real covxy;
    real covyx;
    real covyy;

    if (n <= 1) error("n must be at least 2...");

    s=0.0;
    for (j=1;j<=n;j++) s += x[j];
    xm=s/((real)n);
    s=0.0;
    for (j=1;j<=n;j++) s += y[j];
    ym=s/((real)n);

    covxx = 0;
    for (j=1;j<=n;j++) covxx += (x[j]-xm)*(x[j]-xm);
    covxx /= (real)(n-1);

    covxy = 0;
    for (j=1;j<=n;j++) covxy += (x[j]-xm)*(y[j]-ym);
    covxy /= (real)(n-1);

    covyx = 0;
    for (j=1;j<=n;j++) covyx += (y[j]-ym)*(x[j]-xm);
    covyx /= (real)(n-1);

    covyy = 0;
    for (j=1;j<=n;j++) covyy += (y[j]-ym)*(y[j]-ym);
    covyy /= (real)(n-1);

    cvm[1][1] = covxx;
    cvm[1][2] = covxy;
    cvm[2][1] = covyx;
    cvm[2][2] = covyy;
}

//B Routines below are from NR
//  Brought here to avoid issues when GSLON = 1

//B use as:
//  moment(data,i,&ave,&adev,&sdev,&vrnce,&skew,&curt);
//  Given an array of data[1..n], this routine returns
//  its mean ave, average deviation adev,
//  standard deviation sdev, variance var,
//  skewness skew, and kurtosis curt.
//E
void moment(double data[], int n, double *ave, double *adev, double *sdev,
    double *var, double *skew, double *curt)
{
//    void nrerror(char error_text[]);
    int j;
    double ep=0.0,s,p;

    if (n <= 1) error("n must be at least 2 in moment");
    s=0.0;
    for (j=1;j<=n;j++) s += data[j];
    *ave=s/n;
    *adev=(*var)=(*skew)=(*curt)=0.0;
    for (j=1;j<=n;j++) {
        *adev += fabs(s=data[j]-(*ave));
        ep += s;
        *var += (p=s*s);
        *skew += (p *= s);
        *curt += (p *= s);
    }
    *adev /= n;
    *var=(*var-ep*ep/n)/(n-1);
    *sdev=sqrt(*var);
    if (*var) {
        *skew /= (n*(*var)*(*sdev));
        *curt=(*curt)/(n*(*var)*(*var))-3.0;
    } else error("No skew/kurtosis when variance = 0 (in moment)");
}


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void gaussj(double **a, int n, double **b, int m)
{
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv,temp;

    indxc=ivector(1,n);
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if (ipiv[j] != 1)
                for (k=1;k<=n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
            for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a[icol][icol] == 0.0) error("gaussj: Singular Matrix");
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=1;l<=n;l++) a[icol][l] *= pivinv;
        for (l=1;l<=m;l++) b[icol][l] *= pivinv;
        for (ll=1;ll<=n;ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
            }
    }
    for (l=n;l>=1;l--) {
        if (indxr[l] != indxc[l])
            for (k=1;k<=n;k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
}
#undef SWAP

void lubksb(double **a, int n, int *indx, double b[])
{
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}

#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;

    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) error("Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}
#undef TINY

//E Routines below are from NR
