/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2016 Krzysztof M. Gorski, Eric Hivon, Martin Reinecke,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann,
 *                          Reza Ansari & Kenneth M. Ganga
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.sourceforge.net
 *
 *---------------------------------------------------------------------------*/

//=============================================================================
//        1          2          3          4        ^ 5          6          7


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "chealpix.h"

static const double twothird=2.0/3.0;
static const double pi=3.141592653589793238462643383279502884197;
static const double twopi=6.283185307179586476925286766559005768394;
static const double halfpi=1.570796326794896619231321691639751442099;
static const double inv_halfpi=0.6366197723675813430755350534900574;

//B cBalls note:
// defines UTIL_ASSERT() through a helper that exits.
//  The active FITS reader wrappers look much better now,
//  but exported functions like ang2pix_ring() still abort on invalid theta.
//  This is lower priority unless those functions are used
//  from new addon paths, but it is still an unsafe library surface.
static void util_fail_ (const char *file, int line, const char *func,
  const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }

#if defined (__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond,msg) \
  if(!(cond)) util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
#define UTIL_FAIL(msg) \
  util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
//E

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
static double fmodulo (double v1, double v2)
  {
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
  }
/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
static int imodulo (int v1, int v2)
  { int v=v1%v2; return (v>=0) ? v : v+v2; }
static int isqrt(int v)
  { return (int)(sqrt(v+0.5)); }

/* ctab[m] = (short)(
       (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
    | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4)); */
static const short ctab[]={
#define Z(a) a,a+1,a+256,a+257
#define Y(a) Z(a),Z(a+2),Z(a+512),Z(a+514)
#define X(a) Y(a),Y(a+4),Y(a+1024),Y(a+1028)
X(0),X(8),X(2048),X(2056)
#undef X
#undef Y
#undef Z
};
/* utab[m] = (short)(
      (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
    | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7)); */
static const short utab[]={
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
X(0),X(1),X(4),X(5)
#undef X
#undef Y
#undef Z
};

static const int jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
static const int jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

#ifndef __BMI2__

static int xyf2nest (int nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside) +
      (utab[ix&0xff] | (utab[ix>>8]<<16)
    | (utab[iy&0xff]<<1) | (utab[iy>>8]<<17));
  }
static void nest2xyf (int nside, int pix, int *ix, int *iy, int *face_num)
  {
  int npface_=nside*nside, raw;
  *face_num = pix/npface_;
  pix &= (npface_-1);
  raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  *ix = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  pix >>= 1;
  raw = (pix&0x5555) | ((pix&0x55550000)>>15);
  *iy = ctab[raw&0xff] | (ctab[raw>>8]<<4);
  }
  
#else

#include <x86intrin.h>

static int xyf2nest (int nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside) +
    (_pdep_u32(ix,0x55555555u) | _pdep_u32(iy,0xaaaaaaaau));
  }
static void nest2xyf (int nside, int pix, int *ix, int *iy, int *face_num)
  {
  int npface_=nside*nside, raw;
  *face_num = pix/npface_;
  pix &= (npface_-1);
  *ix=_pext_u32(pix,0x55555555u);
  *iy=_pext_u32(pix,0xaaaaaaaau);
  }

#endif

static inline int special_div (int a, int b)
  {
  int t=(a>=(b<<1));
  a-=t*(b<<1);
  return (t<<1)+(a>=b);
  }

static int xyf2ring (int nside_, int ix, int iy, int face_num)
  {
  int nl4 = 4*nside_;
  int jr = (jrll[face_num]*nside_) - ix - iy  - 1, jp;

  int nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = 12*nside_*nside_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    int ncap_=2*nside_*(nside_-1);
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }
static void ring2xyf (int nside_, int pix, int *ix, int *iy, int *face_num)
  {
  int iring, iphi, kshift, nr;
  int ncap_=2*nside_*(nside_-1);
  int npix_=12*nside_*nside_;
  int nl2 = 2*nside_;

  if (pix<ncap_) /* North Polar cap */
    {
    iring = (1+isqrt(1+2*pix))>>1; /* counted from North pole */
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    *face_num = special_div(iphi-1,nr);
    }
  else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
    int ip = pix - ncap_;
    iring = (ip/(4*nside_)) + nside_; /* counted from North pole */
    iphi  = (ip%(4*nside_)) + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
    unsigned int ire = iring-nside_+1;
    unsigned int irm = nl2+2-ire;
    int ifm = (iphi - ire/2 + nside_ -1) / nside_;
    int ifp = (iphi - irm/2 + nside_ -1) / nside_;
    *face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else /* South Polar cap */
    {
    int ip = npix_ - pix;
    iring = (1+isqrt(2*ip-1))>>1; /* counted from South pole */
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    *face_num=8+special_div(iphi-1,nr);
    }

  int irt = iring - (jrll[*face_num]*nside_) + 1;
  int ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy =(-(ipt+irt))>>1;
  }

static int ang2pix_nest_z_phi (long nside_, double z, double phi)
  {
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */
  int face_num, ix, iy;

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*(z*0.75);
    int jp = (int)(temp1-temp2); /* index of  ascending edge line */
    int jm = (int)(temp1+temp2); /* index of descending edge line */
    int ifp = jp/nside_;  /* in {0,4} */
    int ifm = jm/nside_;
    face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));

    ix = jm & (nside_-1);
    iy = nside_ - (jp & (nside_-1)) - 1;
    }
  else /* polar region, za > 2/3 */
    {
    int ntt = (int)tt, jp, jm;
    double tp, tmp;
    if (ntt>=4) ntt=3;
    tp = tt-ntt;
    tmp = nside_*sqrt(3*(1-za));

    jp = (int)(tp*tmp); /* increasing edge line index */
    jm = (int)((1.0-tp)*tmp); /* decreasing edge line index */
    if (jp>=nside_) jp = nside_-1; /* for points too close to the boundary */
    if (jm>=nside_) jm = nside_-1;
    if (z >= 0)
      {
      face_num = ntt;  /* in {0,3} */
      ix = nside_ - jm - 1;
      iy = nside_ - jp - 1;
      }
    else
      {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
      }
    }

  return xyf2nest(nside_,ix,iy,face_num);
  }

static int ang2pix_ring_z_phi (long nside_, double z, double phi)
  {
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*z*0.75;
    int jp = (int)(temp1-temp2); /* index of  ascending edge line */
    int jm = (int)(temp1+temp2); /* index of descending edge line */

    /* ring number counted from z=2/3 */
    int ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
    int kshift = 1-(ir&1); /* kshift=1 if ir even, 0 otherwise */

    int ip = (jp+jm-nside_+kshift+1)/2; /* in {0,4n-1} */
    ip = imodulo(ip,4*nside_);

    return nside_*(nside_-1)*2 + (ir-1)*4*nside_ + ip;
    }
  else  /* North & South polar caps */
    {
    double tp = tt-(int)(tt);
    double tmp = nside_*sqrt(3*(1-za));

    int jp = (int)(tp*tmp); /* increasing edge line index */
    int jm = (int)((1.0-tp)*tmp); /* decreasing edge line index */

    int ir = jp+jm+1; /* ring number counted from the closest pole */
    int ip = (int)(tt*ir); /* in {0,4*ir-1} */
    ip = imodulo(ip,4*ir);

    if (z>0)
      return 2*ir*(ir-1) + ip;
    else
      return 12*nside_*nside_ - 2*ir*(ir+1) + ip;
    }
  }

static void pix2ang_ring_z_phi (int nside_, int pix, double *z, double *phi)
  {
  long ncap_=nside_*(nside_-1)*2;
  long npix_=12*nside_*nside_;
  double fact2_ = 4./npix_;
  if (pix<ncap_) /* North Polar cap */
    {
    int iring = (1+isqrt(1+2*pix))>>1; /* counted from North pole */
    int iphi  = (pix+1) - 2*iring*(iring-1);

    *z = 1.0 - (iring*iring)*fact2_;
    *phi = (iphi-0.5) * halfpi/iring;
    }
  else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
    double fact1_  = (nside_<<1)*fact2_;
    int ip  = pix - ncap_;
    int iring = ip/(4*nside_) + nside_; /* counted from North pole */
    int iphi  = ip%(4*nside_) + 1;
    /* 1 if iring+nside is odd, 1/2 otherwise */
    double fodd = ((iring+nside_)&1) ? 1 : 0.5;

    int nl2 = 2*nside_;
    *z = (nl2-iring)*fact1_;
    *phi = (iphi-fodd) * pi/nl2;
    }
  else /* South Polar cap */
    {
    int ip = npix_ - pix;
    int iring = (1+isqrt(2*ip-1))>>1; /* counted from South pole */
    int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

    *z = -1.0 + (iring*iring)*fact2_;
    *phi = (iphi-0.5) * halfpi/iring;
    }
  }

static void pix2ang_nest_z_phi (int nside_, int pix, double *z, double *phi)
  {
  int nl4 = nside_*4;
  int npix_=12*nside_*nside_;
  double fact2_ = 4./npix_;
  int face_num, ix, iy, nr, kshift;

  nest2xyf(nside_,pix,&ix,&iy,&face_num);
  int jr = (jrll[face_num]*nside_) - ix - iy - 1;

  if (jr<nside_)
    {
    nr = jr;
    *z = 1 - nr*nr*fact2_;
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    *z = nr*nr*fact2_ - 1;
    kshift = 0;
    }
  else
    {
    double fact1_ = (nside_<<1)*fact2_;
    nr = nside_;
    *z = (2*nside_-jr)*fact1_;
    kshift = (jr-nside_)&1;
    }

  int jp = (jpll[face_num]*nr + ix -iy + 1 + kshift) / 2;
  if (jp>nl4) jp-=nl4;
  if (jp<1) jp+=nl4;

  *phi = (jp-(kshift+1)*0.5)*(halfpi/nr);
  }

void ang2vec(double theta, double phi, double *vec)
  {
  double sz = sin(theta);
  vec[0] = sz * cos(phi);
  vec[1] = sz * sin(phi);
  vec[2] = cos(theta);
  }

void vec2ang(const double *vec, double *theta, double *phi)
  {
  *theta = atan2(sqrt(vec[0]*vec[0]+vec[1]*vec[1]),vec[2]);
  *phi = atan2 (vec[1],vec[0]);
  if (*phi<0.) *phi += twopi;
  }

long npix2nside(long npix)
  {
  long res = isqrt(npix/12);
  return (res*res*12==npix) ? res : -1;
  }

long nside2npix(const long nside)
  { return 12*nside*nside; }

void ang2pix_ring(long nside, double theta, double phi, long *ipix)
  {
  UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta out of range");
  *ipix=ang2pix_ring_z_phi (nside,cos(theta),phi);
  }
void ang2pix_nest(long nside, double theta, double phi, long *ipix)
  {
  UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta out of range");
  *ipix=ang2pix_nest_z_phi (nside,cos(theta),phi);
  }
void vec2pix_ring(long nside, const double *vec, long *ipix)
  {
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  *ipix=ang2pix_ring_z_phi (nside,vec[2]/vlen, atan2(vec[1],vec[0]));
  }
void vec2pix_nest(long nside, const double *vec, long *ipix)
  {
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  *ipix=ang2pix_nest_z_phi (nside,vec[2]/vlen, atan2(vec[1],vec[0]));
  }
void pix2ang_ring(long nside, long ipix, double *theta, double *phi)
  {
  double z;
  pix2ang_ring_z_phi (nside,ipix,&z,phi);
  *theta=acos(z);
  }
void pix2ang_nest(long nside, long ipix, double *theta, double *phi)
  {
  double z;
  pix2ang_nest_z_phi (nside,ipix,&z,phi);
  *theta=acos(z);
  }
void pix2vec_ring(long nside, long ipix, double *vec)
  {
  double z, phi;
  pix2ang_ring_z_phi (nside,ipix,&z,&phi);
  double stheta=sqrt((1.-z)*(1.+z));
  vec[0]=stheta*cos(phi);
  vec[1]=stheta*sin(phi);
  vec[2]=z;
  }
void pix2vec_nest(long nside, long ipix, double *vec)
  {
  double z, phi;
  pix2ang_nest_z_phi (nside,ipix,&z,&phi);
  double stheta=sqrt((1.-z)*(1.+z));
  vec[0]=stheta*cos(phi);
  vec[1]=stheta*sin(phi);
  vec[2]=z;
  }
void nest2ring(long nside, long ipnest, long *ipring)
  {
  int ix, iy, face_num;
  if ((nside&(nside-1))!=0) { *ipring=-1; return; }
  nest2xyf (nside, ipnest, &ix, &iy, &face_num);
  *ipring = xyf2ring (nside, ix, iy, face_num);
  }
void ring2nest(long nside, long ipring, long *ipnest)
  {
  int ix, iy, face_num;
  if ((nside&(nside-1))!=0) { *ipnest=-1; return; }
  ring2xyf (nside, ipring, &ix, &iy, &face_num);
  *ipnest = xyf2nest (nside, ix, iy, face_num);
  }

/* 64bit functions */

static int64_t imodulo64 (int64_t v1, int64_t v2)
  { int64_t v=v1%v2; return (v>=0) ? v : v+v2; }
static long isqrt64(int64_t v)
  {
  int64_t res = sqrt(v+0.5);
  if (v<((int64_t)(1)<<50)) return (long)res;
  if (res*res>v)
    --res;
  else if ((res+1)*(res+1)<=v)
    ++res;
  return (long)res;
  }

#ifndef __BMI2__
static int64_t spread_bits64 (int v)
  {
  return  (int64_t)(utab[ v     &0xff])
       | ((int64_t)(utab[(v>> 8)&0xff])<<16)
       | ((int64_t)(utab[(v>>16)&0xff])<<32)
       | ((int64_t)(utab[(v>>24)&0xff])<<48);
  }

static int64_t compress_bits64 (int64_t v)
  {
  int64_t raw = v&0x5555555555555555ull;
  raw|=raw>>15;
  return ctab[ raw     &0xff]      | (ctab[(raw>> 8)&0xff]<< 4)
      | (ctab[(raw>>32)&0xff]<<16) | (ctab[(raw>>40)&0xff]<<20);
  }

static int64_t xyf2nest64 (int64_t nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside) + spread_bits64(ix) + (spread_bits64(iy)<<1);
  }

static void nest2xyf64 (int64_t nside, int64_t pix, int *ix, int *iy,
  int *face_num)
  {
  int64_t npface_=nside*nside;
  *face_num = pix/npface_;
  pix &= (npface_-1);
  *ix = compress_bits64(pix);
  *iy = compress_bits64(pix>>1);
  }

#else

static int64_t xyf2nest64 (int64_t nside, int ix, int iy, int face_num)
  {
  return (face_num*nside*nside)
    + _pdep_u64(ix,0x5555555555555555ull) + _pdep_u64(iy,0xaaaaaaaaaaaaaaaaull);
  }

static void nest2xyf64 (int64_t nside, int64_t pix, int *ix, int *iy,
  int *face_num)
  {
  int64_t npface_=nside*nside;
  *face_num = pix/npface_;
  pix &= (npface_-1);
  *ix = _pext_u64(pix,0x5555555555555555ull);
  *iy = _pext_u64(pix,0xaaaaaaaaaaaaaaaaull);
  }

#endif

static inline int64_t special_div64 (int64_t a, int64_t b)
  {
  int64_t t=(a>=(b<<1));
  a-=t*(b<<1);
  return (t<<1)+(a>=b);
  }

static int64_t xyf2ring64 (int64_t nside_, int ix, int iy, int face_num)
  {
  int64_t nl4 = 4*nside_;
  int64_t jr = (jrll[face_num]*nside_) - ix - iy  - 1, jp;

  int64_t nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = 12*nside_*nside_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    int64_t ncap_=2*nside_*(nside_-1);
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  return n_before + jp - 1;
  }
static void ring2xyf64 (int64_t nside_, int64_t pix, int *ix, int *iy,
  int *face_num)
  {
  int64_t iring, iphi, kshift, nr;
  int64_t ncap_=2*nside_*(nside_-1);
  int64_t npix_=12*nside_*nside_;
  int64_t nl2 = 2*nside_;

  if (pix<ncap_) /* North Polar cap */
    {
    iring = (1+isqrt64(1+2*pix))>>1; /* counted from North pole */
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    *face_num=special_div64(iphi-1,nr);
    }
  else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
    int64_t ip = pix - ncap_;
    iring = (ip/(4*nside_)) + nside_; /* counted from North pole */
    iphi  = (ip%(4*nside_)) + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
    int64_t ire = iring-nside_+1;
    int64_t irm = nl2+2-ire;
    int64_t ifm = (iphi - ire/2 + nside_ -1) / nside_;
    int64_t ifp = (iphi - irm/2 + nside_ -1) / nside_;
    *face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else /* South Polar cap */
    {
    int64_t ip = npix_ - pix;
    iring = (1+isqrt64(2*ip-1))>>1; /* counted from South pole */
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    *face_num=8+special_div64(iphi-1,nr);
    }

  int64_t irt = iring - (jrll[*face_num]*nside_) + 1;
  int64_t ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy =(-(ipt+irt))>>1;
  }

static int64_t ang2pix_nest_z_phi64 (int64_t nside_, double z, double s,
  double phi)
  {
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */
  int face_num, ix, iy;

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*(z*0.75);
    int64_t jp = (int64_t)(temp1-temp2); /* index of  ascending edge line */
    int64_t jm = (int64_t)(temp1+temp2); /* index of descending edge line */
    int64_t ifp = jp/nside_;  /* in {0,4} */
    int64_t ifm = jm/nside_;
    face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));

    ix = jm & (nside_-1);
    iy = nside_ - (jp & (nside_-1)) - 1;
    }
  else /* polar region, za > 2/3 */
    {
    int ntt = (int)tt, jp, jm;
    double tp, tmp;
    if (ntt>=4) ntt=3;
    tp = tt-ntt;
    if (s>-2.)
      tmp = nside_*s/sqrt((1.+za)/3.);
    else
      tmp = nside_*sqrt(3*(1-za));

    jp = (int64_t)(tp*tmp); /* increasing edge line index */
    jm = (int64_t)((1.0-tp)*tmp); /* decreasing edge line index */
    if (jp>=nside_) jp = nside_-1; /* for points too close to the boundary */
    if (jm>=nside_) jm = nside_-1;
    if (z >= 0)
      {
      face_num = ntt;  /* in {0,3} */
      ix = nside_ - jm - 1;
      iy = nside_ - jp - 1;
      }
    else
      {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
      }
    }

  return xyf2nest64(nside_,ix,iy,face_num);
  }

static int64_t ang2pix_ring_z_phi64 (int64_t nside_, double z, double s,
  double phi)
  {
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*z*0.75;
    int64_t jp = (int64_t)(temp1-temp2); /* index of  ascending edge line */
    int64_t jm = (int64_t)(temp1+temp2); /* index of descending edge line */

    /* ring number counted from z=2/3 */
    int64_t ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
    int kshift = 1-(ir&1); /* kshift=1 if ir even, 0 otherwise */

    int64_t ip = (jp+jm-nside_+kshift+1)/2; /* in {0,4n-1} */
    ip = imodulo64(ip,4*nside_);

    return nside_*(nside_-1)*2 + (ir-1)*4*nside_ + ip;
    }
  else  /* North & South polar caps */
    {
    double tp = tt-(int)(tt);
    double tmp = (s>-2.) ? nside_*s/sqrt((1.+za)/3.) : nside_*sqrt(3*(1-za));

    int64_t jp = (int64_t)(tp*tmp); /* increasing edge line index */
    int64_t jm = (int64_t)((1.0-tp)*tmp); /* decreasing edge line index */

    int64_t ir = jp+jm+1; /* ring number counted from the closest pole */
    int64_t ip = (int64_t)(tt*ir); /* in {0,4*ir-1} */
    ip = imodulo64(ip,4*ir);

    if (z>0)
      return 2*ir*(ir-1) + ip;
    else
      return 12*nside_*nside_ - 2*ir*(ir+1) + ip;
    }
  }

static void pix2ang_ring_z_phi64 (int64_t nside_, int64_t pix,
  double *z, double *s, double *phi)
  {
  int64_t ncap_=nside_*(nside_-1)*2;
  int64_t npix_=12*nside_*nside_;
  double fact2_  = 4./npix_;
  *s=-5;
  if (pix<ncap_) /* North Polar cap */
    {
    int64_t iring = (1+isqrt64(1+2*pix))>>1; /* from N pole */
    int64_t iphi  = (pix+1) - 2*iring*(iring-1);
    double tmp=(iring*iring)*fact2_;

    *z = 1.0 - tmp;
    if (*z>0.99) *s=sqrt(tmp*(2.-tmp));
    *phi = (iphi-0.5) * halfpi/iring;
    }
  else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
    double fact1_  = (nside_<<1)*fact2_;
    int64_t ip  = pix - ncap_;
    int64_t iring = ip/(4*nside_) + nside_; /* counted from North pole */
    int64_t iphi  = ip%(4*nside_) + 1;
    /* 1 if iring+nside is odd, 1/2 otherwise */
    double fodd = ((iring+nside_)&1) ? 1 : 0.5;

    int64_t nl2 = 2*nside_;
    *z = (nl2-iring)*fact1_;
    *phi = (iphi-fodd) * pi/nl2;
    }
  else /* South Polar cap */
    {
    int64_t ip = npix_ - pix;
    int64_t iring = (1+isqrt64(2*ip-1))>>1; /* from S pole */
    int64_t iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

    double tmp=(iring*iring)*fact2_;
    *z = tmp - 1.0;
    if (*z<-0.99) *s=sqrt(tmp*(2.-tmp));
    *phi = (iphi-0.5) * halfpi/iring;
    }
  }

static void pix2ang_nest_z_phi64 (int64_t nside_, int64_t pix, double *z,
  double *s, double *phi)
  {
  int64_t nl4 = nside_*4;
  int64_t npix_=12*nside_*nside_;
  double fact2_ = 4./npix_;
  int face_num, ix, iy;
  int64_t jr, nr, kshift, jp;
  *s=-5;

  nest2xyf64(nside_,pix,&ix,&iy,&face_num);
  jr = (jrll[face_num]*nside_) - ix - iy - 1;

  if (jr<nside_)
    {
    double tmp;
    nr = jr;
    tmp=(nr*nr)*fact2_;
    *z = 1 - tmp;
    if (*z>0.99) *s=sqrt(tmp*(2.-tmp));
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    double tmp;
    nr = nl4-jr;
    tmp=(nr*nr)*fact2_;
    *z = tmp - 1;
    if (*z<-0.99) *s=sqrt(tmp*(2.-tmp));
    kshift = 0;
    }
  else
    {
    double fact1_ = (nside_<<1)*fact2_;
    nr = nside_;
    *z = (2*nside_-jr)*fact1_;
    kshift = (jr-nside_)&1;
    }

  jp = (jpll[face_num]*nr + ix -iy + 1 + kshift) / 2;
  if (jp>nl4) jp-=nl4;
  if (jp<1) jp+=nl4;

  *phi = (jp-(kshift+1)*0.5)*(halfpi/nr);
  }

long npix2nside64(int64_t npix)
  {
  int64_t res = isqrt64(npix/12.);
  return (res*res*12==npix) ? (long)res : -1;
  }

int64_t nside2npix64(int64_t nside)
  { return 12*nside*nside; }

void ang2pix_ring64(int64_t nside, double theta, double phi, int64_t *ipix)
  {
  UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta out of range");
  double cth=cos(theta), sth=(fabs(cth)>0.99) ? sin(theta) : -5;
  *ipix=ang2pix_ring_z_phi64 (nside,cth,sth,phi);
  }
void ang2pix_nest64(int64_t nside, double theta, double phi, int64_t *ipix)
  {
  UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta out of range");
  double cth=cos(theta), sth=(fabs(cth)>0.99) ? sin(theta) : -5;
  *ipix=ang2pix_nest_z_phi64 (nside,cth,sth,phi);
  }
void vec2pix_ring64(int64_t nside, const double *vec, int64_t *ipix)
  {
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  double cth = vec[2]/vlen;
  double sth=(fabs(cth)>0.99) ? sqrt(vec[0]*vec[0]+vec[1]*vec[1])/vlen : -5;
  *ipix=ang2pix_ring_z_phi64 (nside,cth,sth,atan2(vec[1],vec[0]));
  }
void vec2pix_nest64(int64_t nside, const double *vec, int64_t *ipix)
  {
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  double cth = vec[2]/vlen;
  double sth=(fabs(cth)>0.99) ? sqrt(vec[0]*vec[0]+vec[1]*vec[1])/vlen : -5;
  *ipix=ang2pix_nest_z_phi64 (nside,cth,sth,atan2(vec[1],vec[0]));
  }
void pix2ang_ring64(int64_t nside, int64_t ipix, double *theta, double *phi)
  {
  double z,s;
  pix2ang_ring_z_phi64 (nside,ipix,&z,&s,phi);
  *theta= (s<-2) ? acos(z) : atan2(s,z);
  }
void pix2ang_nest64(int64_t nside, int64_t ipix, double *theta, double *phi)
  {
  double z,s;
  pix2ang_nest_z_phi64 (nside,ipix,&z,&s,phi);
  *theta= (s<-2) ? acos(z) : atan2(s,z);
  }
void pix2vec_ring64(int64_t nside, int64_t ipix, double *vec)
  {
  double z, phi, stheta;
  pix2ang_ring_z_phi64 (nside,ipix,&z,&stheta,&phi);
  if (stheta<-2) stheta=sqrt((1.-z)*(1.+z));
  vec[0]=stheta*cos(phi);
  vec[1]=stheta*sin(phi);
  vec[2]=z;
  }
void pix2vec_nest64(int64_t nside, int64_t ipix, double *vec)
  {
  double z, phi, stheta;
  pix2ang_nest_z_phi64 (nside,ipix,&z,&stheta,&phi);
  if (stheta<-2) stheta=sqrt((1.-z)*(1.+z));
  vec[0]=stheta*cos(phi);
  vec[1]=stheta*sin(phi);
  vec[2]=z;
  }
void nest2ring64(int64_t nside, int64_t ipnest, int64_t *ipring)
  {
  int ix, iy, face_num;
  if ((nside&(nside-1))!=0) { *ipring=-1; return; }
  nest2xyf64 (nside, ipnest, &ix, &iy, &face_num);
  *ipring = xyf2ring64 (nside, ix, iy, face_num);
  }
void ring2nest64(int64_t nside, int64_t ipring, int64_t *ipnest)
  {
  int ix, iy, face_num;
  if ((nside&(nside-1))!=0) { *ipnest=-1; return; }
  ring2xyf64 (nside, ipring, &ix, &iy, &face_num);
  *ipnest = xyf2nest64 (nside, ix, iy, face_num);
  }

#include "fitsio.h"

static void setCoordSysHP(char coordsys,char *coordsys9)
  {
  strcpy(coordsys9,"C       ");
  if (coordsys=='G')
    strcpy (coordsys9,"G       ");
  else if (coordsys=='E')
    strcpy (coordsys9,"E       ");
  else if ((coordsys!='C')&&(coordsys!='Q'))
    fprintf(stderr, "%s (%d): System Cordinates are not correct"
                    "(Galactic,Ecliptic,Celestial=Equatorial). "
                    " Celestial system was set.\n", __FILE__, __LINE__);
  }


//B added by cBalls

#define CHEALPIX_STATUS_BAD_ARG   (-10001)
#define CHEALPIX_STATUS_BAD_HDU   (-10002)
#define CHEALPIX_STATUS_BAD_SHAPE (-10003)
#define CHEALPIX_STATUS_ALLOC     (-10004)
#define CHEALPIX_STATUS_BAD_ORDER (-10005)

static const char *chealpix_status_message(int status)
{
    switch (status) {
    case CHEALPIX_STATUS_BAD_ARG:
        return "bad argument";
    case CHEALPIX_STATUS_BAD_HDU:
        return "unexpected FITS HDU";
    case CHEALPIX_STATUS_BAD_SHAPE:
        return "invalid HEALPix FITS shape";
    case CHEALPIX_STATUS_ALLOC:
        return "allocation failure";
    case CHEALPIX_STATUS_BAD_ORDER:
        return "unsupported HEALPix ordering";
    default:
        return "unknown HEALPix/FITS error";
    }
}

static void chealpix_report_status(const char *where, int status)
{
    if (status == 0) return;

    if (status > 0) {
        fprintf(stderr, "%s: CFITSIO status=%d\n", where, status);
        fits_report_error(stderr, status);
    } else {
        fprintf(stderr, "%s: %s (%d)\n",
                where, chealpix_status_message(status), status);
    }
}

static int chealpix_is_space(char value)
{
    return value == ' ' || value == '\t' || value == '\r' || value == '\n';
}

static int chealpix_ordering_equals(const char *ordering,
                                    const char *expected)
{
    unsigned char actual;
    unsigned char wanted;

    if (ordering == NULL || expected == NULL)
        return 0;

    while (chealpix_is_space(*ordering))
        ordering++;

    while (*ordering != '\0' && *expected != '\0') {
        actual = (unsigned char)*ordering++;
        wanted = (unsigned char)*expected++;
        if (actual >= 'a' && actual <= 'z')
            actual = (unsigned char)(actual - 'a' + 'A');
        if (wanted >= 'a' && wanted <= 'z')
            wanted = (unsigned char)(wanted - 'a' + 'A');
        if (actual != wanted)
            return 0;
    }

    if (*expected != '\0')
        return 0;
    while (chealpix_is_space(*ordering))
        ordering++;

    return *ordering == '\0';
}

int healpix_map_to_ring_status(float **map_io, long nside,
                               const char *ordering)
{
    float *ring_map;
    long ipnest;
    long ipring;
    long npix;

    if (map_io == NULL || *map_io == NULL || ordering == NULL || nside <= 0)
        return CHEALPIX_STATUS_BAD_ARG;
    if (nside > LONG_MAX / nside / 12)
        return CHEALPIX_STATUS_BAD_SHAPE;

    npix = nside2npix(nside);
    if (chealpix_ordering_equals(ordering, "RING"))
        return 0;
    if (!chealpix_ordering_equals(ordering, "NESTED"))
        return CHEALPIX_STATUS_BAD_ORDER;
    if ((nside & (nside - 1)) != 0)
        return CHEALPIX_STATUS_BAD_SHAPE;
    if ((size_t)npix > SIZE_MAX / sizeof(*ring_map))
        return CHEALPIX_STATUS_BAD_SHAPE;

    ring_map = (float *)malloc((size_t)npix * sizeof(*ring_map));
    if (ring_map == NULL)
        return CHEALPIX_STATUS_ALLOC;

    for (ipnest = 0; ipnest < npix; ipnest++) {
        nest2ring(nside, ipnest, &ipring);
        if (ipring < 0 || ipring >= npix) {
            free(ring_map);
            return CHEALPIX_STATUS_BAD_SHAPE;
        }
        ring_map[ipring] = (*map_io)[ipnest];
    }

    free(*map_io);
    *map_io = ring_map;
    return 0;
}

int read_healpix_map_status(const char *infile, long *nside,
                             char *coordsys, char *ordering,
                             float **map_out)
{
    long naxes = 0, *naxis = NULL, npix = 0;
    int status = 0, close_status = 0, hdutype = 0, nfound = 0, anynul = 0;
    float nulval = HEALPIX_NULLVAL, *map = NULL;
    fitsfile *fptr = NULL;

    if (!infile || !nside || !coordsys || !ordering || !map_out)
        return CHEALPIX_STATUS_BAD_ARG;
    *map_out = NULL;

    fits_open_file(&fptr, infile, READONLY, &status);
    if (status) goto fail;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    if (status) goto fail;
    if (hdutype != BINARY_TBL) { status = CHEALPIX_STATUS_BAD_HDU; goto fail; }

    fits_read_key_lng(fptr, "NAXIS", &naxes, NULL, &status);
    if (status) goto fail;
    if (naxes < 2) { status = CHEALPIX_STATUS_BAD_SHAPE; goto fail; }

    naxis = (long *)malloc((size_t)naxes * sizeof(long));
    if (!naxis) { status = CHEALPIX_STATUS_ALLOC; goto fail; }

    fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status);
    if (status) goto fail;
    if (nfound != naxes) { status = CHEALPIX_STATUS_BAD_SHAPE; goto fail; }

    fits_read_key_lng(fptr, "NSIDE", nside, NULL, &status);
    if (status) goto fail;
    if (*nside <= 0) { status = CHEALPIX_STATUS_BAD_SHAPE; goto fail; }

    npix = 12 * (*nside) * (*nside);
    if (naxis[1] == 0 || (npix % naxis[1]) != 0) {
        status = CHEALPIX_STATUS_BAD_SHAPE;
        goto fail;
    }

    status = 0;
    if (fits_read_key(fptr, TSTRING, "COORDSYS", coordsys, NULL, &status)) {
        strcpy(coordsys, "C");
        status = 0;
    }
    if (fits_read_key(fptr, TSTRING, "ORDERING", ordering, NULL, &status)) {
        strcpy(ordering, "RING");
        status = 0;
    }

    map = (float *)malloc((size_t)npix * sizeof(float));
    if (!map) { status = CHEALPIX_STATUS_ALLOC; goto fail; }

    fits_read_col(fptr, TFLOAT, 1, 1, 1, npix, &nulval, map, &anynul, &status);
    if (status) goto fail;

    fits_close_file(fptr, &close_status);
    fptr = NULL;
    if (close_status) { status = close_status; goto fail; }

    free(naxis);
    *map_out = map;
    return 0;

fail:
    if (fptr) {
        close_status = 0;
        fits_close_file(fptr, &close_status);
    }
    free(naxis);
    free(map);
    return status ? status : CHEALPIX_STATUS_BAD_SHAPE;
}

int get_fits_size_status(const char *filename, long *nside,
                         char *ordering, long *npix_out)
{
    fitsfile *fptr = NULL;
    int status = 0, close_status = 0, hdutype = 0;
    long obs_npix = 0;

    if (!filename || !nside || !ordering || !npix_out)
        return CHEALPIX_STATUS_BAD_ARG;
    *npix_out = 0;

    fits_open_file(&fptr, filename, READONLY, &status);
    if (status) goto done;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    if (status) goto done;

    fits_read_key(fptr, TSTRING, "ORDERING", ordering, NULL, &status);
    if (status) goto done;
    fits_read_key(fptr, TLONG, "NSIDE", nside, NULL, &status);
    if (status) goto done;

    if (fits_read_key(fptr, TLONG, "OBS_NPIX", &obs_npix, NULL, &status)) {
        obs_npix = 12 * (*nside) * (*nside);
        status = 0;
    }

done:
    if (fptr) {
        fits_close_file(fptr, &close_status);
        if (!status && close_status) status = close_status;
    }
    if (!status) *npix_out = obs_npix;
    return status;
}

int write_healpix_map_status(const float *signal, long nside,
                             const char *filename, char nest,
                             const char *coordsys)
{
    fitsfile *fptr = NULL;
    int status = 0, close_status = 0, hdutype = 0;
    long naxes[] = {0, 0};
    long npix = 12L * nside * nside;
    char order[9], coordsys9[9];
    char *ttype[] = {"SIGNAL"}, *tform[] = {"1E"}, *tunit[] = {" "};

    if (!signal || !filename || !coordsys || strlen(coordsys) < 1 || nside <= 0)
        return CHEALPIX_STATUS_BAD_ARG;

    fits_create_file(&fptr, filename, &status);
    if (status) goto done;
    fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);
    fits_write_date(fptr, &status);
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    if (status) goto done;

    fits_create_tbl(fptr, BINARY_TBL, npix, 1, ttype, tform, tunit,
                    "BINTABLE", &status);
    strcpy(order, nest ? "NESTED  " : "RING    ");
    setCoordSysHP(coordsys[0], coordsys9);

    fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX",
                   "HEALPIX Pixelisation", &status);
    fits_write_key(fptr, TSTRING, "ORDERING", order,
                   "Pixel ordering scheme, either RING or NESTED", &status);
    fits_write_key(fptr, TLONG, "NSIDE", &nside,
                   "Resolution parameter for HEALPIX", &status);
    fits_write_key(fptr, TSTRING, "COORDSYS", coordsys9,
                   "Pixelisation coordinate system", &status);
    fits_write_col(fptr, TFLOAT, 1, 1, 1, npix, (void *)signal, &status);

done:
    if (fptr) {
        fits_close_file(fptr, &close_status);
        if (!status && close_status) status = close_status;
    }
    return status;
}

float *read_healpix_map(const char *infile, long *nside, char *coordsys,
                        char *ordering)
{
    float *map = NULL;
    int status = read_healpix_map_status(infile, nside, coordsys,
                                         ordering, &map);

    if (status != 0) {
        chealpix_report_status("read_healpix_map", status);
        return NULL;
    }

    return map;
}

long get_fits_size(const char *filename, long *nside, char *ordering)
{
    long npix = -1;
    int status = get_fits_size_status(filename, nside, ordering, &npix);

    if (status != 0) {
        chealpix_report_status("get_fits_size", status);
        return -1;
    }

    return npix;
}

void write_healpix_map(const float *signal, long nside, const char *filename,
                       char nest, const char *coordsys)
{
    int status = write_healpix_map_status(signal, nside, filename,
                                          nest, coordsys);

    if (status != 0)
        chealpix_report_status("write_healpix_map", status);
}

//B activate later if necessary
//#undef CHEALPIX_STATUS_BAD_ARG
//#undef CHEALPIX_STATUS_BAD_HDU
//#undef CHEALPIX_STATUS_BAD_SHAPE
//#undef CHEALPIX_STATUS_ALLOC
//E

//E
