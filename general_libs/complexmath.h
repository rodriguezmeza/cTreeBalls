
#ifndef _complexmath_h
#define _complexmath_h

#define CMul(cr, ci, ar, ai, br, bi)                                    \
{                                                                       \
    (cr) = (ar)*(br) - (ai)*(bi);                                       \
    (ci) = (ai)*(br) + (ar)*(bi);                                       \
}

#define C3Mul(dr, di, ar, ai, br, bi, cr, ci)                           \
{                                                                       \
    (dr) = (ar)*(br)*(cr) - (ai)*(bi)*(cr)                              \
            - (ai)*(br)*(ci) - (ar)*(bi)*(ci);                          \
    (di) = (ai)*(br)*(cr) + (ar)*(bi)*(cr)                              \
            + (ar)*(br)*(ci) - (ai)*(bi)*(ci);                          \
}

#endif // ! _complexmath_h
