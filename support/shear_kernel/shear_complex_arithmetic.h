#ifndef CBALLS_SHEAR_COMPLEX_ARITHMETIC_H
#define CBALLS_SHEAR_COMPLEX_ARITHMETIC_H
/* Private numerical implementation: include after the engine types.
 * Arithmetic and iteration order are preserved; see engineering_contracts.rst. */
static shear_complex shear_make(real re, real im)
{
    shear_complex value;
    value.re = re;
    value.im = im;
    return value;
}

static shear_complex shear_add(shear_complex a, shear_complex b)
{
    return shear_make(a.re + b.re, a.im + b.im);
}

static shear_complex shear_mul(shear_complex a, shear_complex b);
static shear_complex shear_scale(shear_complex value, real scale);

static void shear_accumulate_complex(shear_complex *target,
                                     shear_complex increment)
{
#ifdef SHEAR_HAVE_VECTOR2
    shear_vector2 value;
    const shear_vector2 addend = {increment.re, increment.im};

    memcpy(&value, target, sizeof(value));
    value += addend;
    memcpy(target, &value, sizeof(value));
#else
    target->re += increment.re;
    target->im += increment.im;
#endif
}

static void shear_accumulate_complex_product(shear_complex *target,
                                             shear_complex left,
                                             shear_complex right)
{
#ifdef SHEAR_HAVE_VECTOR2
    shear_vector2 value;
    const shear_vector2 right_vector = {right.re, right.im};
    const shear_vector2 rotated_right = {-right.im, right.re};

    memcpy(&value, target, sizeof(value));
    value += right_vector*left.re + rotated_right*left.im;
    memcpy(target, &value, sizeof(value));
#else
    shear_accumulate_complex(target, shear_mul(left, right));
#endif
}

static void shear_accumulate_complex_scaled(shear_complex *target,
                                            shear_complex value,
                                            real scale)
{
#ifdef SHEAR_HAVE_VECTOR2
    shear_vector2 total;
    const shear_vector2 increment = {value.re, value.im};

    memcpy(&total, target, sizeof(total));
    total += increment*scale;
    memcpy(target, &total, sizeof(total));
#else
    shear_accumulate_complex(target, shear_scale(value, scale));
#endif
}

static shear_complex shear_sub(shear_complex a, shear_complex b)
{
    return shear_make(a.re - b.re, a.im - b.im);
}

static shear_complex shear_mul(shear_complex a, shear_complex b)
{
    return shear_make(a.re*b.re - a.im*b.im,
                      a.re*b.im + a.im*b.re);
}

static shear_complex shear_conj(shear_complex value)
{
    return shear_make(value.re, -value.im);
}

static shear_complex shear_scale(shear_complex value, real scale)
{
    return shear_make(value.re*scale, value.im*scale);
}

static real shear_abs2(shear_complex value)
{
    return value.re*value.re + value.im*value.im;
}


#endif
