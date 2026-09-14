#ifndef DUAL_NODE_SLEEF_LOG_H
#define DUAL_NODE_SLEEF_LOG_H

#include <stddef.h>

#include <sleef.h>

#if defined(__AVX2__) || defined(__SSE2__)
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

static inline int dual_node_sleef_double_lanes(void)
{
#if defined(__AVX2__)
    return 4;
#elif defined(__SSE2__) || defined(__aarch64__)
    return 2;
#else
    return 1;
#endif
}

static inline int dual_node_sleef_float_lanes(void)
{
#if defined(__AVX2__)
    return 8;
#elif defined(__SSE2__) || defined(__aarch64__)
    return 4;
#else
    return 1;
#endif
}

static inline void dual_node_sleef_log_double(
        double *output, const double *input, size_t count)
{
    size_t i = 0;

#if defined(__AVX2__)
    for (; i + 4 <= count; i += 4) {
        const __m256d value = _mm256_loadu_pd(input + i);
        _mm256_storeu_pd(output + i, Sleef_logd4_u10avx2(value));
    }
#elif defined(__SSE2__)
    for (; i + 2 <= count; i += 2) {
        const __m128d value = _mm_loadu_pd(input + i);
        _mm_storeu_pd(output + i, Sleef_logd2_u10sse2(value));
    }
#elif defined(__aarch64__)
    for (; i + 2 <= count; i += 2) {
        const float64x2_t value = vld1q_f64(input + i);
        vst1q_f64(output + i, Sleef_logd2_u10advsimd(value));
    }
#endif
    for (; i < count; i++) output[i] = Sleef_log_u10(input[i]);
}

static inline void dual_node_sleef_log_float(
        float *output, const float *input, size_t count)
{
    size_t i = 0;

#if defined(__AVX2__)
    for (; i + 8 <= count; i += 8) {
        const __m256 value = _mm256_loadu_ps(input + i);
        _mm256_storeu_ps(output + i, Sleef_logf8_u10avx2(value));
    }
#elif defined(__SSE2__)
    for (; i + 4 <= count; i += 4) {
        const __m128 value = _mm_loadu_ps(input + i);
        _mm_storeu_ps(output + i, Sleef_logf4_u10sse2(value));
    }
#elif defined(__aarch64__)
    for (; i + 4 <= count; i += 4) {
        const float32x4_t value = vld1q_f32(input + i);
        vst1q_f32(output + i, Sleef_logf4_u10advsimd(value));
    }
#endif
    for (; i < count; i++) output[i] = Sleef_logf_u10(input[i]);
}

#endif
