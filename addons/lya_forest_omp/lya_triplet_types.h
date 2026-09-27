/* Private per-pivot segment summaries. Each node has one forest and one
 * (radial, polar) bin; leaves retain the measured pixel geometry. */
typedef struct {
    size_t begin, end, left, right;
    REAL lo[3], hi[3];
    long double weight, weighted_delta;
} lya_segment;
