/* Private per-pivot segment summaries. Each node has one forest and one
 * (radial, polar) bin; leaves retain the measured pixel geometry. */
typedef struct {
    size_t begin, end, left, right;
    REAL lo[3], hi[3];
    long double weight, weighted_delta;
} lya_segment;

/* Per-pivot hierarchy of forest-segment moments, independent of cell kernels. */
typedef struct {
    lya_segment bounds;
    size_t source,left,right,count,bin,first,second;
    INTEGER forest_lo,forest_hi;
    uint64_t forests[4];
} lya_los_moment;
typedef struct {size_t source,bin; REAL angle;} lya_los_moment_entry;
