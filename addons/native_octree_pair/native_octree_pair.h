#ifndef CBALLS_NATIVE_OCTREE_PAIR_H
#define CBALLS_NATIVE_OCTREE_PAIR_H

/* Include globaldefs.h before this header. */
typedef struct {
    int rank_count;
    bool (*task_owned)(struct cmdline_data *, INTEGER);
    bool (*publish)(struct cmdline_data *);
    int (*consensus)(struct cmdline_data *, int, const char *);
    int (*reduce_reals)(struct cmdline_data *, real *, size_t);
    int (*reduce_integers)(struct cmdline_data *, INTEGER *, size_t);
} cballs_native_pair_parallel;

typedef struct {
    /* BALLS4 historically correlates Weight*Kappa even without weights-norm. */
    bool weighted_signal;
    /* Its legacy no-one-ball option also disables pair-cell aggregation. */
    bool no_one_ball_is_exact;
    /* Preserve BALLS4's ordered auto-pair output with options=asymmetric. */
    bool honor_asymmetric;
    /* Preserve BALLS4's historical pair-density/CF normalization. */
    bool balls4_density_normalization;
    const cballs_native_pair_parallel *parallel;
} cballs_native_pair_policy;

int cballs_native_octree_pair_search(
        struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *,
        INTEGER, INTEGER *, int, int, const cballs_native_pair_policy *);

#endif /* CBALLS_NATIVE_OCTREE_PAIR_H */
