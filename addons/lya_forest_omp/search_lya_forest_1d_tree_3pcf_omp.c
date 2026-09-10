/* Exact radial Lyman-alpha 3PCF using a balanced one-dimensional tree.
 *
 * For each pivot, the interval tree accumulates weighted neighbor moments in
 * signed radial-lag bins.  Ordered neighbor-pair moments are formed from their
 * outer products.  Pixels in the pivot forest are removed first, and ordered
 * pairs from every remaining individual forest are subtracted with a sparse
 * forest/bin accumulator.  The result is identical to enumerating triplets of
 * three distinct forests, without the quadratic neighbor-pair loop.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"
#include "lya_forest_parallel.h"

#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>

#ifndef LYA1D_TREE3_LEAF_SIZE
#define LYA1D_TREE3_LEAF_SIZE 8
#endif
#if LYA1D_TREE3_LEAF_SIZE < 1
#error "LYA1D_TREE3_LEAF_SIZE must be positive"
#endif

#define LYA1D_TREE3_NO_INDEX SIZE_MAX

typedef struct {
    size_t first;
    size_t end;
    size_t left;
    size_t right;
    REAL rmin;
    REAL rmax;
    long double sum_w;
    long double sum_w_delta;
    long double sum_w2;
    long double sum_w_delta2;
} lya1d_tree3_node;

typedef struct {
    INTEGER forest_id;
    size_t head;
    uint64_t generation;
} lya1d_tree3_forest_slot;

typedef struct {
    int bin;
    size_t next;
    uint64_t count;
    long double sum_w;
    long double sum_w_delta;
} lya1d_tree3_moment;

typedef struct {
    long double *hist_num;
    long double *hist_den;
    long double *all_w;
    long double *all_w_delta;
    long double *all_w2;
    long double *all_w_delta2;
    uint64_t *all_count;
    long double *pivot_w;
    long double *pivot_w_delta;
    long double *pivot_w2;
    long double *pivot_w_delta2;
    uint64_t *pivot_count;
    long double *same_num;
    long double *same_den;
    uint64_t *same_count;
    lya1d_tree3_forest_slot *forest_slots;
    size_t forest_capacity;
    size_t forest_count;
    uint64_t generation;
    lya1d_tree3_moment *moments;
    size_t moment_capacity;
    size_t moment_count;
    uint64_t candidate_visits;
    uint64_t ordered_triplets;
    uint64_t node_visits;
    uint64_t bulk_nodes;
    int failed;
} lya1d_tree3_worker;

local int lya1d_tree3_size_mul(size_t left, size_t right, size_t *result)
{
    if (left != 0 && right > SIZE_MAX / left) return FAILURE;
    *result = left * right;
    return SUCCESS;
}

local int lya1d_tree3_radial_order(const void *left_ptr,
                                    const void *right_ptr)
{
    bodyptr left = *(bodyptr const *)left_ptr;
    bodyptr right = *(bodyptr const *)right_ptr;

    if (LyaDistance(left) < LyaDistance(right)) return -1;
    if (LyaDistance(left) > LyaDistance(right)) return 1;
    if (Id(left) < Id(right)) return -1;
    if (Id(left) > Id(right)) return 1;
    return 0;
}

local int lya1d_tree3_signed_bin(REAL value, REAL maximum,
                                 int bins_per_side)
{
    int bin;
    int total_bins = 2 * bins_per_side;

    if (!isfinite(value) || !(value > -maximum && value < maximum))
        return -1;
    bin = (int)((value + maximum) / maximum * (REAL)bins_per_side);
    return bin >= 0 && bin < total_bins ? bin : -1;
}

local int lya1d_tree3_node_count(size_t count, size_t *result)
{
    size_t left_count;
    size_t right_count;
    size_t middle;

    if (count <= (size_t)LYA1D_TREE3_LEAF_SIZE) {
        *result = 1;
        return SUCCESS;
    }
    middle = count / 2;
    if (lya1d_tree3_node_count(middle, &left_count) == FAILURE
        || lya1d_tree3_node_count(count - middle, &right_count) == FAILURE
        || left_count > SIZE_MAX - right_count - 1)
        return FAILURE;
    *result = 1 + left_count + right_count;
    return SUCCESS;
}

local size_t lya1d_tree3_build(lya1d_tree3_node *nodes, size_t *used,
                               bodyptr *order, size_t first, size_t end)
{
    size_t node_index = (*used)++;
    lya1d_tree3_node *node = nodes + node_index;
    size_t i;

    memset(node, 0, sizeof(*node));
    node->first = first;
    node->end = end;
    node->left = LYA1D_TREE3_NO_INDEX;
    node->right = LYA1D_TREE3_NO_INDEX;
    node->rmin = LyaDistance(order[first]);
    node->rmax = LyaDistance(order[end - 1]);

    if (end - first <= (size_t)LYA1D_TREE3_LEAF_SIZE) {
        for (i = first; i < end; i++) {
            long double weight = (long double)Weight(order[i]);
            long double weighted_delta = weight * (long double)Kappa(order[i]);

            node->sum_w += weight;
            node->sum_w_delta += weighted_delta;
            node->sum_w2 += weight * weight;
            node->sum_w_delta2 += weighted_delta * weighted_delta;
        }
    } else {
        size_t middle = first + (end - first) / 2;
        lya1d_tree3_node *left;
        lya1d_tree3_node *right;

        node->left = lya1d_tree3_build(nodes, used, order, first, middle);
        node->right = lya1d_tree3_build(nodes, used, order, middle, end);
        left = nodes + node->left;
        right = nodes + node->right;
        node->sum_w = left->sum_w + right->sum_w;
        node->sum_w_delta = left->sum_w_delta + right->sum_w_delta;
        node->sum_w2 = left->sum_w2 + right->sum_w2;
        node->sum_w_delta2 = left->sum_w_delta2
                            + right->sum_w_delta2;
    }
    return node_index;
}

local uint64_t lya1d_tree3_hash(INTEGER forest_id)
{
    uint64_t value = (uint64_t)forest_id;

    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return value;
}

local int lya1d_tree3_worker_init(lya1d_tree3_worker *worker,
                                   size_t bins, size_t matrix_bins,
                                   ErrorMsg error_message)
{
    memset(worker, 0, sizeof(*worker));
    worker->forest_capacity = 16;
    worker->moment_capacity = 16;

    if (cballs_calloc_checked((void **)&worker->hist_num, matrix_bins,
                              sizeof(*worker->hist_num),
                              "radial tree 3PCF numerator", error_message,
                              _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->hist_den, matrix_bins,
                                 sizeof(*worker->hist_den),
                                 "radial tree 3PCF denominator", error_message,
                                 _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->all_w, bins,
                                 sizeof(*worker->all_w),
                                 "radial tree 3PCF weight moments",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->all_w_delta, bins,
                                 sizeof(*worker->all_w_delta),
                                 "radial tree 3PCF field moments",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->all_w2, bins,
                                 sizeof(*worker->all_w2),
                                 "radial tree 3PCF squared weights",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->all_w_delta2, bins,
                                 sizeof(*worker->all_w_delta2),
                                 "radial tree 3PCF squared fields",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->all_count, bins,
                                 sizeof(*worker->all_count),
                                 "radial tree 3PCF neighbor counts",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pivot_w, bins,
                                 sizeof(*worker->pivot_w),
                                 "radial tree pivot-forest weights",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pivot_w_delta, bins,
                                 sizeof(*worker->pivot_w_delta),
                                 "radial tree pivot-forest fields",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pivot_w2, bins,
                                 sizeof(*worker->pivot_w2),
                                 "radial tree pivot-forest squared weights",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pivot_w_delta2, bins,
                                 sizeof(*worker->pivot_w_delta2),
                                 "radial tree pivot-forest squared fields",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pivot_count, bins,
                                 sizeof(*worker->pivot_count),
                                 "radial tree pivot-forest counts",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->same_num, matrix_bins,
                                 sizeof(*worker->same_num),
                                 "radial tree same-forest numerator",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->same_den, matrix_bins,
                                 sizeof(*worker->same_den),
                                 "radial tree same-forest denominator",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->same_count, matrix_bins,
                                 sizeof(*worker->same_count),
                                 "radial tree same-forest counts",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->forest_slots,
                                 worker->forest_capacity,
                                 sizeof(*worker->forest_slots),
                                 "radial tree forest hash", error_message,
                                 _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->moments,
                                 worker->moment_capacity,
                                 sizeof(*worker->moments),
                                 "radial tree sparse moments", error_message,
                                 _ERRORMSGSIZE_) == FAILURE)
        goto fail;
    return SUCCESS;

fail:
    free(worker->hist_num);
    free(worker->hist_den);
    free(worker->all_w);
    free(worker->all_w_delta);
    free(worker->all_w2);
    free(worker->all_w_delta2);
    free(worker->all_count);
    free(worker->pivot_w);
    free(worker->pivot_w_delta);
    free(worker->pivot_w2);
    free(worker->pivot_w_delta2);
    free(worker->pivot_count);
    free(worker->same_num);
    free(worker->same_den);
    free(worker->same_count);
    free(worker->forest_slots);
    free(worker->moments);
    memset(worker, 0, sizeof(*worker));
    return FAILURE;
}

local void lya1d_tree3_worker_free(lya1d_tree3_worker *worker)
{
    free(worker->hist_num);
    free(worker->hist_den);
    free(worker->all_w);
    free(worker->all_w_delta);
    free(worker->all_w2);
    free(worker->all_w_delta2);
    free(worker->all_count);
    free(worker->pivot_w);
    free(worker->pivot_w_delta);
    free(worker->pivot_w2);
    free(worker->pivot_w_delta2);
    free(worker->pivot_count);
    free(worker->same_num);
    free(worker->same_den);
    free(worker->same_count);
    free(worker->forest_slots);
    free(worker->moments);
    memset(worker, 0, sizeof(*worker));
}

local int lya1d_tree3_grow_hash(lya1d_tree3_worker *worker)
{
    lya1d_tree3_forest_slot *replacement;
    size_t new_capacity;
    size_t i;

    if (worker->forest_capacity > SIZE_MAX / 2) return FAILURE;
    new_capacity = 2 * worker->forest_capacity;
    replacement = calloc(new_capacity, sizeof(*replacement));
    if (replacement == NULL) return FAILURE;

    for (i = 0; i < worker->forest_capacity; i++) {
        lya1d_tree3_forest_slot *old = worker->forest_slots + i;
        size_t slot;

        if (old->generation != worker->generation) continue;
        slot = (size_t)(lya1d_tree3_hash(old->forest_id)
                        & (uint64_t)(new_capacity - 1));
        while (replacement[slot].generation == worker->generation)
            slot = (slot + 1) & (new_capacity - 1);
        replacement[slot] = *old;
    }
    free(worker->forest_slots);
    worker->forest_slots = replacement;
    worker->forest_capacity = new_capacity;
    return SUCCESS;
}

local int lya1d_tree3_grow_moments(lya1d_tree3_worker *worker)
{
    lya1d_tree3_moment *replacement;
    size_t new_capacity;

    if (worker->moment_capacity > SIZE_MAX / 2
        || 2 * worker->moment_capacity > SIZE_MAX / sizeof(*replacement))
        return FAILURE;
    new_capacity = 2 * worker->moment_capacity;
    replacement = malloc(new_capacity * sizeof(*replacement));
    if (replacement == NULL) return FAILURE;
    memcpy(replacement, worker->moments,
           worker->moment_count * sizeof(*replacement));
    free(worker->moments);
    worker->moments = replacement;
    worker->moment_capacity = new_capacity;
    return SUCCESS;
}

local lya1d_tree3_forest_slot *lya1d_tree3_forest(
    lya1d_tree3_worker *worker, INTEGER forest_id)
{
    size_t slot;

    if (worker->forest_count >= worker->forest_capacity / 2
        && lya1d_tree3_grow_hash(worker) == FAILURE) {
        worker->failed = TRUE;
        return NULL;
    }
    slot = (size_t)(lya1d_tree3_hash(forest_id)
                    & (uint64_t)(worker->forest_capacity - 1));
    while (worker->forest_slots[slot].generation == worker->generation) {
        if (worker->forest_slots[slot].forest_id == forest_id)
            return worker->forest_slots + slot;
        slot = (slot + 1) & (worker->forest_capacity - 1);
    }
    worker->forest_slots[slot].generation = worker->generation;
    worker->forest_slots[slot].forest_id = forest_id;
    worker->forest_slots[slot].head = LYA1D_TREE3_NO_INDEX;
    worker->forest_count++;
    return worker->forest_slots + slot;
}

local int lya1d_tree3_add_same_forest(lya1d_tree3_worker *worker,
                                      INTEGER forest_id, int bin,
                                      long double weight,
                                      long double weighted_delta,
                                      size_t bins)
{
    lya1d_tree3_forest_slot *forest =
        lya1d_tree3_forest(worker, forest_id);
    size_t moment_index;
    size_t matching = LYA1D_TREE3_NO_INDEX;

    if (forest == NULL) return FAILURE;
    for (moment_index = forest->head;
         moment_index != LYA1D_TREE3_NO_INDEX;
         moment_index = worker->moments[moment_index].next) {
        lya1d_tree3_moment *moment = worker->moments + moment_index;
        size_t forward = (size_t)moment->bin * bins + (size_t)bin;
        size_t reverse = (size_t)bin * bins + (size_t)moment->bin;
        long double numerator = moment->sum_w_delta * weighted_delta;
        long double denominator = moment->sum_w * weight;

        if (forward == reverse) {
            if (moment->count > (UINT64_MAX - worker->same_count[forward]) / 2) {
                worker->failed = TRUE;
                return FAILURE;
            }
            worker->same_num[forward] += 2.0L * numerator;
            worker->same_den[forward] += 2.0L * denominator;
            worker->same_count[forward] += 2 * moment->count;
        } else {
            if (UINT64_MAX - worker->same_count[forward] < moment->count
                || UINT64_MAX - worker->same_count[reverse] < moment->count) {
                worker->failed = TRUE;
                return FAILURE;
            }
            worker->same_num[forward] += numerator;
            worker->same_den[forward] += denominator;
            worker->same_count[forward] += moment->count;
            worker->same_num[reverse] += numerator;
            worker->same_den[reverse] += denominator;
            worker->same_count[reverse] += moment->count;
        }
        if (moment->bin == bin) matching = moment_index;
    }

    if (matching == LYA1D_TREE3_NO_INDEX) {
        lya1d_tree3_moment *moment;

        if (worker->moment_count == worker->moment_capacity
            && lya1d_tree3_grow_moments(worker) == FAILURE) {
            worker->failed = TRUE;
            return FAILURE;
        }
        matching = worker->moment_count++;
        moment = worker->moments + matching;
        moment->bin = bin;
        moment->next = forest->head;
        moment->count = 0;
        moment->sum_w = 0.0L;
        moment->sum_w_delta = 0.0L;
        forest->head = matching;
    }
    if (worker->moments[matching].count == UINT64_MAX) {
        worker->failed = TRUE;
        return FAILURE;
    }
    worker->moments[matching].count++;
    worker->moments[matching].sum_w += weight;
    worker->moments[matching].sum_w_delta += weighted_delta;
    return SUCCESS;
}

local void lya1d_tree3_add_node(const lya1d_tree3_node *node, int bin,
                                lya1d_tree3_worker *worker)
{
    size_t count = node->end - node->first;

    worker->all_w[bin] += node->sum_w;
    worker->all_w_delta[bin] += node->sum_w_delta;
    worker->all_w2[bin] += node->sum_w2;
    worker->all_w_delta2[bin] += node->sum_w_delta2;
    worker->all_count[bin] += (uint64_t)count;
    worker->bulk_nodes++;
}

local int lya1d_tree3_query(const lya1d_tree3_node *nodes, size_t node_index,
                            bodyptr *order, REAL pivot_distance,
                            REAL maximum, int bins_per_side,
                            lya1d_tree3_worker *worker)
{
    const lya1d_tree3_node *node = nodes + node_index;
    REAL lag_min = node->rmin - pivot_distance;
    REAL lag_max = node->rmax - pivot_distance;
    int first_bin;
    int last_bin;
    size_t i;

    worker->node_visits++;
    if (lag_max <= -maximum || lag_min >= maximum) return SUCCESS;
    first_bin = lya1d_tree3_signed_bin(lag_min, maximum, bins_per_side);
    last_bin = lya1d_tree3_signed_bin(lag_max, maximum, bins_per_side);
    if (first_bin >= 0 && first_bin == last_bin) {
        lya1d_tree3_add_node(node, first_bin, worker);
        return SUCCESS;
    }

    if (node->left == LYA1D_TREE3_NO_INDEX) {
        for (i = node->first; i < node->end; i++) {
            bodyptr body = order[i];
            long double weight;
            long double weighted_delta;
            int bin = lya1d_tree3_signed_bin(
                LyaDistance(body) - pivot_distance, maximum, bins_per_side);

            if (bin < 0) continue;
            weight = (long double)Weight(body);
            weighted_delta = weight * (long double)Kappa(body);
            worker->all_w[bin] += weight;
            worker->all_w_delta[bin] += weighted_delta;
            worker->all_w2[bin] += weight * weight;
            worker->all_w_delta2[bin] += weighted_delta * weighted_delta;
            worker->all_count[bin]++;
        }
        return SUCCESS;
    }
    if (lya1d_tree3_query(nodes, node->left, order, pivot_distance,
                          maximum, bins_per_side, worker) == FAILURE
        || lya1d_tree3_query(nodes, node->right, order, pivot_distance,
                             maximum, bins_per_side, worker) == FAILURE)
        return FAILURE;
    return SUCCESS;
}

local int lya1d_tree3_count_product(uint64_t left, uint64_t right,
                                    uint64_t *result)
{
    if (left != 0 && right > UINT64_MAX / left) return FAILURE;
    *result = left * right;
    return SUCCESS;
}

local void lya1d_tree3_reset_pivot(lya1d_tree3_worker *worker,
                                   size_t bins, size_t matrix_bins)
{
    memset(worker->all_w, 0, bins * sizeof(*worker->all_w));
    memset(worker->all_w_delta, 0, bins * sizeof(*worker->all_w_delta));
    memset(worker->all_w2, 0, bins * sizeof(*worker->all_w2));
    memset(worker->all_w_delta2, 0, bins * sizeof(*worker->all_w_delta2));
    memset(worker->all_count, 0, bins * sizeof(*worker->all_count));
    memset(worker->pivot_w, 0, bins * sizeof(*worker->pivot_w));
    memset(worker->pivot_w_delta, 0,
           bins * sizeof(*worker->pivot_w_delta));
    memset(worker->pivot_w2, 0, bins * sizeof(*worker->pivot_w2));
    memset(worker->pivot_w_delta2, 0,
           bins * sizeof(*worker->pivot_w_delta2));
    memset(worker->pivot_count, 0, bins * sizeof(*worker->pivot_count));
    memset(worker->same_num, 0, matrix_bins * sizeof(*worker->same_num));
    memset(worker->same_den, 0, matrix_bins * sizeof(*worker->same_den));
    memset(worker->same_count, 0,
           matrix_bins * sizeof(*worker->same_count));
    worker->forest_count = 0;
    worker->moment_count = 0;
    worker->generation++;
    if (worker->generation == 0) {
        memset(worker->forest_slots, 0,
               worker->forest_capacity * sizeof(*worker->forest_slots));
        worker->generation = 1;
    }
}

local int lya1d_tree3_accumulate_pivot(struct cmdline_data *cmd,
                                       const lya1d_tree3_node *nodes,
                                       size_t root, bodyptr *order,
                                       size_t pivot_index,
                                       size_t first, size_t end,
                                       size_t bins, size_t matrix_bins,
                                       lya1d_tree3_worker *worker)
{
    bodyptr pivot = order[pivot_index];
    long double pivot_weight = (long double)Weight(pivot);
    long double pivot_field = pivot_weight * (long double)Kappa(pivot);
    size_t i;
    size_t b1;
    size_t b2;

    lya1d_tree3_reset_pivot(worker, bins, matrix_bins);
    if (lya1d_tree3_query(nodes, root, order, LyaDistance(pivot),
                          cmd->lya3RMax, cmd->lya3RBins, worker) == FAILURE)
        return FAILURE;

    for (i = first; i < end; i++) {
        bodyptr neighbor = order[i];
        long double weight = (long double)Weight(neighbor);
        long double weighted_delta = weight * (long double)Kappa(neighbor);
        int bin = lya1d_tree3_signed_bin(
            LyaDistance(neighbor) - LyaDistance(pivot),
            cmd->lya3RMax, cmd->lya3RBins);

        if (bin < 0) continue;
        if (i != pivot_index) worker->candidate_visits++;
        if (LyaForestId(neighbor) == LyaForestId(pivot)) {
            worker->pivot_w[bin] += weight;
            worker->pivot_w_delta[bin] += weighted_delta;
            worker->pivot_w2[bin] += weight * weight;
            worker->pivot_w_delta2[bin] += weighted_delta * weighted_delta;
            worker->pivot_count[bin]++;
        } else if (lya1d_tree3_add_same_forest(
                       worker, LyaForestId(neighbor), bin, weight,
                       weighted_delta, bins) == FAILURE) {
            return FAILURE;
        }
    }

    for (b1 = 0; b1 < bins; b1++) {
        if (worker->pivot_count[b1] > worker->all_count[b1]) {
            worker->failed = TRUE;
            return FAILURE;
        }
        worker->all_w[b1] -= worker->pivot_w[b1];
        worker->all_w_delta[b1] -= worker->pivot_w_delta[b1];
        worker->all_w2[b1] -= worker->pivot_w2[b1];
        worker->all_w_delta2[b1] -= worker->pivot_w_delta2[b1];
        worker->all_count[b1] -= worker->pivot_count[b1];
    }

    for (b1 = 0; b1 < bins; b1++) {
        for (b2 = 0; b2 < bins; b2++) {
            size_t index = b1 * bins + b2;
            uint64_t all_pairs;
            uint64_t cross_pairs;
            long double pair_num = worker->all_w_delta[b1]
                                 * worker->all_w_delta[b2];
            long double pair_den = worker->all_w[b1] * worker->all_w[b2];
            long double tolerance;

            if (lya1d_tree3_count_product(worker->all_count[b1],
                                          worker->all_count[b2],
                                          &all_pairs) == FAILURE) {
                worker->failed = TRUE;
                return FAILURE;
            }
            if (b1 == b2) {
                if (all_pairs < worker->all_count[b1]) {
                    worker->failed = TRUE;
                    return FAILURE;
                }
                all_pairs -= worker->all_count[b1];
                pair_num -= worker->all_w_delta2[b1];
                pair_den -= worker->all_w2[b1];
            }
            if (worker->same_count[index] > all_pairs) {
                worker->failed = TRUE;
                return FAILURE;
            }
            cross_pairs = all_pairs - worker->same_count[index];
            pair_num -= worker->same_num[index];
            pair_den -= worker->same_den[index];
            if (cross_pairs == 0) {
                pair_num = 0.0L;
                pair_den = 0.0L;
            }
            tolerance = 128.0L * LDBL_EPSILON
                      * (fabsl(worker->all_w[b1] * worker->all_w[b2])
                         + fabsl(worker->same_den[index])
                         + fabsl(b1 == b2 ? worker->all_w2[b1] : 0.0L)
                         + 1.0L);
            if (pair_den < 0.0L) {
                if (-pair_den <= tolerance) {
                    pair_den = 0.0L;
                } else {
                    worker->failed = TRUE;
                    return FAILURE;
                }
            }
            if (pair_den == 0.0L) pair_num = 0.0L;
            if (UINT64_MAX - worker->ordered_triplets < cross_pairs) {
                worker->failed = TRUE;
                return FAILURE;
            }
            worker->ordered_triplets += cross_pairs;
            worker->hist_num[index] += pivot_field * pair_num;
            worker->hist_den[index] += pivot_weight * pair_den;
        }
    }
    return SUCCESS;
}

local void lya1d_tree3_reset_block(lya1d_tree3_worker *worker,
                                   size_t matrix_bins)
{
    memset(worker->hist_num, 0, matrix_bins * sizeof(*worker->hist_num));
    memset(worker->hist_den, 0, matrix_bins * sizeof(*worker->hist_den));
    worker->candidate_visits = 0;
    worker->ordered_triplets = 0;
    worker->node_visits = 0;
    worker->bulk_nodes = 0;
    worker->failed = FALSE;
}

local int lya1d_tree3_commit(lya1d_tree3_worker *worker,
                             long double *num, long double *den,
                             size_t matrix_bins, uint64_t counters[4])
{
    size_t i;
    uint64_t block_counters[4] = {
        worker->candidate_visits, worker->ordered_triplets,
        worker->node_visits, worker->bulk_nodes
    };

    for (i = 0; i < 4; i++)
        if (UINT64_MAX - counters[i] < block_counters[i]) return FAILURE;
    for (i = 0; i < matrix_bins; i++) {
        num[i] += worker->hist_num[i];
        den[i] += worker->hist_den[i];
    }
    for (i = 0; i < 4; i++) counters[i] += block_counters[i];
    return SUCCESS;
}

local int lya1d_tree3_write(struct cmdline_data *cmd,
                            struct global_data *gd,
                            const long double *num,
                            const long double *den,
                            uint64_t ordered_triplets)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream;
    int bins_per_side = cmd->lya3RBins;
    int total_bins = 2 * bins_per_side;
    int output_empty = cballs_opt_lya_output_empty_bins(cmd);
    int b1;
    int b2;
    int write_failed;

    if (format_checked(path, sizeof(path), "radial tree 3PCF path",
                       "%s_lya1d%s", gd->fpfnamehistZetaMFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 3PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open radial tree 3PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }

    fprintf(stream, "# Exact radial Lyman-alpha forest 3PCF, interval tree\n");
    fprintf(stream, "# distinct-forest ordered triplets: %" PRIu64 "\n",
            ordered_triplets);
    fprintf(stream, "# transverse separation is ignored; lags are signed about the pivot\n");
    fprintf(stream, "# zero-denominator policy: zeta=0; empty bins are %s\n",
            output_empty ? "included" : "omitted");
    fprintf(stream, "# columns: b1 b2 lag1 lag2 zeta numerator denominator\n");
    for (b1 = 0; b1 < total_bins; b1++) {
        for (b2 = 0; b2 < total_bins; b2++) {
            size_t index = (size_t)b1 * (size_t)total_bins + (size_t)b2;
            long double zeta;
            REAL lag1;
            REAL lag2;

            if (!output_empty && den[index] == 0.0L) continue;
            lag1 = -cmd->lya3RMax
                 + ((REAL)b1 + 0.5) * cmd->lya3RMax / (REAL)bins_per_side;
            lag2 = -cmd->lya3RMax
                 + ((REAL)b2 + 0.5) * cmd->lya3RMax / (REAL)bins_per_side;
            zeta = den[index] == 0.0L ? 0.0L : num[index] / den[index];
            fprintf(stream, "%d %d %.17g %.17g %.17g %.17g %.17g\n",
                    b1, b2, lag1, lag2, (double)zeta,
                    (double)num[index], (double)den[index]);
        }
    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing radial tree 3PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

global int searchcalc_lya_forest_1d_tree_3pcf_omp(
    struct cmdline_data *cmd, struct global_data *gd,
    bodyptr table, INTEGER nbody)
{
    bodyptr *order = NULL;
    size_t *first = NULL;
    size_t *end = NULL;
    lya1d_tree3_node *nodes = NULL;
    long double *num = NULL;
    long double *den = NULL;
    size_t count = (size_t)nbody;
    size_t active_count = 0;
    size_t bins = 0;
    size_t matrix_bins = 0;
    size_t node_capacity = 0;
    size_t nodes_used = 0;
    size_t root = 0;
    size_t block_count = 0;
    size_t left_cursor = 0;
    size_t right_cursor = 0;
    size_t i;
    uint64_t counters[4] = {0, 0, 0, 0};
    int allocation_failed = FALSE;
    int runtime_failed = FALSE;
    ErrorMsg worker_error = "";
    int status = FAILURE;
    double cpustart = CPUTIME;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial tree Ly-alpha 3PCF requires 3D catalog storage");
    return FAILURE;
#endif
#ifndef OPENMPCODE
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial tree Ly-alpha 3PCF requires OPENMPMACHINE=1");
    return FAILURE;
#endif

    if (nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree Ly-alpha 3PCF received an empty catalog");
        return FAILURE;
    }
    if (cmd->lya3RBins < 1 || cmd->lya3RBins > INT_MAX / 2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree Ly-alpha 3PCF bin count is invalid");
        return FAILURE;
    }
    bins = 2 * (size_t)cmd->lya3RBins;
    if (lya1d_tree3_size_mul(bins, bins, &matrix_bins) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree Ly-alpha 3PCF histogram dimensions overflow");
        return FAILURE;
    }

    if (cballs_calloc_checked((void **)&order, count, sizeof(*order),
                              "radial tree 3PCF sorted index",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    for (i = 0; i < count; i++) {
        bodyptr body = table + i;

        if (Update(body) == FALSE || Mask(body) != MASK_NODE_VALID) continue;
        order[active_count++] = body;
    }
    if (cballs_calloc_checked((void **)&num, matrix_bins, sizeof(*num),
                              "global radial tree 3PCF numerator",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&den, matrix_bins, sizeof(*den),
                                 "global radial tree 3PCF denominator",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    if (active_count == 0) {
        status = SUCCESS;
        goto setup_done;
    }

    qsort(order, active_count, sizeof(*order), lya1d_tree3_radial_order);
    if (cballs_calloc_checked((void **)&first, active_count, sizeof(*first),
                              "radial tree 3PCF lower limits",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&end, active_count, sizeof(*end),
                                 "radial tree 3PCF upper limits",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    for (i = 0; i < active_count; i++) {
        REAL pivot_distance = LyaDistance(order[i]);

        while (left_cursor < i
               && pivot_distance - LyaDistance(order[left_cursor])
                  >= cmd->lya3RMax)
            left_cursor++;
        if (right_cursor < i + 1) right_cursor = i + 1;
        while (right_cursor < active_count
               && LyaDistance(order[right_cursor]) - pivot_distance
                  < cmd->lya3RMax)
            right_cursor++;
        first[i] = left_cursor;
        end[i] = right_cursor;
    }
    if (lya1d_tree3_node_count(active_count, &node_capacity) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 3PCF node count overflows size_t");
        goto setup_done;
    }
    if (cballs_calloc_checked((void **)&nodes, node_capacity,
                              sizeof(*nodes), "radial 3PCF interval tree",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    root = lya1d_tree3_build(nodes, &nodes_used, order, 0, active_count);
    if (nodes_used != node_capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 3PCF construction produced an inconsistent node count");
        goto setup_done;
    }
    block_count = 1 + (active_count - 1)
                      / (size_t)LYA1D_OMP_PIVOT_BLOCK_SIZE;
    status = SUCCESS;

setup_done:
    status = lya_parallel_consensus(cmd, status,
                                    "Ly-alpha radial-tree 3PCF setup");
    if (status == FAILURE) goto cleanup;
    if (active_count == 0) goto postprocess;

    ThreadCount(cmd, gd, (INTEGER)active_count, 0);
    const size_t first_task = lya_parallel_first(cmd);
    const size_t task_stride = lya_parallel_stride(cmd);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: exact 1D interval-tree 3PCF; pixels=%zu nodes=%zu "
        "bins=%zu pivot_blocks=%zu\n",
        cmd->searchMethod, active_count, nodes_used, bins, block_count);

#pragma omp parallel shared(allocation_failed,runtime_failed,worker_error,num,den,counters)
    {
        lya1d_tree3_worker worker;
        ErrorMsg local_error = "";
        int worker_ready = lya1d_tree3_worker_init(
            &worker, bins, matrix_bins, local_error) == SUCCESS;
        size_t block_index;

        if (!worker_ready) {
#pragma omp critical(lya1d_tree3_failure)
            {
                if (!allocation_failed)
                    snprintf(worker_error, sizeof(worker_error), "%s",
                             local_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier
#pragma omp for schedule(static,1) ordered
        for (block_index = first_task; block_index < block_count;
             block_index += task_stride) {
            size_t block_first = block_index
                               * (size_t)LYA1D_OMP_PIVOT_BLOCK_SIZE;
            size_t block_end = MIN(
                block_first + (size_t)LYA1D_OMP_PIVOT_BLOCK_SIZE,
                active_count);
            size_t pivot_index;

            if (worker_ready) lya1d_tree3_reset_block(&worker, matrix_bins);
            if (worker_ready && !allocation_failed) {
                for (pivot_index = block_first;
                     pivot_index < block_end; pivot_index++) {
                    if (lya1d_tree3_accumulate_pivot(
                            cmd, nodes, root, order, pivot_index,
                            first[pivot_index], end[pivot_index], bins,
                            matrix_bins, &worker) == FAILURE)
                        break;
                }
                if (worker.failed) {
#pragma omp critical(lya1d_tree3_runtime_failure)
                    {
                        if (!runtime_failed)
                            snprintf(worker_error, sizeof(worker_error),
                                     "radial tree 3PCF worker allocation, counter, or numerical failure");
                        runtime_failed = TRUE;
                    }
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !allocation_failed && !worker.failed
                    && lya1d_tree3_commit(&worker, num, den, matrix_bins,
                                          counters) == FAILURE) {
#pragma omp critical(lya1d_tree3_runtime_failure)
                    {
                        if (!runtime_failed)
                            snprintf(worker_error, sizeof(worker_error),
                                     "radial tree 3PCF global counter overflow");
                        runtime_failed = TRUE;
                    }
                }
            }
        }
        if (worker_ready) lya1d_tree3_worker_free(&worker);
    }

postprocess:
    if (allocation_failed || runtime_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 3PCF worker failed: %.2000s",
                 worker_error[0] != '\0' ? worker_error : "unknown error");
    }
    status = lya_parallel_consensus(
        cmd, allocation_failed || runtime_failed ? FAILURE : SUCCESS,
        "Ly-alpha radial-tree 3PCF workers");
    if (status == FAILURE) goto cleanup;
    if (lya_parallel_reduce_long_doubles(cmd, num, matrix_bins) == FAILURE
        || lya_parallel_reduce_long_doubles(cmd, den, matrix_bins) == FAILURE
        || lya_parallel_reduce_uint64(cmd, counters, 4) == FAILURE) {
        status = FAILURE;
        goto cleanup;
    }
    if (!lya_parallel_publish(cmd)) {
        status = SUCCESS;
        goto publication;
    }
    status = FAILURE;

#ifdef LONGINT
    if (counters[0] > (uint64_t)LONG_MAX
        || gd->nbbcalc > LONG_MAX - (INTEGER)counters[0]) {
#else
    if (counters[0] > (uint64_t)INT_MAX
        || gd->nbbcalc > INT_MAX - (INTEGER)counters[0]) {
#endif
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 3PCF candidate count exceeds the INTEGER range");
        goto publication;
    }
    gd->nbbcalc += (INTEGER)counters[0];
    if (!cballs_opt_no_out_hist(cmd)
        && lya1d_tree3_write(cmd, gd, num, den, counters[1]) == FAILURE)
        goto publication;

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: radial_candidates=%" PRIu64 " ordered_triplets=%" PRIu64
        " node_visits=%" PRIu64 " bulk_nodes=%" PRIu64 " CPU=%g\n",
        cmd->searchMethod, counters[0], counters[1], counters[2], counters[3],
        gd->cpusearch);
    status = SUCCESS;

publication:
    status = lya_parallel_consensus(cmd, status,
                                    "Ly-alpha radial-tree 3PCF output");
cleanup:
    gd->cpusearch = CPUTIME - cpustart;
    free(order);
    free(first);
    free(end);
    free(nodes);
    free(num);
    free(den);
    return status;
}
