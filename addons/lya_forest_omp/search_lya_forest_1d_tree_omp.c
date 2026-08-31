/* Exact radial Lyman-alpha 2PCF using a balanced one-dimensional tree.
 *
 * The tree aggregates a node pair only when the complete separation interval
 * belongs to one histogram bin.  Same-forest pairs are accumulated separately
 * and subtracted in long-double precision, leaving the required cross-forest
 * estimator.  Fixed task blocks and ordered commits make the result independent
 * of the OpenMP thread count.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"
#include "lya_forest_parallel.h"

#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>

#ifndef LYA1D_TREE_LEAF_SIZE
#define LYA1D_TREE_LEAF_SIZE 4
#endif
#ifndef LYA1D_TREE_TASK_PAIR_LIMIT
#define LYA1D_TREE_TASK_PAIR_LIMIT 262144
#endif
#ifndef LYA1D_TREE_REDUCTION_BLOCK_SIZE
#define LYA1D_TREE_REDUCTION_BLOCK_SIZE 16
#endif
#ifndef LYA1D_TREE_FOREST_PIVOT_CHUNK
#define LYA1D_TREE_FOREST_PIVOT_CHUNK 256
#endif

#if LYA1D_TREE_LEAF_SIZE < 1
#error "LYA1D_TREE_LEAF_SIZE must be positive"
#endif
#if LYA1D_TREE_TASK_PAIR_LIMIT < 1
#error "LYA1D_TREE_TASK_PAIR_LIMIT must be positive"
#endif
#if LYA1D_TREE_REDUCTION_BLOCK_SIZE < 1
#error "LYA1D_TREE_REDUCTION_BLOCK_SIZE must be positive"
#endif
#if LYA1D_TREE_FOREST_PIVOT_CHUNK < 1
#error "LYA1D_TREE_FOREST_PIVOT_CHUNK must be positive"
#endif

#define LYA1D_TREE_NO_CHILD SIZE_MAX
#define LYA1D_TREE_UNRESOLVED_BIN (-1)
#define LYA1D_TREE_REJECTED_PAIR (-2)

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
} lya1d_tree_node;

typedef struct {
    size_t first_node;
    size_t second_node;
} lya1d_tree_task;

typedef struct {
    size_t forest_end;
    size_t pivot_first;
    size_t pivot_end;
} lya1d_forest_task;

typedef struct {
    long double *num;
    long double *den;
    uint64_t *pairs;
    size_t *touched;
    unsigned char *marked;
    size_t touched_count;
    uint64_t node_visits;
    uint64_t bulk_pairs;
    int failed;
} lya1d_tree_worker;

local int lya1d_tree_radial_order(const void *left_ptr,
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

local int lya1d_tree_forest_order(const void *left_ptr,
                                  const void *right_ptr)
{
    bodyptr left = *(bodyptr const *)left_ptr;
    bodyptr right = *(bodyptr const *)right_ptr;

    if (LyaForestId(left) < LyaForestId(right)) return -1;
    if (LyaForestId(left) > LyaForestId(right)) return 1;
    return lya1d_tree_radial_order(left_ptr, right_ptr);
}

local int lya1d_tree_positive_bin(REAL value, REAL maximum, int bins)
{
    int bin;

    if (!isfinite(value) || value < 0.0 || value >= maximum) return -1;
    bin = (int)(value / maximum * (REAL)bins);
    return bin >= 0 && bin < bins ? bin : -1;
}

local int lya1d_tree_node_count(size_t count, size_t *result)
{
    size_t left_count;
    size_t right_count;
    size_t middle;

    if (count <= (size_t)LYA1D_TREE_LEAF_SIZE) {
        *result = 1;
        return SUCCESS;
    }
    middle = count / 2;
    if (lya1d_tree_node_count(middle, &left_count) == FAILURE
        || lya1d_tree_node_count(count - middle, &right_count) == FAILURE
        || left_count > SIZE_MAX - right_count - 1)
        return FAILURE;
    *result = 1 + left_count + right_count;
    return SUCCESS;
}

local size_t lya1d_tree_build(lya1d_tree_node *nodes, size_t *used,
                              bodyptr *order, size_t first, size_t end)
{
    size_t node_index = (*used)++;
    lya1d_tree_node *node = nodes + node_index;
    size_t i;

    memset(node, 0, sizeof(*node));
    node->first = first;
    node->end = end;
    node->left = LYA1D_TREE_NO_CHILD;
    node->right = LYA1D_TREE_NO_CHILD;
    node->rmin = LyaDistance(order[first]);
    node->rmax = LyaDistance(order[end - 1]);

    if (end - first <= (size_t)LYA1D_TREE_LEAF_SIZE) {
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
        lya1d_tree_node *left;
        lya1d_tree_node *right;

        node->left = lya1d_tree_build(nodes, used, order, first, middle);
        node->right = lya1d_tree_build(nodes, used, order, middle, end);
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

local int lya1d_tree_pair_bin(const lya1d_tree_node *first,
                              const lya1d_tree_node *second,
                              REAL maximum, int bins)
{
    REAL dmin;
    REAL dmax;
    int first_bin;
    int last_bin;

    if (first == second) {
        dmin = 0.0;
        dmax = first->rmax - first->rmin;
    } else {
        dmin = second->rmin - first->rmax;
        dmax = second->rmax - first->rmin;
        if (dmin < 0.0) dmin = 0.0;
    }
    if (dmin >= maximum) return LYA1D_TREE_REJECTED_PAIR;
    first_bin = lya1d_tree_positive_bin(dmin, maximum, bins);
    last_bin = lya1d_tree_positive_bin(dmax, maximum, bins);
    if (first_bin >= 0 && first_bin == last_bin) return first_bin;
    return LYA1D_TREE_UNRESOLVED_BIN;
}

local int lya1d_tree_is_leaf(const lya1d_tree_node *node)
{
    return node->left == LYA1D_TREE_NO_CHILD;
}

local int lya1d_tree_small_task(const lya1d_tree_node *first,
                                const lya1d_tree_node *second)
{
    size_t first_count = first->end - first->first;
    size_t second_count = second->end - second->first;
    size_t limit = (size_t)LYA1D_TREE_TASK_PAIR_LIMIT;

    if (first == second) second_count = first_count;
    return second_count == 0 || first_count <= limit / second_count;
}

local int lya1d_tree_append_task(lya1d_tree_task *tasks, size_t capacity,
                                 size_t *used, size_t first, size_t second)
{
    if (*used == SIZE_MAX || (tasks != NULL && *used >= capacity))
        return FAILURE;
    if (tasks != NULL) {
        tasks[*used].first_node = first;
        tasks[*used].second_node = second;
    }
    (*used)++;
    return SUCCESS;
}

local int lya1d_tree_collect_tasks(const lya1d_tree_node *nodes,
                                   size_t first_index, size_t second_index,
                                   REAL maximum, int bins,
                                   lya1d_tree_task *tasks, size_t capacity,
                                   size_t *used)
{
    const lya1d_tree_node *first = nodes + first_index;
    const lya1d_tree_node *second = nodes + second_index;
    int pair_bin = lya1d_tree_pair_bin(first, second, maximum, bins);
    int first_leaf = lya1d_tree_is_leaf(first);
    int second_leaf = lya1d_tree_is_leaf(second);
    size_t first_count = first->end - first->first;
    size_t second_count = second->end - second->first;
    REAL first_span = first->rmax - first->rmin;
    REAL second_span = second->rmax - second->rmin;

    if (pair_bin == LYA1D_TREE_REJECTED_PAIR) return SUCCESS;
    if (pair_bin >= 0 || lya1d_tree_small_task(first, second)
        || (first_leaf && second_leaf))
        return lya1d_tree_append_task(tasks, capacity, used,
                                      first_index, second_index);

    if (first_index == second_index) {
        if (lya1d_tree_collect_tasks(nodes, first->left, first->left,
                                     maximum, bins, tasks, capacity, used)
            == FAILURE
            || lya1d_tree_collect_tasks(nodes, first->left, first->right,
                                        maximum, bins, tasks, capacity, used)
               == FAILURE
            || lya1d_tree_collect_tasks(nodes, first->right, first->right,
                                        maximum, bins, tasks, capacity, used)
               == FAILURE)
            return FAILURE;
        return SUCCESS;
    }

    if (!first_leaf
        && (second_leaf || first_count > second_count
            || (first_count == second_count && first_span >= second_span))) {
        if (lya1d_tree_collect_tasks(nodes, first->left, second_index,
                                     maximum, bins, tasks, capacity, used)
            == FAILURE
            || lya1d_tree_collect_tasks(nodes, first->right, second_index,
                                        maximum, bins, tasks, capacity, used)
               == FAILURE)
            return FAILURE;
    } else {
        if (lya1d_tree_collect_tasks(nodes, first_index, second->left,
                                     maximum, bins, tasks, capacity, used)
            == FAILURE
            || lya1d_tree_collect_tasks(nodes, first_index, second->right,
                                        maximum, bins, tasks, capacity, used)
               == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

local int lya1d_tree_worker_init(lya1d_tree_worker *worker, size_t bins,
                                 ErrorMsg error_message)
{
    memset(worker, 0, sizeof(*worker));
    if (cballs_calloc_checked((void **)&worker->num, bins,
                              sizeof(*worker->num),
                              "radial tree worker numerator",
                              error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->den, bins,
                                 sizeof(*worker->den),
                                 "radial tree worker denominator",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->pairs, bins,
                                 sizeof(*worker->pairs),
                                 "radial tree worker pair counts",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->touched, bins,
                                 sizeof(*worker->touched),
                                 "radial tree worker touched bins",
                                 error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_calloc_checked((void **)&worker->marked, bins,
                                 sizeof(*worker->marked),
                                 "radial tree worker touched flags",
                                 error_message, _ERRORMSGSIZE_) == FAILURE)
        goto fail;
    return SUCCESS;

fail:
    free(worker->num);
    free(worker->den);
    free(worker->pairs);
    free(worker->touched);
    free(worker->marked);
    memset(worker, 0, sizeof(*worker));
    return FAILURE;
}

local void lya1d_tree_worker_free(lya1d_tree_worker *worker)
{
    free(worker->num);
    free(worker->den);
    free(worker->pairs);
    free(worker->touched);
    free(worker->marked);
    memset(worker, 0, sizeof(*worker));
}

local int lya1d_tree_count_product(size_t first, size_t second,
                                   int same_node, uint64_t *result)
{
    if (same_node) {
        if (first < 2) {
            *result = 0;
            return SUCCESS;
        }
        if (first > UINT64_MAX / (first - 1)) return FAILURE;
        *result = (uint64_t)(first * (first - 1) / 2);
        return SUCCESS;
    }
    if (first != 0 && second > UINT64_MAX / first) return FAILURE;
    *result = (uint64_t)(first * second);
    return SUCCESS;
}

local int lya1d_tree_worker_add(lya1d_tree_worker *worker, int bin,
                                long double numerator,
                                long double denominator, uint64_t pairs)
{
    if (UINT64_MAX - worker->pairs[bin] < pairs) {
        worker->failed = TRUE;
        return FAILURE;
    }
    if (!worker->marked[bin]) {
        worker->marked[bin] = 1;
        worker->touched[worker->touched_count++] = (size_t)bin;
    }
    worker->num[bin] += numerator;
    worker->den[bin] += denominator;
    worker->pairs[bin] += pairs;
    return SUCCESS;
}

local int lya1d_tree_add_nodes(const lya1d_tree_node *first,
                               const lya1d_tree_node *second, int bin,
                               lya1d_tree_worker *worker)
{
    size_t first_count = first->end - first->first;
    size_t second_count = second->end - second->first;
    uint64_t pairs;
    long double numerator;
    long double denominator;

    if (lya1d_tree_count_product(first_count, second_count, first == second,
                                 &pairs) == FAILURE) {
        worker->failed = TRUE;
        return FAILURE;
    }
    if (first == second) {
        numerator = 0.5L * (first->sum_w_delta * first->sum_w_delta
                            - first->sum_w_delta2);
        denominator = 0.5L * (first->sum_w * first->sum_w
                              - first->sum_w2);
    } else {
        numerator = first->sum_w_delta * second->sum_w_delta;
        denominator = first->sum_w * second->sum_w;
    }
    if (lya1d_tree_worker_add(worker, bin, numerator, denominator, pairs)
        == FAILURE)
        return FAILURE;
    if (UINT64_MAX - worker->bulk_pairs < pairs) {
        worker->failed = TRUE;
        return FAILURE;
    }
    worker->bulk_pairs += pairs;
    return SUCCESS;
}

local int lya1d_tree_add_direct(bodyptr first, bodyptr second, REAL maximum,
                                int bins, lya1d_tree_worker *worker)
{
    REAL separation = rabs(LyaDistance(second) - LyaDistance(first));
    int bin = lya1d_tree_positive_bin(separation, maximum, bins);
    long double first_weight;
    long double second_weight;

    if (bin < 0) return SUCCESS;
    first_weight = (long double)Weight(first);
    second_weight = (long double)Weight(second);
    return lya1d_tree_worker_add(
        worker, bin,
        first_weight * (long double)Kappa(first)
            * second_weight * (long double)Kappa(second),
        first_weight * second_weight, 1);
}

local int lya1d_tree_walk(const lya1d_tree_node *nodes, bodyptr *order,
                          size_t first_index, size_t second_index,
                          REAL maximum, int bins,
                          lya1d_tree_worker *worker)
{
    const lya1d_tree_node *first = nodes + first_index;
    const lya1d_tree_node *second = nodes + second_index;
    int pair_bin = lya1d_tree_pair_bin(first, second, maximum, bins);
    int first_leaf = lya1d_tree_is_leaf(first);
    int second_leaf = lya1d_tree_is_leaf(second);
    size_t first_count = first->end - first->first;
    size_t second_count = second->end - second->first;
    REAL first_span = first->rmax - first->rmin;
    REAL second_span = second->rmax - second->rmin;
    size_t i;
    size_t j;

    worker->node_visits++;
    if (pair_bin == LYA1D_TREE_REJECTED_PAIR) return SUCCESS;
    if (pair_bin >= 0)
        return lya1d_tree_add_nodes(first, second, pair_bin, worker);

    if (first_leaf && second_leaf) {
        if (first_index == second_index) {
            for (i = first->first; i < first->end; i++)
                for (j = i + 1; j < first->end; j++)
                    if (lya1d_tree_add_direct(order[i], order[j], maximum,
                                              bins, worker) == FAILURE)
                        return FAILURE;
        } else {
            for (i = first->first; i < first->end; i++)
                for (j = second->first; j < second->end; j++)
                    if (lya1d_tree_add_direct(order[i], order[j], maximum,
                                              bins, worker) == FAILURE)
                        return FAILURE;
        }
        return SUCCESS;
    }

    if (first_index == second_index) {
        if (lya1d_tree_walk(nodes, order, first->left, first->left,
                            maximum, bins, worker) == FAILURE
            || lya1d_tree_walk(nodes, order, first->left, first->right,
                               maximum, bins, worker) == FAILURE
            || lya1d_tree_walk(nodes, order, first->right, first->right,
                               maximum, bins, worker) == FAILURE)
            return FAILURE;
        return SUCCESS;
    }

    if (!first_leaf
        && (second_leaf || first_count > second_count
            || (first_count == second_count && first_span >= second_span))) {
        if (lya1d_tree_walk(nodes, order, first->left, second_index,
                            maximum, bins, worker) == FAILURE
            || lya1d_tree_walk(nodes, order, first->right, second_index,
                               maximum, bins, worker) == FAILURE)
            return FAILURE;
    } else {
        if (lya1d_tree_walk(nodes, order, first_index, second->left,
                            maximum, bins, worker) == FAILURE
            || lya1d_tree_walk(nodes, order, first_index, second->right,
                               maximum, bins, worker) == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

local void lya1d_tree_worker_commit(lya1d_tree_worker *worker,
                                    long double *num, long double *den,
                                    uint64_t *pairs)
{
    size_t i;

    for (i = 0; i < worker->touched_count; i++) {
        size_t bin = worker->touched[i];

        num[bin] += worker->num[bin];
        den[bin] += worker->den[bin];
        pairs[bin] += worker->pairs[bin];
        worker->num[bin] = 0.0L;
        worker->den[bin] = 0.0L;
        worker->pairs[bin] = 0;
        worker->marked[bin] = 0;
    }
    worker->touched_count = 0;
}

local int lya1d_tree_add_forest_task(const lya1d_forest_task *task,
                                     bodyptr *forest_order, REAL maximum,
                                     int bins, lya1d_tree_worker *worker)
{
    size_t pivot;

    for (pivot = task->pivot_first; pivot < task->pivot_end; pivot++) {
        bodyptr first = forest_order[pivot];
        size_t second;

        for (second = pivot + 1; second < task->forest_end; second++) {
            bodyptr next = forest_order[second];

            if (LyaDistance(next) - LyaDistance(first) >= maximum) break;
            if (lya1d_tree_add_direct(first, next, maximum, bins, worker)
                == FAILURE)
                return FAILURE;
        }
    }
    return SUCCESS;
}

local int lya1d_tree_write_2pcf(struct cmdline_data *cmd,
                                struct global_data *gd,
                                const long double *num,
                                const long double *den,
                                const uint64_t *pairs,
                                uint64_t pair_count)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream;
    int bin;
    int write_failed;

    if (format_checked(path, sizeof(path), "radial tree 2PCF path",
                       "%s_lya1d%s", gd->fpfnamehistXi2pcfFileName,
                       EXTFILES) != 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree 2PCF output path is too long");
        return FAILURE;
    }
    stream = fopen(path, "w");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "cannot open radial tree 2PCF output '%s': %s",
                 path, strerror(errno));
        return FAILURE;
    }

    fprintf(stream, "# Exact radial Lyman-alpha forest 2PCF, interval tree\n");
    fprintf(stream, "# distinct-forest unordered pairs: %" PRIu64 "\n",
            pair_count);
    fprintf(stream, "# transverse separation is ignored\n");
    fprintf(stream, "# columns: bin radial_separation xi numerator denominator\n");
    for (bin = 0; bin < cmd->lya2RpBins; bin++) {
        long double numerator = pairs[bin] == 0 ? 0.0L : num[bin];
        long double denominator = pairs[bin] == 0 ? 0.0L : den[bin];
        long double correlation = denominator == 0.0L
                                ? 0.0L : numerator / denominator;
        REAL separation = ((REAL)bin + 0.5) * cmd->lya2RpMax
                        / (REAL)cmd->lya2RpBins;

        fprintf(stream, "%d %.17g %.17g %.17g %.17g\n", bin, separation,
                (double)correlation, (double)numerator, (double)denominator);
    }
    write_failed = ferror(stream);
    if (fclose(stream) != 0) write_failed = TRUE;
    if (write_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "failed writing radial tree 2PCF output '%s'", path);
        return FAILURE;
    }
    return SUCCESS;
}

global int searchcalc_lya_forest_1d_tree_omp(struct cmdline_data *cmd,
                                             struct global_data *gd,
                                             bodyptr table, INTEGER nbody)
{
    bodyptr *radial_order = NULL;
    bodyptr *forest_order = NULL;
    lya1d_tree_node *nodes = NULL;
    lya1d_tree_task *tree_tasks = NULL;
    lya1d_forest_task *forest_tasks = NULL;
    long double *all_num = NULL;
    long double *all_den = NULL;
    long double *same_num = NULL;
    long double *same_den = NULL;
    uint64_t *all_pairs = NULL;
    uint64_t *same_pairs = NULL;
    uint64_t *cross_pairs = NULL;
    size_t count = (size_t)nbody;
    size_t active_count = 0;
    size_t node_capacity = 0;
    size_t nodes_used = 0;
    size_t root = 0;
    size_t tree_task_count = 0;
    size_t forest_task_count = 0;
    size_t tree_block_count;
    size_t forest_block_count;
    size_t forest_count = 0;
    size_t i;
    uint64_t pair_count = 0;
    uint64_t node_visits = 0;
    uint64_t bulk_pairs = 0;
    int allocation_failed = FALSE;
    int runtime_failed = FALSE;
    ErrorMsg worker_error = "";
    int status = FAILURE;
    double cpustart = CPUTIME;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial tree Ly-alpha correlations require 3D catalog storage");
    return FAILURE;
#endif
#ifndef OPENMPCODE
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "radial tree Ly-alpha addon requires OPENMPMACHINE=1");
    return FAILURE;
#endif

    if (nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree Ly-alpha search received an empty catalog");
        return FAILURE;
    }

    if (cballs_calloc_checked((void **)&radial_order, count,
                              sizeof(*radial_order),
                              "radial tree sorted index", cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&forest_order, count,
                                 sizeof(*forest_order),
                                 "radial tree forest index", cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;

    for (i = 0; i < count; i++) {
        bodyptr body = table + i;

        if (Update(body) == FALSE || Mask(body) != MASK_NODE_VALID) continue;
        radial_order[active_count] = body;
        forest_order[active_count] = body;
        active_count++;
    }

    if (cballs_calloc_checked((void **)&all_num, (size_t)cmd->lya2RpBins,
                              sizeof(*all_num), "radial tree numerator",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&all_den,
                                 (size_t)cmd->lya2RpBins, sizeof(*all_den),
                                 "radial tree denominator", cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&same_num,
                                 (size_t)cmd->lya2RpBins, sizeof(*same_num),
                                 "same-forest tree numerator",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&same_den,
                                 (size_t)cmd->lya2RpBins, sizeof(*same_den),
                                 "same-forest tree denominator",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&all_pairs,
                                 (size_t)cmd->lya2RpBins, sizeof(*all_pairs),
                                 "radial tree pair counts", cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&same_pairs,
                                 (size_t)cmd->lya2RpBins, sizeof(*same_pairs),
                                 "same-forest pair counts", cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE
        || cballs_calloc_checked((void **)&cross_pairs,
                                 (size_t)cmd->lya2RpBins, sizeof(*cross_pairs),
                                 "cross-forest pair counts", cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;

    if (active_count == 0) goto setup_ready;
    qsort(radial_order, active_count, sizeof(*radial_order),
          lya1d_tree_radial_order);
    qsort(forest_order, active_count, sizeof(*forest_order),
          lya1d_tree_forest_order);

    if (lya1d_tree_node_count(active_count, &node_capacity) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree node count overflows size_t");
        goto setup_done;
    }
    if (cballs_calloc_checked((void **)&nodes, node_capacity,
                              sizeof(*nodes), "radial interval tree",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    root = lya1d_tree_build(nodes, &nodes_used, radial_order,
                            0, active_count);
    if (nodes_used != node_capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree construction produced an inconsistent node count");
        goto setup_done;
    }

    if (lya1d_tree_collect_tasks(nodes, root, root, cmd->lya2RpMax,
                                 cmd->lya2RpBins, NULL, 0,
                                 &tree_task_count) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree task count overflows size_t");
        goto setup_done;
    }
    if (tree_task_count > 0
        && cballs_calloc_checked((void **)&tree_tasks, tree_task_count,
                                 sizeof(*tree_tasks), "radial tree tasks",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    i = 0;
    if (lya1d_tree_collect_tasks(nodes, root, root, cmd->lya2RpMax,
                                 cmd->lya2RpBins, tree_tasks,
                                 tree_task_count, &i) == FAILURE
        || i != tree_task_count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree task construction failed");
        goto setup_done;
    }

    for (i = 0; i < active_count;) {
        size_t forest_end = i + 1;
        size_t pivot;

        while (forest_end < active_count
               && LyaForestId(forest_order[forest_end])
                  == LyaForestId(forest_order[i]))
            forest_end++;
        forest_count++;
        if (forest_end - i > 1) {
            for (pivot = i; pivot + 1 < forest_end;
                 pivot += (size_t)LYA1D_TREE_FOREST_PIVOT_CHUNK) {
                if (forest_task_count == SIZE_MAX) {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "same-forest task count overflows size_t");
                    goto setup_done;
                }
                forest_task_count++;
            }
        }
        i = forest_end;
    }
    if (forest_task_count > 0
        && cballs_calloc_checked((void **)&forest_tasks, forest_task_count,
                                 sizeof(*forest_tasks),
                                 "same-forest correction tasks",
                                 cmd->error_message,
                                 sizeof(cmd->error_message)) == FAILURE)
        goto setup_done;
    {
        size_t forest_first;
        size_t task_index = 0;

        for (forest_first = 0; forest_first < active_count;) {
            size_t forest_end = forest_first + 1;
            size_t pivot;

            while (forest_end < active_count
                   && LyaForestId(forest_order[forest_end])
                      == LyaForestId(forest_order[forest_first]))
                forest_end++;
            for (pivot = forest_first; pivot + 1 < forest_end;
                 pivot += (size_t)LYA1D_TREE_FOREST_PIVOT_CHUNK) {
                size_t pivot_end = MIN(
                    pivot + (size_t)LYA1D_TREE_FOREST_PIVOT_CHUNK,
                    forest_end - 1);

                forest_tasks[task_index].forest_end = forest_end;
                forest_tasks[task_index].pivot_first = pivot;
                forest_tasks[task_index].pivot_end = pivot_end;
                task_index++;
            }
            forest_first = forest_end;
        }
        if (task_index != forest_task_count) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "same-forest task construction failed");
            goto setup_done;
        }
    }

setup_ready:
    tree_block_count = tree_task_count == 0 ? 0
        : 1 + (tree_task_count - 1)
              / (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE;
    forest_block_count = forest_task_count == 0 ? 0
        : 1 + (forest_task_count - 1)
              / (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE;

    status = SUCCESS;
setup_done:
    status = lya_parallel_consensus(cmd, status, "Ly-alpha interval-tree setup");
    if (status == FAILURE) goto cleanup;
    if (active_count == 0) goto postprocess;

    ThreadCount(cmd, gd, (INTEGER)active_count, 0);
    const size_t first_task = lya_parallel_first(cmd);
    const size_t task_stride = lya_parallel_stride(cmd);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n%s: exact 1D interval tree; pixels=%zu forests=%zu nodes=%zu "
        "tree_tasks=%zu correction_tasks=%zu\n",
        cmd->searchMethod, active_count, forest_count, nodes_used,
        tree_task_count, forest_task_count);

#pragma omp parallel shared(allocation_failed,runtime_failed,worker_error,all_num,all_den,all_pairs,same_num,same_den,same_pairs,node_visits,bulk_pairs)
    {
        lya1d_tree_worker worker;
        ErrorMsg local_error = "";
        int worker_ready = lya1d_tree_worker_init(
            &worker, (size_t)cmd->lya2RpBins, local_error) == SUCCESS;
        size_t block_index;

        if (!worker_ready) {
#pragma omp critical(lya1d_tree_failure)
            {
                if (!allocation_failed)
                    snprintf(worker_error, sizeof(worker_error), "%s",
                             local_error);
                allocation_failed = TRUE;
            }
        }

#pragma omp barrier
#pragma omp for schedule(dynamic,1) ordered
        for (block_index = first_task; block_index < tree_block_count;
             block_index += task_stride) {
            size_t task_first = block_index
                              * (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE;
            size_t task_end = MIN(
                task_first + (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE,
                tree_task_count);
            size_t task;
            uint64_t block_visits = 0;
            uint64_t block_bulk_pairs = 0;

            if (worker_ready && !allocation_failed && !worker.failed) {
                worker.node_visits = 0;
                worker.bulk_pairs = 0;
                for (task = task_first; task < task_end; task++) {
                    if (lya1d_tree_walk(
                            nodes, radial_order,
                            tree_tasks[task].first_node,
                            tree_tasks[task].second_node,
                            cmd->lya2RpMax, cmd->lya2RpBins,
                            &worker) == FAILURE)
                        break;
                }
                block_visits = worker.node_visits;
                block_bulk_pairs = worker.bulk_pairs;
                if (worker.failed) {
#pragma omp critical(lya1d_tree_runtime_failure)
                    {
                        if (!runtime_failed)
                            snprintf(worker_error, sizeof(worker_error),
                                     "radial tree pair counter overflow");
                        runtime_failed = TRUE;
                    }
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !allocation_failed && !worker.failed) {
                    lya1d_tree_worker_commit(&worker, all_num, all_den,
                                             all_pairs);
                    node_visits += block_visits;
                    bulk_pairs += block_bulk_pairs;
                }
            }
        }

#pragma omp for schedule(dynamic,1) ordered
        for (block_index = first_task; block_index < forest_block_count;
             block_index += task_stride) {
            size_t task_first = block_index
                              * (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE;
            size_t task_end = MIN(
                task_first + (size_t)LYA1D_TREE_REDUCTION_BLOCK_SIZE,
                forest_task_count);
            size_t task;

            if (worker_ready && !allocation_failed && !worker.failed) {
                for (task = task_first; task < task_end; task++) {
                    if (lya1d_tree_add_forest_task(
                            forest_tasks + task, forest_order,
                            cmd->lya2RpMax, cmd->lya2RpBins,
                            &worker) == FAILURE)
                        break;
                }
                if (worker.failed) {
#pragma omp critical(lya1d_tree_runtime_failure)
                    {
                        if (!runtime_failed)
                            snprintf(worker_error, sizeof(worker_error),
                                     "same-forest pair counter overflow");
                        runtime_failed = TRUE;
                    }
                }
            }

#pragma omp ordered
            {
                if (worker_ready && !allocation_failed && !worker.failed)
                    lya1d_tree_worker_commit(&worker, same_num, same_den,
                                             same_pairs);
            }
        }
        if (worker_ready) lya1d_tree_worker_free(&worker);
    }

postprocess:
    if (allocation_failed || runtime_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree worker failed: %.2000s",
                 worker_error[0] != '\0' ? worker_error : "unknown error");
    }

    status = lya_parallel_consensus(cmd,
        allocation_failed || runtime_failed ? FAILURE : SUCCESS,
        "Ly-alpha interval-tree workers");
    if (status == FAILURE) goto cleanup;
    const size_t bins = (size_t)cmd->lya2RpBins;
    uint64_t counters[2] = {node_visits, bulk_pairs};
    if (lya_parallel_reduce_long_doubles(cmd, all_num, bins) == FAILURE
        || lya_parallel_reduce_long_doubles(cmd, all_den, bins) == FAILURE
        || lya_parallel_reduce_long_doubles(cmd, same_num, bins) == FAILURE
        || lya_parallel_reduce_long_doubles(cmd, same_den, bins) == FAILURE
        || lya_parallel_reduce_uint64(cmd, all_pairs, bins) == FAILURE
        || lya_parallel_reduce_uint64(cmd, same_pairs, bins) == FAILURE
        || lya_parallel_reduce_uint64(cmd, counters, 2) == FAILURE) {
        status = FAILURE;
        goto cleanup;
    }
    if (!lya_parallel_publish(cmd)) {
        status = SUCCESS;
        goto publication;
    }
    node_visits = counters[0];
    bulk_pairs = counters[1];
    status = FAILURE;

    for (i = 0; i < (size_t)cmd->lya2RpBins; i++) {
        long double all_denominator = all_den[i];
        long double tolerance;

        if (same_pairs[i] > all_pairs[i]) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "same-forest pair count exceeds all-pair count in bin %zu",
                     i);
            goto publication;
        }
        cross_pairs[i] = all_pairs[i] - same_pairs[i];
        all_num[i] -= same_num[i];
        all_den[i] -= same_den[i];
        if (cross_pairs[i] == 0) {
            all_num[i] = 0.0L;
            all_den[i] = 0.0L;
        }
        tolerance = 64.0L * LDBL_EPSILON
                  * (fabsl(all_denominator) + fabsl(same_den[i]) + 1.0L);
        if (all_den[i] < 0.0L) {
            if (-all_den[i] <= tolerance) {
                all_den[i] = 0.0L;
            } else {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "radial tree produced a negative denominator in bin %zu",
                         i);
                goto publication;
            }
        }
        if (all_den[i] == 0.0L) all_num[i] = 0.0L;
        if (UINT64_MAX - pair_count < cross_pairs[i]) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "radial tree total pair count overflows uint64_t");
            goto publication;
        }
        pair_count += cross_pairs[i];
    }

#ifdef LONGINT
    if (pair_count > (uint64_t)LONG_MAX
        || gd->nbbcalc > LONG_MAX - (INTEGER)pair_count) {
#else
    if (pair_count > (uint64_t)INT_MAX
        || gd->nbbcalc > INT_MAX - (INTEGER)pair_count) {
#endif
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "radial tree pair count exceeds the build INTEGER range");
        goto publication;
    }
    gd->nbbcalc += (INTEGER)pair_count;
    if (!cballs_opt_no_out_hist(cmd)
        && lya1d_tree_write_2pcf(cmd, gd, all_num, all_den,
                                cross_pairs, pair_count) == FAILURE)
        goto publication;

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: pairs=%" PRIu64 " node_visits=%" PRIu64
        " bulk_pairs=%" PRIu64 " CPU=%g\n",
        cmd->searchMethod, pair_count, node_visits, bulk_pairs, gd->cpusearch);
    status = SUCCESS;

publication:
    status = lya_parallel_consensus(cmd, status, "Ly-alpha interval-tree output");
cleanup:
    gd->cpusearch = CPUTIME - cpustart;
    free(radial_order);
    free(forest_order);
    free(nodes);
    free(tree_tasks);
    free(forest_tasks);
    free(all_num);
    free(all_den);
    free(same_num);
    free(same_den);
    free(all_pairs);
    free(same_pairs);
    free(cross_pairs);
    return status;
}
