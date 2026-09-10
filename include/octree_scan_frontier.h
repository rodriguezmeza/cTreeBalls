#ifndef CTREEBALLS_OCTREE_SCAN_FRONTIER_H
#define CTREEBALLS_OCTREE_SCAN_FRONTIER_H

#include <limits.h>
#include <stdint.h>

/*
 * Expand the disjoint BALLS4 scan-level cover into a spatially ordered list
 * of body pivots.  Exact-body estimators can use this order for locality and
 * load balancing without replacing bodies by aggregate cell pivots.
 */
static inline int cballs_octree_scan_collect_bodies(
        struct cmdline_data *cmd, nodeptr node,
        bodyptr first, bodyptr finish, bool read_mask,
        bodyptr *pivots, size_t capacity, size_t *count)
{
    nodeptr child;

    if (node == NULL || (read_mask && Mask(node) == MASK_NODE_MASKED))
        return SUCCESS;
    if (Type(node) == CELL) {
        for (child = More(node); child != Next(node); child = Next(child))
            if (cballs_octree_scan_collect_bodies(
                    cmd, child, first, finish, read_mask,
                    pivots, capacity, count) == FAILURE)
                return FAILURE;
        return SUCCESS;
    }
    if (Type(node) != BODY && Type(node) != BODY3)
        return SUCCESS;
    if ((bodyptr)node < first || (bodyptr)node >= finish)
        return SUCCESS;
    if (read_mask && Mask(node) != MASK_NODE_VALID)
        return SUCCESS;
    if (*count >= capacity) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: scan-level body frontier overflow",
                 cmd->searchMethod);
        return FAILURE;
    }
    pivots[(*count)++] = (bodyptr)node;
    return SUCCESS;
}

static inline int cballs_octree_scan_body_order(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr first, bodyptr finish, int catalog,
        bodyptr **result, INTEGER *result_count)
{
    const bool read_mask = cballs_opt_read_mask(cmd);
    const ptrdiff_t difference = finish - first;
    size_t capacity;
    size_t count = 0;
    size_t expected = 0;
    bodyptr *pivots;

    if (result == NULL || result_count == NULL || first == NULL
        || finish < first || catalog < 0 || catalog >= gd->ninfiles
        || nodetablescanlevB4[catalog] == NULL
        || gd->nnodescanlevTableB4[catalog] <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: scan-level body frontier is unavailable",
                 cmd->searchMethod);
        return FAILURE;
    }
    *result = NULL;
    *result_count = 0;
    if (difference < 0 || (uintmax_t)difference > SIZE_MAX) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: scan-level body frontier size is invalid",
                 cmd->searchMethod);
        return FAILURE;
    }
    capacity = (size_t)difference;
    if (capacity > SIZE_MAX / sizeof(*pivots)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: scan-level body frontier allocation overflow",
                 cmd->searchMethod);
        return FAILURE;
    }
    pivots = capacity ? malloc(capacity * sizeof(*pivots)) : NULL;
    if (capacity && pivots == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: cannot allocate scan-level body frontier",
                 cmd->searchMethod);
        return FAILURE;
    }
    for (bodyptr pivot = first; pivot < finish; pivot++)
        if (!read_mask || Mask(pivot) == MASK_NODE_VALID)
            expected++;
    for (INTEGER task = 0;
         task < gd->nnodescanlevTableB4[catalog]; task++) {
        if (cballs_octree_scan_collect_bodies(
                cmd, nodetablescanlevB4[catalog][task], first, finish,
                read_mask, pivots, capacity, &count) == FAILURE) {
            free(pivots);
            return FAILURE;
        }
    }
    if (count != expected
#ifdef LONGINT
        || count > (size_t)LONG_MAX
#else
        || count > (size_t)INT_MAX
#endif
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: scan-level body frontier covers %zu pivots, expected %zu",
                 cmd->searchMethod, count, expected);
        free(pivots);
        return FAILURE;
    }
    *result = pivots;
    *result_count = (INTEGER)count;
    return SUCCESS;
}

#endif
