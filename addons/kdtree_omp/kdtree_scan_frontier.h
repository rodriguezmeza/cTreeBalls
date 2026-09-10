#ifndef CBALLS_KDTREE_SCAN_FRONTIER_H
#define CBALLS_KDTREE_SCAN_FRONTIER_H

#ifndef KDTREE_SCAN_FRONTIER_TARGET
#define KDTREE_SCAN_FRONTIER_TARGET ((INTEGER)256)
#endif

typedef struct {
    INTEGER count;
    bool tree_order;
} kdtree_scan_frontier;

static inline kdtree_scan_frontier kdtree_scan_frontier_make(
        const ballxptr tree, INTEGER pivot_count, bool complete_auto_catalog)
{
    kdtree_scan_frontier frontier = {1, FALSE};

#ifdef BALLS4SCANLEV
    frontier.tree_order = complete_auto_catalog && tree != NULL
        && tree->ntab != NULL && tree->nsplit > 0;
    while (frontier.count < KDTREE_SCAN_FRONTIER_TARGET) {
        const INTEGER limit = frontier.tree_order
            ? (INTEGER)tree->nsplit : pivot_count;

        if (frontier.count > limit / 2) break;
        frontier.count *= 2;
    }
#else
    (void)tree;
    (void)pivot_count;
    (void)complete_auto_catalog;
#endif
    return frontier;
}

static inline void kdtree_scan_frontier_range(
        const kdtree_scan_frontier *frontier, const ballxptr tree,
        INTEGER pivot_count, INTEGER pivot_offset, INTEGER task,
        INTEGER *first, INTEGER *end)
{
    if (frontier->tree_order) {
        const ballnode *node = &tree->ntab[frontier->count + task];

        *first = (INTEGER)node->first;
        *end = (INTEGER)node->last + 1;
        return;
    }

    {
        const INTEGER quotient = pivot_count / frontier->count;
        const INTEGER remainder = pivot_count % frontier->count;
        const INTEGER first_extra = MIN(task, remainder);
        const INTEGER end_task = task + 1;
        const INTEGER end_extra = MIN(end_task, remainder);

        *first = pivot_offset + task * quotient + first_extra;
        *end = pivot_offset + end_task * quotient + end_extra;
    }
}

static inline bodyptr kdtree_scan_frontier_body(
        const kdtree_scan_frontier *frontier, const ballxptr tree,
        bodyptr body_table, INTEGER index)
{
    return frontier->tree_order ? tree->bptr[index] : body_table + index;
}

#endif /* CBALLS_KDTREE_SCAN_FRONTIER_H */
