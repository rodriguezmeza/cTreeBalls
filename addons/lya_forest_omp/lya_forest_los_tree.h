#ifndef LYA_FOREST_LOS_TREE_H
#define LYA_FOREST_LOS_TREE_H

/* Internal neighbor backend. The owner keeps the catalog and octree alive. */
typedef struct lya_los_index lya_los_index;
typedef struct {
    size_t *seen;
    size_t *forests;
    size_t epoch;
    unsigned long long octree_nodes, forest_skips, forest_hits;
    unsigned long long radial_nodes, pixel_tests;
} lya_los_workspace;

typedef int (*lya_los_visit)(bodyptr pixel, REAL distance, void *context,
                             ErrorMsg error_message);

int lya_los_build(lya_los_index **result, bodyptr table, INTEGER count,
                  nodeptr root, ErrorMsg error_message);
void lya_los_free(lya_los_index *index);
size_t lya_los_forest_count(const lya_los_index *index);
int lya_los_workspace_init(const lya_los_index *index,
                            lya_los_workspace *workspace,
                            ErrorMsg error_message);
void lya_los_workspace_free(lya_los_workspace *workspace);
int lya_los_query(const lya_los_index *index, lya_los_workspace *workspace,
                  bodyptr pivot, REAL cutoff, lya_los_visit visit,
                  void *context, ErrorMsg error_message);

#endif
