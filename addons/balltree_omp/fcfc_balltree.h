/*
 * FCFC-style ball tree for cTreeBalls.
 *
 * The construction strategy is adapted from FCFC by Cheng Zhao:
 * https://github.com/cheng-zhao/FCFC
 * Copyright (c) 2020--2022 Cheng Zhao, used under the MIT license.
 */

#ifndef _fcfc_balltree_h
#define _fcfc_balltree_h

typedef struct {
    vector center;                 /* centre of the enclosing sphere */
    cballs_storage_real radius;   /* conservative enclosing radius */
    vector cmpos;                  /* mass-weighted aggregate position */
    cballs_storage_real aggregate_radius; /* opening radius around cmpos */
    real kappa;                    /* mean scalar field */
    real weight;                   /* total point weight */
    real kappa_sum;                /* unweighted scalar-field sum */
    real kappa_sq_sum;             /* sum of squared scalar fields */
    real field_weight_sum;         /* sum of Weight(point) */
    real field_weight_sq_sum;      /* sum of squared point weights */
    real weighted_kappa_sum;       /* sum of Weight(point) * Kappa(point) */
    real weighted_kappa_sq_sum;    /* sum of (Weight(point) * Kappa(point))^2 */
    INTEGER first;
    INTEGER last;
    INTEGER left;
    INTEGER right;
} fcfc_ballnode;

#ifdef SINGLEP
typedef struct {
    vector pos;
    real kappa;
    real weighted_kappa;
} fcfc_ballpoint;
#endif

typedef struct {
    INTEGER npoint;
    INTEGER nnode;
    INTEGER capacity;
    int max_depth;
    bodyptr *bptr;
    fcfc_ballnode *nodes;
#ifdef SINGLEP
    fcfc_ballpoint *packed_points;
#endif
} fcfc_balltree, *fcfc_balltreeptr;

int fcfc_balltree_build(struct cmdline_data *, struct global_data *,
                        bodyptr, INTEGER, int, fcfc_balltreeptr *);
int fcfc_balltree_frontier(struct cmdline_data *, fcfc_balltreeptr,
                           INTEGER, INTEGER **, INTEGER *);
void fcfc_balltree_free(fcfc_balltreeptr);

#endif /* !_fcfc_balltree_h */
