/* Deterministic spatial construction for the native spherical shear trees.
 * Included privately by treeload.c, after its legacy helpers are declared. */
#ifndef CBALLS_NATIVE_SHEAR_OCTREE_BUILD_H
#define CBALLS_NATIVE_SHEAR_OCTREE_BUILD_H

#ifdef OPENMPCODE
#include <omp.h>
#endif
#include <time.h>

#define SHEAR_OCTREE_BUILD_CUTOFF 8192

typedef struct {
    INTEGER cells;
    int failure;                 /* 1 allocation, 2 depth, 3 duplicate */
} shear_octree_build_result;

static double shear_octree_build_wall(void)
{
    struct timespec stamp;
    if (clock_gettime(CLOCK_MONOTONIC, &stamp) != 0) return 0.0;
    return (double)stamp.tv_sec + 1.e-9*(double)stamp.tv_nsec;
}

static bool shear_octree_build_method(const struct global_data *gd)
{
    bool eligible = FALSE;
#ifdef OCTREESHEARSPHEREOMP
    eligible |= gd->searchMethod_int == OCTREESHEARSPHEREMETHOD;
#endif
#ifdef OCTREESHEARSPHERE2BALLSOMP
    eligible |= gd->searchMethod_int == OCTREESHEARSPHERE2BALLSOMPMETHOD;
#endif
#ifdef OCTREESHEARSPHERE2BALLSMPI
    eligible |= gd->searchMethod_int == OCTREESHEARSPHERE2BALLSMPIMETHOD;
#endif
    return eligible;
}

static bool shear_octree_bulk_eligible(const struct cmdline_data *cmd,
        const struct global_data *gd, INTEGER count)
{
    /* Random duplicate-position perturbations retain their original serial
     * insertion and RNG order. Small catalogs avoid partition scratch. */
    return shear_octree_build_method(gd) && count >= 2*SHEAR_OCTREE_BUILD_CUTOFF
        && !cballs_opt_no_check_equal_positions(cmd);
}

/* Once a branch is small, insert directly into that private subtree. This
 * avoids partitioning/copying its pointer range at every remaining level. */
static void shear_octree_insert_branch(struct cmdline_data *cmd, cellptr root,
        bodyptr *points, INTEGER count, real size, int depth,
        shear_octree_build_result *result)
{
    for (INTEGER i = 0; i < count; i++) {
        bodyptr point = points[i];
        cellptr parent = root;
        real parent_size = size;
        int level = depth, octant = subindex(point, parent);
        while (Subp(parent)[octant] != NULL) {
            if (Type(Subp(parent)[octant]) != CELL) {
                compute_vector separation;
                real distance2;
                cellptr child = NULL;
                DOTPSUBV(distance2, separation, Pos(point), Pos(Subp(parent)[octant]));
                if (distance2 == 0.0) { result->failure = 3; return; }
                if (level + 1 >= MAXLEVEL) { result->failure = 2; return; }
                if (cballs_calloc_checked((void **)&child, 1, sizeof(*child),
                        "spherical octree cell", NULL, 0) == FAILURE) {
                    result->failure = 1; return;
                }
                Type(child) = CELL;
                if (cballs_opt_read_mask(cmd)) Mask(child) = FALSE;
                for (int axis = 0; axis < NDIM; axis++)
                    Pos(child)[axis] = Pos(parent)[axis]
                        + (Pos(point)[axis] < Pos(parent)[axis] ? -parent_size : parent_size)/4;
                Subp(child)[subindex((bodyptr)Subp(parent)[octant], child)] = Subp(parent)[octant];
                Subp(parent)[octant] = (nodeptr)child;
                result->cells++;
            }
            parent = (cellptr)Subp(parent)[octant];
            octant = subindex(point, parent);
            parent_size /= 2;
            level++;
        }
        Subp(parent)[octant] = (nodeptr)point;
    }
}

static void shear_octree_partition_build(struct cmdline_data *cmd, cellptr parent,
        bodyptr *points, bodyptr *scratch, INTEGER count, real size, int depth,
        bool parallel, shear_octree_build_result *result)
{
    INTEGER counts[NSUB] = {0}, offsets[NSUB+1] = {0}, cursor[NSUB];
    shear_octree_build_result children[NSUB] = {{0}};
    result->cells = 1;
    result->failure = 0;
    if (depth >= MAXLEVEL) { result->failure = 2; return; }
    if (count > 1) {
        compute_vector separation;
        real distance2;
        DOTPSUBV(distance2, separation, Pos(points[0]), Pos(points[1]));
        if (distance2 == 0.0) { result->failure = 3; return; }
    }
    Nb(parent) = count;
    if (count <= SHEAR_OCTREE_BUILD_CUTOFF) {
        shear_octree_insert_branch(cmd, parent, points, count, size, depth, result);
        return;
    }
    for (INTEGER i = 0; i < count; i++) counts[subindex(points[i], parent)]++;
    for (int octant = 0; octant < NSUB; octant++) {
        offsets[octant+1] = offsets[octant] + counts[octant];
        cursor[octant] = offsets[octant];
    }
    /* A stable partition preserves input order within every spatial cell. */
    for (INTEGER i = 0; i < count; i++)
        scratch[cursor[subindex(points[i], parent)]++] = points[i];
    memcpy(points, scratch, (size_t)count*sizeof(*points));
    for (int octant = 0; octant < NSUB; octant++) {
        const INTEGER first = offsets[octant], members = counts[octant];
        cellptr child = NULL;
        if (members == 0) continue;
        if (members == 1) { Subp(parent)[octant] = (nodeptr)points[first]; continue; }
        /* Checked allocators return to this worker; no longjmp crosses an
         * OpenMP task. Attach immediately so ordinary cleanup owns the cell. */
        if (cballs_calloc_checked((void **)&child, 1, sizeof(*child),
                "spherical octree cell", NULL, 0) == FAILURE) {
            children[octant].failure = 1;
            continue;
        }
        Type(child) = CELL;
        if (cballs_opt_read_mask(cmd)) Mask(child) = FALSE;
        for (int axis = 0; axis < NDIM; axis++)
            Pos(child)[axis] = Pos(parent)[axis]
                + (Pos(points[first])[axis] < Pos(parent)[axis] ? -size : size)/4;
        Subp(parent)[octant] = (nodeptr)child;
#ifdef OPENMPCODE
        if (parallel && members >= SHEAR_OCTREE_BUILD_CUTOFF) {
#pragma omp task shared(children) firstprivate(child, first, members, octant, size, depth, points, scratch, cmd)
            shear_octree_partition_build(cmd, child, points+first, scratch+first,
                members, size/2, depth+1, TRUE, &children[octant]);
        } else
#endif
        {
            shear_octree_partition_build(cmd, child, points+first, scratch+first,
                members, size/2, depth+1, FALSE, &children[octant]);
        }
    }
#ifdef OPENMPCODE
    if (parallel) {
#pragma omp taskwait
    }
#endif
    /* Single-writer counters and deterministic failure selection. */
    for (int octant = 0; octant < NSUB; octant++) {
        result->cells += children[octant].cells;
        if (!result->failure) result->failure = children[octant].failure;
    }
}

static int shear_octree_bulk_build(struct cmdline_data *cmd,
        struct global_data *gd, bodyptr table, INTEGER count, int ifile)
{
    bodyptr *points = NULL, *scratch = NULL;
    shear_octree_build_result result = {0};
    /* Check the combined temporary-storage requirement before allocating. */
    if (cballs_malloc_checked((void **)&points, (size_t)count,
            2*sizeof(*points), "spherical octree partition buffers",
            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) return FAILURE;
    scratch = points + count;
    for (INTEGER i = 0; i < count; i++) points[i] = nthBody(table, i);
#ifdef OPENMPCODE
    if (cmd->numthreads > 1 && !omp_in_parallel()) {
#pragma omp parallel num_threads(cmd->numthreads)
        {
#pragma omp single nowait
            shear_octree_partition_build(cmd, roottable[ifile], points, scratch,
                count, gd->rSizeTable[ifile], 0, TRUE, &result);
        }
    } else
#endif
    {
        shear_octree_partition_build(cmd, roottable[ifile], points, scratch,
            count, gd->rSizeTable[ifile], 0, FALSE, &result);
    }
    free(points);
    gd->ncellTable[ifile] += result.cells - 1;  /* Root already accounted for. */
    gd->bytes_tot_cells += (result.cells - 1)*sizeof(cell);
    if (result.failure) {
        const char *message = result.failure == 1 ? "memory allocation failed"
            : result.failure == 2 ? "tree level out of range"
            : "two bodies have same position; consider option 'no-check-two-bodies-eq-pos'";
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "spherical octree build: %s", message);
        return FAILURE;
    }
    return SUCCESS;
}
#endif
