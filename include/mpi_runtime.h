#ifndef CBALLS_MPI_RUNTIME_H
#define CBALLS_MPI_RUNTIME_H

/* Context-local views of MPI_COMM_WORLD. MPI itself remains process-owned. */
typedef struct { int active, rank, size; } cballs_mpi_engine_state;
enum {
    CBALLS_MPI_FCFC_BALLTREE_2BALLS_MPI,
    CBALLS_MPI_FCFC_BALLTREE_2BALLS_MPI_3PCF,
    CBALLS_MPI_FCFC_BALLTREE_MPI,
    CBALLS_MPI_FCFC_BALLTREE_SHEAR_SPHERE_2BALLS_MPI,
    CBALLS_MPI_FCFC_KDTREE_2BALLS_MPI,
    CBALLS_MPI_FCFC_KDTREE_MPI,
    CBALLS_MPI_FCFC_KDTREE_SHEAR_SPHERE_2BALLS_MPI,
    CBALLS_MPI_LYA_FOREST_MPI,
    CBALLS_MPI_FCFC_OCTREE_2BALLS_MPI,
    CBALLS_MPI_CB3D_MPI,
    CBALLS_MPI_FCFC_OCTREE_BALLS4_MPI,
    CBALLS_MPI_FCFC_OCTREE_GGG_MPI,
    CBALLS_MPI_FCFC_OCTREE_SHEAR_SPHERE_2BALLS_MPI,
    CBALLS_MPI_STARTUP,
    CBALLS_MPI_ENGINE_COUNT
};
struct cmdline_data;
struct global_data;
int cballs_mpi_agree_catalogs(struct cmdline_data *, struct global_data *);
int cballs_mpi_bootstrap(struct cmdline_data *, struct global_data *, const char *);
int cballs_mpi_context_consensus(struct cmdline_data *, int, const char *);
int cballs_mpi_agree_text(struct cmdline_data *, const char *, const char *);
void cballs_mpi_local_begin(void);
void cballs_mpi_local_end(void);
void cballs_mpi_test_boundary(const char *);
const char *cballs_mpi_boundary_trace(void);
#ifdef CBALLS_MPI_ENABLED
int cballs_mpi_shared_prepare(cballs_mpi_engine_state *, struct cmdline_data *, struct global_data *, int);
int cballs_mpi_shared_finalize(cballs_mpi_engine_state *, struct cmdline_data *);
int cballs_mpi_shared_error(struct cmdline_data *, const char *, int);
int cballs_mpi_shared_consensus(cballs_mpi_engine_state *, struct cmdline_data *, int, int, const char *);
#endif
#endif
