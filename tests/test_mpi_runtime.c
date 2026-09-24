/* Real two-rank tests: borrowed/owned MPI, rank-local errors, and shutdown. */
#define global
#include "globaldefs.h"
#ifdef CBALLS_MPI_ENABLED
#include <mpi.h>
#include <assert.h>
int main(int argc, char **argv)
{
    struct cmdline_data cmd = {0};
    struct global_data gd = {0};
    cballs_mpi_engine_state a = {0}, b = {0};
    int borrowed = argc > 1 && strcmp(argv[1], "borrowed") == 0;
    int provided, finalized, rank, size;
    cmd.searchMethod = "mpi-runtime-test";
    if (borrowed) assert(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) == MPI_SUCCESS);
    assert(cballs_mpi_shared_prepare(&a, &cmd, &gd, TRUE) == SUCCESS);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size == 2 && a.size == size && a.rank == rank);
    assert(cballs_mpi_shared_prepare(&b, &cmd, &gd, TRUE) == SUCCESS);
    if (rank == 1) strcpy(cmd.error_message, "rank-local injected allocation failure");
    assert(cballs_mpi_shared_consensus(&a, &cmd, TRUE, rank == 1 ? FAILURE : SUCCESS,
                                      "allocation") == FAILURE);
    assert(strcmp(cmd.error_message, "rank-local injected allocation failure") == 0);
    cmd.error_message[0] = '\0';
    assert(cballs_mpi_shared_consensus(&b, &cmd, TRUE, SUCCESS, "recovery") == SUCCESS);
    assert(cballs_mpi_shared_finalize(&a, &cmd) == SUCCESS);
    MPI_Finalized(&finalized);
    assert(finalized == !borrowed);
    assert(cballs_mpi_shared_finalize(&b, &cmd) == SUCCESS);
    if (borrowed) {
        /* Cleanup cannot invalidate a communicator supplied by the host. */
        assert(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
        assert(cballs_mpi_shared_prepare(&a, &cmd, &gd, TRUE) == SUCCESS);
        MPI_Finalize();
    }
    assert(cballs_mpi_shared_prepare(&a, &cmd, &gd, TRUE) == FAILURE);
    assert(strstr(cmd.error_message, "after MPI_Finalize") != NULL);
    if (rank == 0) puts("PASS: shared MPI ownership, failure consensus, recovery, finalized rejection");
    return 0;
}
#else
int main(void) { puts("SKIP: MPI is not compiled"); return 77; }
#endif
