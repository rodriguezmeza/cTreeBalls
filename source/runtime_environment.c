/* Runtime environment help. Keep the allowlist alongside the native readers;
 * never dump the full process environment (which may contain credentials). */
#include "globaldefs.h"
#include "environment_info.h"
#include "resource_checks.h"
#include "cballs_make_info.h"

typedef struct {
    const char *name;
    const char *fallback;
    const char *description;
} cballs_environment_entry;

/* Render values as data, including empty strings, without terminal escapes or
 * newlines supplied by the environment. Bound output for accidentally huge values. */
static void environment_value(const char *value)
{
    size_t i;
    if (value == NULL) {
        fputs("<unset>", stdout);
        return;
    }
    putchar('"');
    for (i = 0; value[i] && i < 1024; ++i) {
        unsigned char c = (unsigned char)value[i];
        switch (c) {
        case '\\': fputs("\\\\", stdout); break;
        case '"': fputs("\\\"", stdout); break;
        case '\n': fputs("\\n", stdout); break;
        case '\r': fputs("\\r", stdout); break;
        case '\t': fputs("\\t", stdout); break;
        default:
            if (c < 32 || c >= 127) printf("\\x%02x", (unsigned)c);
            else putchar(c);
        }
    }
    putchar('"');
    if (value[i]) fputs(" [truncated after 1024 bytes]", stdout);
}

static void environment_group(const char *title,
        const cballs_environment_entry *entries, size_t count)
{
    printf("\n%s\n", title);
    for (size_t i = 0; i < count; ++i) {
        printf("  %s = ", entries[i].name);
        environment_value(getenv(entries[i].name));
        printf("\n    When unset: %s. %s\n", entries[i].fallback,
               entries[i].description);
    }
}

static void environment_unset(const cballs_environment_entry *entries, size_t count)
{
    fputs("  unset", stdout);
    for (size_t i = 0; i < count; ++i) printf(" %s", entries[i].name);
    putchar('\n');
}

#define ENTRY_COUNT(entries) (sizeof(entries) / sizeof((entries)[0]))

void cballs_print_environment(void)
{
    static const cballs_environment_entry native[] = {
        {"CBALLS_MEMORY_BUDGET_MB", CBALLS_STRINGIFY(CBALLS_DEFAULT_MEMORY_BUDGET_MIB) " MiB",
         "Positive decimal integer, 1 MiB = 1048576 bytes. Checked allocation/histogram plans per MPI rank; not a process RSS cap. Empty, zero and malformed values are invalid."},
        {"CBALLS_SHEAR_PROFILE", "disabled",
         "Extra shear tree/search timers and counters; ANY defined value enables profiling, including 0 and the empty string. Unset for timing without this instrumentation."},
        {"CBALLS_SHEAR_PIVOT_TOL", CBALLS_MAKE_CBALLS_SHEAR_PIVOT_TOL " radians (compiled)",
         "Finite [0,3] phase budget for eligible shear-pivot-reuse searches; 0 disables reuse. Requires the options token shear-pivot-reuse and a supporting build/engine; not a fractional correlation-error bound."},
        {"CBALLS_SHEAR_BIN_THETA", CBALLS_MAKE_CBALLS_SHEAR_BIN_THETA " bin widths (compiled)",
         "Finite [0,1] radial-bin tolerance for eligible shear-pivot-reuse searches; 0 uses strict bin containment. Independent of theta."},
        {"CBALLS_SCALAR_PIVOT_TOL", CBALLS_MAKE_CBALLS_SCALAR_PIVOT_TOL " radians (compiled)",
         "Finite [0,3] phase budget for eligible scalar-pivot-reuse 3PCF searches; 0 disables reuse. Requires scalar-pivot-reuse, positive theta, nonperiodic data and no smooth-pivot. Not a fractional correlation-error bound; 3 does not imply a 3-percent error."},
        {"CBALLS_SCALAR_BIN_THETA", CBALLS_MAKE_CBALLS_SCALAR_BIN_THETA " bin widths (compiled)",
         "Finite [0,1] radial-bin tolerance for eligible scalar-pivot-reuse searches; 0 uses strict bin containment. Independent of theta and the dual-node-bin-theta options token; 0.5 is not a 50-percent error bound."},
        {"CBALLS_SHEAR_TASKS", "0 (automatic, engine/thread/rank dependent)",
         "Integer [0,INT_MAX-8]. Requested global shear frontier task target; 0 selects automatic sizing. Clamped by available work and task histogram memory."},
        {"CBALLS_SHEAR_TASK_MEMORY_MB", "256 MiB per rank",
         "Integer [1,min(1048576,SIZE_MAX>>20)] MiB for shear task histograms. At least one histogram is required even if larger than this limit. Does not size Ly-alpha trees or cap total process memory."},
        {"CBALLS_TUTOR_PYTHON", "python3 (also when empty)",
         "IA-tutor only: Python executable name or path, without shell arguments."},
        {"CBALLS_TUTOR_ROOT", "automatic script/source-root discovery (also when empty)",
         "IA-tutor only: directory containing scripts/ia_tutor.py."}
    };
    static const cballs_environment_entry openmp[] = {
        {"OMP_NUM_THREADS", "OpenMP runtime default",
         "Initial OpenMP thread count; cTreeBalls explicitly sets numberThreads (Python --threads), which overrides this initial count."},
        {"OMP_DYNAMIC", "OpenMP runtime default",
         "TRUE/FALSE: permit runtime adjustment of team sizes; FALSE requests fixed sizing subject to runtime limits."},
        {"OMP_THREAD_LIMIT", "OpenMP runtime default",
         "Upper limit on OpenMP threads; can constrain numberThreads."},
        {"OMP_PROC_BIND", "OpenMP runtime default",
         "Thread binding policy (for example FALSE, TRUE, CLOSE or SPREAD); support depends on the runtime."},
        {"OMP_PLACES", "OpenMP runtime default",
         "Thread placement (for example cores or threads); interpreted by the runtime."},
        {"OMP_SCHEDULE", "OpenMP runtime default",
         "Schedule for schedule(runtime) loops only; does not override explicit static/dynamic schedules in the engines."},
        {"OMP_WAIT_POLICY", "OpenMP runtime default",
         "ACTIVE/PASSIVE waiting policy; trades CPU use against wake-up latency."},
        {"OMP_STACKSIZE", "OpenMP runtime default",
         "Worker thread stack size (runtime syntax); separate from heap memory and the cTreeBalls allocation budget."},
        {"OMP_MAX_ACTIVE_LEVELS", "OpenMP runtime default",
         "Maximum active nested parallel levels; does not add parallelism to a serial engine."},
        {"OMP_DISPLAY_ENV", "OpenMP runtime default",
         "Ask the runtime to display its own configuration at initialization (TRUE or VERBOSE where supported)."},
        {"GOMP_CPU_AFFINITY", "GNU OpenMP runtime default",
         "GNU libgomp CPU affinity list; ignored by other OpenMP runtimes."},
        {"GOMP_SPINCOUNT", "GNU OpenMP runtime default",
         "GNU libgomp spin-wait iterations before sleeping; affects CPU usage and latency."},
        {"KMP_AFFINITY", "Intel/LLVM OpenMP runtime default",
         "Vendor affinity policy; runtime-specific, not used by GNU libgomp."},
        {"KMP_BLOCKTIME", "Intel/LLVM OpenMP runtime default",
         "Vendor worker wait time; runtime-specific, not used by GNU libgomp."}
    };
    static const cballs_environment_entry libraries[] = {
        {"OPENBLAS_NUM_THREADS", "BLAS library default", "OpenBLAS thread limit in linked libraries or NumPy; separate from native search threads."},
        {"MKL_NUM_THREADS", "BLAS library default", "Intel MKL thread limit; separate from native search threads."},
        {"MKL_DYNAMIC", "MKL library default", "Permit dynamic MKL team sizing; affects MKL only."},
        {"MKL_THREADING_LAYER", "MKL library default", "Select MKL threading backend; affects MKL only."},
        {"VECLIB_MAXIMUM_THREADS", "Accelerate library default", "Apple Accelerate thread limit where supported; separate from native search threads."},
        {"BLIS_NUM_THREADS", "BLIS library default", "BLIS thread limit; separate from native search threads."},
        {"NUMEXPR_NUM_THREADS", "NumExpr library default", "Python NumExpr worker count; does not configure native cTreeBalls loops."},
        {"GSL_RNG_TYPE", "GSL generator default", "Generator selection for synthetic data when built with USEGSL; not used for a supplied forest catalog."},
        {"GSL_RNG_SEED", "GSL seed default", "Read by GSL setup when enabled; cTreeBalls subsequently overrides the seed with its seed parameter."}
    };
    static const cballs_environment_entry launch[] = {
        {"OMPI_COMM_WORLD_SIZE", "not in an Open MPI launch", "Launcher-managed MPI world size; its presence is checked during native MPI startup."},
        {"OMPI_COMM_WORLD_RANK", "not in an Open MPI launch", "Launcher-managed MPI rank."},
        {"PMI_SIZE", "not in a PMI launch", "Launcher-managed MPI world size; its presence is checked during native MPI startup."},
        {"PMI_RANK", "not in a PMI launch", "Launcher-managed PMI rank."},
        {"PMIX_RANK", "not in a PMIx launch", "Launcher-managed rank; its presence is checked during native MPI startup."},
        {"MV2_COMM_WORLD_SIZE", "not in an MVAPICH launch", "Launcher-managed world size, also used by helper scripts."},
        {"SLURM_NTASKS", "not set by Slurm", "Scheduler-managed task count; used by helper scripts."},
        {"SLURM_PROCID", "not set by Slurm", "Scheduler-managed process rank; used by helper scripts."}
    };
    static const cballs_environment_entry frontend[] = {
        {"CTREEBALLS_SOURCE_ROOT", "Python driver's source-root discovery",
         "Python driver/helper source location; does not change native estimator settings."}
    };
    size_t budget = cballs_memory_budget();

    printf("cTreeBalls runtime environment\n"
           "Usage: cballs options=print-environment (alias: print-env).\n"
           "Read-only report: quoted current values, <unset>, and defaults when unset.\n"
           "Values are environment strings; each consumer validates them when used.\n"
           "Compiled profile: OpenMP=%s, MPI=%s, Ly-alpha=%s, shear=%s, precision=%s.\n",
           CBALLS_MAKE_OPENMPMACHINE, CBALLS_MAKE_BUILD_MPI, CBALLS_MAKE_BUILD_LYA,
           CBALLS_MAKE_BUILD_SHEAR, CBALLS_MAKE_BUILD_PRECISION);
    environment_group("cTreeBalls controls (scope depends on the engine/build)", native, ENTRY_COUNT(native));
    if (budget)
        printf("\n  Parsed native memory budget: %zu MiB (%zu bytes) per rank.\n", budget/1048576, budget);
    else
        fputs("\n  INVALID CBALLS_MEMORY_BUDGET_MB: normal allocations will fail. Set a positive decimal integer or unset it.\n", stdout);
    fputs("  Raising a budget does not provide RAM. Native trees, Python/NumPy inputs and runtime overhead also consume memory.\n"
          "  Ly-alpha kernel/bin-slop controls are parameters, not CBALLS_SHEAR_* variables.\n"
          "  theta, nsmooth and rsmooth are search parameters; support is engine-specific. See print-search-methods.\n", stdout);

    fputs(
        "Accuracy: theta and pivot phase/bin tolerances are geometric controls, not 5% or 10% error guarantees. Compare every nonzero coefficient with abs(fast-exact)/abs(exact), using complex magnitudes; exact-zero coefficients must stay zero. L2 and symmetric differences are different metrics.\n"
        "Speed: measure repeated wall times on an otherwise idle host at matching input, bins, multipoles and measured accuracy. Profile scalar two-ball runs separately with options=no-out-Hist,dual-node-profile; remove profiling for timing trials. CPU time sums worker work.\n"
        "Memory: compact scalar octree/PCA/KD tree node capacity is sized automatically from the required tree shape; this is allocation policy, not smoothing, and does not guarantee lower peak RSS or faster traversal.\n"
        "Python kappa/shear HEALPix inputs: --fits-chunk-pixels bounds coordinate scratch (default 1048576); DES tables use --fits-chunk-rows. Kappa/shear/Ly-alpha drivers and native Python benchmarks accept --transfer-catalog-ownership (opt-in; default copies). Native particle storage and per-rank MPI copies remain.\n"
        "These are Python arguments, not native options or environment controls. See tests/python/README_memory_controls.md.\n", stdout);

    environment_group("OpenMP controls (external runtime; ignored in serial builds)", openmp, ENTRY_COUNT(openmp));
    environment_group("Numerical library / Python controls (only if that library is used)", libraries, ENTRY_COUNT(libraries));
    environment_group("Launch context (diagnostic only; do not clear inside MPI/Slurm jobs)", launch, ENTRY_COUNT(launch));
    environment_group("Python frontend helpers", frontend, ENTRY_COUNT(frontend));

    fputs("\nClearing definitions (bash/zsh)\n"
          "  unset CBALLS_MEMORY_BUDGET_MB\n"
          "    Restores the compiled memory default; does not mean unlimited memory.\n"
          "  unset CBALLS_SHEAR_PROFILE\n"
          "    Disables profiling; setting it to 0 or an empty string does NOT disable it.\n"
          "  Clear all cTreeBalls-owned controls listed above in the current shell:\n", stdout);
    environment_unset(native, ENTRY_COUNT(native));
    fputs("  Optional: also clear external OpenMP overrides (affects other programs started here):\n", stdout);
    environment_unset(openmp, ENTRY_COUNT(openmp));
    fputs("  Optional: also clear numerical-library/Python helper overrides:\n", stdout);
    environment_unset(libraries, ENTRY_COUNT(libraries));
    environment_unset(frontend, ENTRY_COUNT(frontend));
    fputs("  Clear only for one child command, keeping the current shell unchanged:\n"
          "    env -u CBALLS_MEMORY_BUDGET_MB -u CBALLS_SHEAR_PROFILE ./cballs options=print-environment\n"
          "  For a Python run, use the same env -u prefix before python3 and your script.\n"
          "  In Python: import os; os.environ.pop('CBALLS_MEMORY_BUDGET_MB', None)\n"
          "  Change native controls before Run(); do not mutate the environment during a run.\n"
          "  Set OpenMP/BLAS controls before starting Python (libraries may read them at import/initialization); restart an existing interpreter if needed.\n"
          "  NAME=\"\" keeps NAME defined; use unset NAME to remove it. Unsetting restores fallbacks; it does not remove options tokens such as shear-pivot-reuse.\n"
          "  For future shells/jobs, remove the corresponding export from your shell startup file, job script or environment activation configuration.\n"
          "  A child program cannot clear definitions in its parent shell; this report changes nothing.\n"
          "  MPI/scheduler launch variables are excluded from the reset commands.\n"
          "\nSee options=make-info for compiled defaults (including THETA and build-time engine settings),\n"
          "options=print-options for runtime option tokens, and options=print-search-methods for engine support.\n"
          "Append ,no-stop to print-environment to report settings and then continue a calculation.\n", stdout);
}

#undef ENTRY_COUNT
