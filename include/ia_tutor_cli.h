#ifndef CBALLS_IA_TUTOR_CLI_H
#define CBALLS_IA_TUTOR_CLI_H

/* Standalone information mode: no catalog allocation, MPI startup, shell,
 * Python embedding or additional native link dependency is required. */
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

static int cballs_tutor_token(const char *value)
{
    while (*value) {
        while (*value == ',' || *value == ' ' || *value == '\t') value++;
        size_t size = strcspn(value, ", \t");
        if (size == 8 && (!strncmp(value, "IA-tutor", 8) || !strncmp(value, "ia-tutor", 8)))
            return 1;
        value += size;
    }
    return 0;
}

static int cballs_tutor_resource(const char *root, char *script, size_t capacity)
{
    int n = snprintf(script, capacity, "%s/scripts/ia_tutor.py", root);
    return n >= 0 && (size_t)n < capacity && access(script, R_OK) == 0;
}

static void cballs_tutor_executable(const char *argument, char *path)
{
    path[0] = '\0';
#ifdef __APPLE__
    char temporary[PATH_MAX];
    uint32_t size = sizeof(temporary);
    if (_NSGetExecutablePath(temporary, &size) == 0 && realpath(temporary, path)) return;
#elif defined(__linux__)
    ssize_t size = readlink("/proc/self/exe", path, PATH_MAX-1);
    if (size > 0) { path[size] = '\0'; return; }
#endif
    if (strchr(argument, '/') && realpath(argument, path)) return;
    const char *search = getenv("PATH");
    if (search) {
        char candidate[PATH_MAX];
        while (*search) {
            size_t length = strcspn(search, ":");
            if (length < PATH_MAX/2) {
                int n = snprintf(candidate, sizeof(candidate), "%.*s/%s",
                                 (int)(length ? length : 1), length ? search : ".", argument);
                if (n > 0 && (size_t)n < sizeof(candidate) && access(candidate, X_OK) == 0
                    && realpath(candidate, path)) return;
            }
            search += length;
            if (*search == ':') search++;
        }
    }
    snprintf(path, PATH_MAX, "%s", argument);
}

/* Return 1 when this invocation belongs to the tutor, with its exit status.
 * A successful exec never returns. Keep this before the scientific parser. */
static int cballs_ia_tutor_cli(int argc, char **argv, int *status)
{
    int found = 0;
    for (int i=1; i<argc; i++)
        if (!strncmp(argv[i], "options=", 8) && cballs_tutor_token(argv[i]+8)) {
            found = i;
            if (strcmp(argv[i]+8, "IA-tutor") && strcmp(argv[i]+8, "ia-tutor")) {
                fprintf(stderr, "IA-tutor is a separate information mode; use cballs options=IA-tutor alone, followed by tutor arguments.\n");
                *status = 2; return 1;
            }
            break;
        }
    if (!found) return 0;
    const char *sizes[] = {"OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "MV2_COMM_WORLD_SIZE", NULL};
    for (int i=0; sizes[i]; i++) {
        const char *value = getenv(sizes[i]);
        if (value && strtol(value, NULL, 10) > 1) {
            fprintf(stderr, "IA-tutor is a serial conversation; launch it outside mpiexec/srun.\n");
            *status = 2; return 1;
        }
    }
    if (getenv("SLURM_PROCID") && getenv("SLURM_NTASKS")
        && strtol(getenv("SLURM_NTASKS"), NULL, 10) > 1) {
        fprintf(stderr, "IA-tutor needs one terminal, not multiple Slurm tasks.\n");
        *status = 2; return 1;
    }
    char executable[PATH_MAX], root[PATH_MAX], script[PATH_MAX];
    cballs_tutor_executable(argv[0], executable);
    const char *override = getenv("CBALLS_TUTOR_ROOT");
    int located = 0;
    if (override && *override) {
        if (strlen(override) < sizeof(root)) {
            strcpy(root, override);
            located = cballs_tutor_resource(root, script, sizeof(script));
        }
    } else {
        snprintf(root, sizeof(root), "%s", executable);
        for (int level=0; level<5; level++) {
            char *slash = strrchr(root, '/');
            if (!slash) break;
            *slash = '\0';
            if (cballs_tutor_resource(root, script, sizeof(script))) { located = 1; break; }
        }
#ifdef __CTREEBALLSDIR__
        if (!located && strlen(__CTREEBALLSDIR__) < sizeof(root)) {
            strcpy(root, __CTREEBALLSDIR__);
            located = cballs_tutor_resource(root, script, sizeof(script));
        }
#endif
    }
    if (!located) {
        fprintf(stderr, "Cannot find scripts/ia_tutor.py. Set CBALLS_TUTOR_ROOT to the cTreeBalls source directory.\n");
        *status = 2; return 1;
    }
    const char *python = getenv("CBALLS_TUTOR_PYTHON");
    if (!python || !*python) python = "python3";
    char **arguments = calloc((size_t)argc+9, sizeof(*arguments));
    if (!arguments) { *status = 2; return 1; }
    int count = 0;
    arguments[count++] = (char *)python;
    arguments[count++] = "-B";
    arguments[count++] = script;
    arguments[count++] = "--root";
    arguments[count++] = root;
    arguments[count++] = "--cballs";
    arguments[count++] = executable;
    for (int i=1; i<argc; i++) {
        if (i == found || (i == found+1 && !strcmp(argv[i], "--"))) continue;
        arguments[count++] = argv[i];
    }
    arguments[count] = NULL;
    fflush(NULL);
    execvp(python, arguments);
    fprintf(stderr, "Cannot start tutor Python '%s': %s. Set CBALLS_TUTOR_PYTHON to a Python 3.9+ executable (one path, no shell options).\n", python, strerror(errno));
    free(arguments);
    *status = 127;
    return 1;
}
#endif
