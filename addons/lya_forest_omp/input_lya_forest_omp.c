/* Six-column interchange reader for Lyman-alpha forest pixels.
 *
 * Columns: x y z delta weight forest_id
 * Coordinates are observer-centered comoving Cartesian positions.  cTreeBalls
 * recenters Pos() for tree construction, so the original radial distance and
 * line of sight are retained in addon-owned body fields.
 */

#include "globaldefs.h"
#include "lya_forest_defs.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>

#define LYA_INPUT_LINE_MAX 4096

local int lya_ascii_parse_line(struct cmdline_data *cmd, char *line,
                               unsigned long line_number, REAL values[5],
                               INTEGER *forest_id, int *is_data)
{
    char *cursor = line;
    char *end = NULL;
    int i;
    intmax_t parsed_id;

    *is_data = FALSE;
    while (isspace((unsigned char)*cursor)) cursor++;
    if (*cursor == '\0' || *cursor == '#') return SUCCESS;

    for (i = 0; i < 5; i++) {
        double value;
        errno = 0;
        value = strtod(cursor, &end);
        if (end == cursor || errno == ERANGE || !isfinite(value)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "lya-ascii line %lu: column %d is not a finite number",
                     line_number, i + 1);
            return FAILURE;
        }
        values[i] = (REAL)value;
        cursor = end;
    }

    errno = 0;
    parsed_id = strtoimax(cursor, &end, 10);
    if (end == cursor || errno == ERANGE || (intmax_t)(INTEGER)parsed_id != parsed_id) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii line %lu: forest_id is not a valid INTEGER",
                 line_number);
        return FAILURE;
    }
    cursor = end;
    while (isspace((unsigned char)*cursor)) cursor++;
    if (*cursor != '\0' && *cursor != '#') {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii line %lu: expected exactly six columns",
                 line_number);
        return FAILURE;
    }

    if (values[4] < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii line %lu: weight must be non-negative",
                 line_number);
        return FAILURE;
    }

    *forest_id = (INTEGER)parsed_id;
    *is_data = TRUE;
    return SUCCESS;
}

local int lya_ascii_read_line(FILE *stream, char *line, size_t line_size,
                              unsigned long line_number,
                              struct cmdline_data *cmd, int *has_line)
{
    size_t length;

    if (fgets(line, (int)line_size, stream) == NULL) {
        if (ferror(stream)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "lya-ascii: input error after line %lu", line_number);
            return FAILURE;
        }
        *has_line = FALSE;
        return SUCCESS;
    }

    length = strlen(line);
    if (length > 0 && line[length - 1] != '\n' && !feof(stream)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii line %lu exceeds %d bytes",
                 line_number, LYA_INPUT_LINE_MAX - 2);
        return FAILURE;
    }
    *has_line = TRUE;
    return SUCCESS;
}

global int inputdata_lya_ascii(struct cmdline_data *cmd,
                               struct global_data *gd,
                               string filename, int ifile)
{
    FILE *stream = NULL;
    char line[LYA_INPUT_LINE_MAX];
    unsigned long line_number = 0;
    INTEGER count = 0;
    INTEGER index = 0;
    int has_line;
    int is_data;
    REAL values[5];
    INTEGER forest_id;
    bodyptr p;
    REAL xmin, xmax, ymin, ymax, zmin, zmax;
    int status = FAILURE;

#if NDIM != 3
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "lya-ascii requires a three-dimensional build");
    return FAILURE;
#endif

    if (ifile < 0 || ifile >= MAXITEMS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii: invalid catalog index %d", ifile);
        return FAILURE;
    }

    stream = fopen(filename, "r");
    if (stream == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii: cannot open '%s': %s", filename, strerror(errno));
        return FAILURE;
    }

    while (TRUE) {
        line_number++;
        if (lya_ascii_read_line(stream, line, sizeof(line), line_number,
                                cmd, &has_line) == FAILURE)
            goto cleanup;
        if (!has_line) break;
        if (lya_ascii_parse_line(cmd, line, line_number, values,
                                 &forest_id, &is_data) == FAILURE)
            goto cleanup;
        if (is_data) count++;
    }
    if (count < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii: '%s' contains no data rows", filename);
        goto cleanup;
    }

    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)count,
                              sizeof(body), "Lyman-alpha body table",
                              cmd->error_message,
                              sizeof(cmd->error_message)) == FAILURE)
        goto cleanup;
    gd->bodytable_allocated = TRUE;
    gd->nbodyTable[ifile] = count;
    cmd->nbody = count;

    rewind(stream);
    clearerr(stream);
    line_number = 0;
    while (TRUE) {
        REAL distance2;
        REAL distance;
        line_number++;
        if (lya_ascii_read_line(stream, line, sizeof(line), line_number,
                                cmd, &has_line) == FAILURE)
            goto cleanup;
        if (!has_line) break;
        if (lya_ascii_parse_line(cmd, line, line_number, values,
                                 &forest_id, &is_data) == FAILURE)
            goto cleanup;
        if (!is_data) continue;
        if (index >= count) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "lya-ascii: input changed while it was being read");
            goto cleanup;
        }

        p = bodytable[ifile] + index;
        Pos(p)[0] = values[0];
        Pos(p)[1] = values[1];
        Pos(p)[2] = values[2];
        DOTVP(distance2, Pos(p), Pos(p));
        distance = rsqrt(distance2);
        if (!isfinite(distance) || distance <= 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "lya-ascii line %lu: observer distance must be positive",
                     line_number);
            goto cleanup;
        }

        Type(p) = BODY;
        Update(p) = TRUE;
        Update2(p) = TRUE;
        Mask(p) = MASK_NODE_VALID;
        Selected(p) = FALSE;
        Mass(p) = 1.0;
        Kappa(p) = values[3];
        Weight(p) = values[4];
        Id(p) = index + 1;
        LyaForestId(p) = forest_id;
        LyaDistance(p) = distance;
        LyaLOS(p)[0] = Pos(p)[0] / distance;
        LyaLOS(p)[1] = Pos(p)[1] / distance;
        LyaLOS(p)[2] = Pos(p)[2] / distance;
        index++;
    }
    if (index != count) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii: counted %" INTEGER_FMT
                 " rows but read %" INTEGER_FMT, count, index);
        goto cleanup;
    }

    xmin = xmax = Pos(bodytable[ifile])[0];
    ymin = ymax = Pos(bodytable[ifile])[1];
    zmin = zmax = Pos(bodytable[ifile])[2];
    DO_BODY(p, bodytable[ifile], bodytable[ifile] + count) {
        xmin = MIN(xmin, Pos(p)[0]); xmax = MAX(xmax, Pos(p)[0]);
        ymin = MIN(ymin, Pos(p)[1]); ymax = MAX(ymax, Pos(p)[1]);
        zmin = MIN(zmin, Pos(p)[2]); zmax = MAX(zmax, Pos(p)[2]);
    }
    gd->Box[0] = xmax - xmin;
    gd->Box[1] = ymax - ymin;
    gd->Box[2] = zmax - zmin;
    gd->input_comment = "Lyman-alpha x y z delta weight forest_id catalog";
    gd->bytes_tot += count * sizeof(body);
    status = SUCCESS;

cleanup:
    if (stream != NULL && fclose(stream) != 0 && status == SUCCESS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "lya-ascii: failed closing '%s': %s", filename,
                 strerror(errno));
        status = FAILURE;
    }
    if (status == FAILURE && bodytable[ifile] != NULL) {
        free(bodytable[ifile]);
        bodytable[ifile] = NULL;
        gd->nbodyTable[ifile] = 0;
        cmd->nbody = 0;
        gd->bodytable_allocated = FALSE;
    }
    return status;
}

#undef LYA_INPUT_LINE_MAX
