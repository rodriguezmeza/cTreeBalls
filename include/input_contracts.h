#ifndef CBALLS_INPUT_CONTRACTS_H
#define CBALLS_INPUT_CONTRACTS_H

/* Format-independent checks, after conversion and before geometry/overrides.
 * The caller retains ownership and supplies its normal cleanup path. */
static inline int cballs_input_finite(struct cmdline_data *cmd,
                                     const char *filename, size_t row,
                                     const char *field, double value)
{
    if (isfinite(value)) return SUCCESS;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "catalog input '%s': data row %zu, %s: non-finite value",
             filename, row, field);
    return FAILURE;
}

static inline int cballs_input_finite_array(struct cmdline_data *cmd,
                                             const char *filename,
                                             const void *data, size_t count,
                                             size_t item_size,
                                             const char *field)
{
    for (size_t i = 0; i < count; ++i) {
        double value = item_size == sizeof(float)
            ? ((const float *)data)[i] : ((const double *)data)[i];
        if (cballs_input_finite(cmd, filename, i + 1, field, value) == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

static inline int cballs_input_validate_bodies(struct cmdline_data *cmd,
                                                const char *filename,
                                                bodyptr catalog, INTEGER count)
{
    if (catalog == NULL || count < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "catalog input '%s': no bodies selected", filename);
        return FAILURE;
    }
    for (INTEGER i = 0; i < count; ++i) {
        bodyptr p = catalog + i;
        if (Mask(p) != MASK_NODE_MASKED && Mask(p) != MASK_NODE_VALID) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "catalog input '%s': data row %zu, mask must be 0 or 1",
                     filename, (size_t)i + 1);
            return FAILURE;
        }
        for (int k = 0; k < NDIM; ++k)
            if (cballs_input_finite(cmd, filename, (size_t)i + 1,
                                    "position", Pos(p)[k]) == FAILURE)
                return FAILURE;
        if (cballs_input_finite(cmd, filename, (size_t)i + 1,
                                "kappa", Kappa(p)) == FAILURE
            || cballs_input_finite(cmd, filename, (size_t)i + 1,
                                   "weight", Weight(p)) == FAILURE)
            return FAILURE;
#ifdef THREEPCFSHEAR
        if ((scanopt(cmd->options, "pos-and-shear")
             || (cmd->searchMethod && strstr(cmd->searchMethod, "shear")))
            && (cballs_input_finite(cmd, filename, (size_t)i + 1,
                                    "gamma1", Gamma1(p)) == FAILURE
                || cballs_input_finite(cmd, filename, (size_t)i + 1,
                                       "gamma2", Gamma2(p)) == FAILURE))
            return FAILURE;
#endif
    }
    return SUCCESS;
}
#endif
