#include <math.h>
#include <stdio.h>

#include "dual_node_sleef_log.h"

int main(void)
{
    enum { count = 4099 };
    double input[count];
    double output[count];
    float input_float[count];
    float output_float[count];
    double maximum_relative = 0.0;
    double maximum_float_relative = 0.0;

    for (int i = 0; i < count; i++) {
        input[i] = exp(-20.0 + 40.0 * (double)i / (double)(count - 1));
        input_float[i] = (float)input[i];
    }

    if (dual_node_sleef_double_lanes() <= 1
        || dual_node_sleef_float_lanes() <= 1) {
        fprintf(stderr, "SLEEF test was not compiled for a supported SIMD ISA\n");
        return 1;
    }
    dual_node_sleef_log_double(output, input, count);
    dual_node_sleef_log_float(output_float, input_float, count);

    for (int i = 0; i < count; i++) {
        const double reference = log(input[i]);
        const double relative = fabs(output[i] - reference)
            / fmax(1.0, fabs(reference));
        const double reference_float = log((double)input_float[i]);
        const double float_relative = fabs(
            (double)output_float[i] - reference_float)
            / fmax(1.0, fabs(reference_float));

        if (!isfinite(output[i]) || !isfinite(output_float[i])) {
            fprintf(stderr, "SLEEF vector log returned a non-finite value\n");
            return 1;
        }
        if (relative > maximum_relative) maximum_relative = relative;
        if (float_relative > maximum_float_relative)
            maximum_float_relative = float_relative;
    }
    if (!(maximum_relative <= 4.0e-15)
        || !(maximum_float_relative <= 3.0e-6)) {
        fprintf(stderr, "SLEEF vector log relative errors %.17g %.17g\n",
                maximum_relative, maximum_float_relative);
        return 1;
    }
    printf("PASS: SLEEF SIMD log lanes=%d/%d relative errors %.3g/%.3g\n",
           dual_node_sleef_double_lanes(), dual_node_sleef_float_lanes(),
           maximum_relative, maximum_float_relative);
    return 0;
}
