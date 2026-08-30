#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "stdinc.h"
#include "common_defs.h"

#define CHECK(condition, message)                    \
    do {                                             \
        if (!(condition)) {                          \
            fprintf(stderr, "FAIL: %s\n", message); \
            return EXIT_FAILURE;                    \
        }                                            \
    } while (0)

static int parse_line(const char *line, char *name, size_t name_size,
                      char *value, size_t value_size, int *is_data,
                      char *errmsg, size_t errmsg_size)
{
    return parse_parameter_line_checked(line, name, name_size,
                                        value, value_size, is_data,
                                        "parameters.test", 7,
                                        errmsg, errmsg_size);
}

int main(void)
{
    char errmsg[512] = {0};
    char line[16];
    char name[64];
    char value[256];
    bool bool_value = FALSE;
    FILE *input;
    int has_line;
    int is_data;
    int i;
    unsigned long line_number;

    CHECK(parse_line("  sizeHistN\t =\t32 # bins\n",
                     name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == SUCCESS,
          "ordinary assignment was rejected");
    CHECK(is_data == TRUE && strcmp(name, "sizeHistN") == 0
          && strcmp(value, "32") == 0,
          "ordinary assignment was parsed incorrectly");

    CHECK(parse_line("outfile =   # use the default\n",
                     name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == SUCCESS
          && is_data == FALSE,
          "intentional empty value was not preserved");

    CHECK(parse_line("preScript = \"echo # keep % keep\" # comment\n",
                     name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == SUCCESS,
          "quoted comment characters were rejected");
    CHECK(is_data == TRUE && strcmp(name, "preScript") == 0
          && strcmp(value, "\"echo # keep % keep\"") == 0,
          "quoted comment characters were truncated");

    CHECK(parse_line("legacy-bare-option\n",
                     name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == SUCCESS
          && is_data == FALSE,
          "legacy bare line should remain ignorable");

    errmsg[0] = '\0';
    CHECK(parse_line("   = 1\n", name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == FAILURE,
          "missing parameter name was accepted");
    CHECK(strstr(errmsg, "no parameter name") != NULL,
          "missing-name diagnostic is unclear");

    errmsg[0] = '\0';
    CHECK(parse_line("preScript = \"unterminated\n",
                     name, sizeof(name), value, sizeof(value),
                     &is_data, errmsg, sizeof(errmsg)) == FAILURE,
          "unterminated quote was accepted");
    CHECK(strstr(errmsg, "unterminated quote") != NULL,
          "unterminated-quote diagnostic is unclear");

    CHECK(parse_bool_checked("YES", &bool_value, errmsg, sizeof(errmsg),
                             "usePeriodic") == SUCCESS
          && bool_value == TRUE,
          "valid true boolean was rejected");
    CHECK(parse_bool_checked("false", &bool_value, errmsg, sizeof(errmsg),
                             "usePeriodic") == SUCCESS
          && bool_value == FALSE,
          "valid false boolean was rejected");
    errmsg[0] = '\0';
    CHECK(parse_bool_checked("truejunk", &bool_value,
                             errmsg, sizeof(errmsg),
                             "usePeriodic") == FAILURE,
          "boolean with trailing junk was accepted");
    CHECK(strstr(errmsg, "invalid boolean value") != NULL,
          "invalid-boolean diagnostic is unclear");

    input = tmpfile();
    CHECK(input != NULL, "could not create overlong-line input");
    for (i = 0; i < 40; i++)
        fputc('x', input);
    fputc('\n', input);
    rewind(input);
    line_number = 0;
    errmsg[0] = '\0';
    CHECK(read_parameter_line_checked(input, line, sizeof(line),
                                      &has_line, &line_number,
                                      "parameters.test",
                                      errmsg, sizeof(errmsg)) == FAILURE,
          "overlong physical line was silently split");
    CHECK(strstr(errmsg, "exceeds") != NULL,
          "overlong-line diagnostic is unclear");
    fclose(input);

    input = tmpfile();
    CHECK(input != NULL, "could not create embedded-NUL input");
    fputs("name = value", input);
    fputc('\0', input);
    fputs("hidden\n", input);
    rewind(input);
    line_number = 0;
    errmsg[0] = '\0';
    CHECK(read_parameter_line_checked(input, line, sizeof(line),
                                      &has_line, &line_number,
                                      "parameters.test",
                                      errmsg, sizeof(errmsg)) == FAILURE,
          "embedded NUL was silently accepted");
    CHECK(strstr(errmsg, "embedded NUL") != NULL,
          "embedded-NUL diagnostic is unclear");
    fclose(input);

    input = tmpfile();
    CHECK(input != NULL, "could not create exact-width input");
    fputs("123456789012345\nnext = 1\n", input);
    rewind(input);
    line_number = 0;
    CHECK(read_parameter_line_checked(input, line, sizeof(line),
                                      &has_line, &line_number,
                                      "parameters.test",
                                      errmsg, sizeof(errmsg)) == SUCCESS
          && has_line == TRUE && strcmp(line, "123456789012345") == 0,
          "exact-width physical line was rejected");
    CHECK(read_parameter_line_checked(input, line, sizeof(line),
                                      &has_line, &line_number,
                                      "parameters.test",
                                      errmsg, sizeof(errmsg)) == SUCCESS
          && has_line == TRUE && strcmp(line, "next = 1\n") == 0,
          "reader lost synchronization after an exact-width line");
    fclose(input);

    puts("PASS: parameter parser rejects malformed input safely");
    return EXIT_SUCCESS;
}
