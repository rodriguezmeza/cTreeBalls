"""Public source distributions use external GSL/CFITSIO, not omitted sources.

Only the staged release tree is rewritten. Atomic replacement is required:
setuptools may create that tree using hard links to the developer's files.
This module uses the standard library so archive inspection needs no compiler.
"""
import os
from pathlib import Path
import re
import stat
import tempfile

EXTERNAL_DEFAULTS = {
    'Makefile_settings': ('GSLINTERNAL', '0'),
    'addons/Makefile_addons_settings': ('CFITSIOLIBON', '0'),
}


def setting_pattern(key):
    return re.compile(r'^[ \t]*(?:override[ \t]+|export[ \t]+)?' + re.escape(key)
                      + r'[ \t]*(?::=|\?=|\+=|=)[ \t]*([^\n]*)$', re.MULTILINE)


def check_external_defaults(contents):
    """Check declared archive defaults without executing code from an archive."""
    result = {}
    for name, (key, expected) in EXTERNAL_DEFAULTS.items():
        matches = setting_pattern(key).findall(contents[name])
        # Make preserves whitespace before a comment and at the end of a value.
        # The build uses exact ifeq comparisons, so "0 # comment" is not "0".
        values = [value.split('#', 1)[0] for value in matches]
        if values != [expected]:
            raise AssertionError(f'external dependency profile requires exactly one '
                                 f'{key}={expected} in {name}; found {values}')
        result[key] = expected
    return result


def normalize_release_tree(root):
    """Detach and rewrite staged defaults; leave every other setting untouched."""
    root = Path(root)
    staged = []
    for name, (key, value) in EXTERNAL_DEFAULTS.items():
        path = root / name
        original = path.read_text(encoding='utf-8')
        pattern = setting_pattern(key)
        if len(pattern.findall(original)) != 1:
            raise RuntimeError(f'cannot package ambiguous/missing {key} in {name}')
        text = pattern.sub(f'{key} = {value}', original)
        staged.append((path, text))
    # Validate both before changing either staged file.
    check_external_defaults({str(path.relative_to(root)): text for path, text in staged})
    for path, text in staged:
        mode = stat.S_IMODE(path.stat().st_mode)
        fd, temporary = tempfile.mkstemp(prefix='.release-profile-', dir=path.parent)
        try:
            with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as stream:
                stream.write(text)
            os.chmod(temporary, mode)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
