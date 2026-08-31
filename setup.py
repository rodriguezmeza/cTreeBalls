#
# cyballs python module setup
# written by: Mario A. Rodriguez-Meza
# date: 15.02.2026
#

#B new imports
from setuptools import setup, Extension
from Cython.Distutils import build_ext as cython_build_ext

import os
#B to parse some constants from C headers
import re
#E
#B for using same macros as C
import shlex
#E
import subprocess
import sys

import numpy as nm
#E

#B for using with pip -m install...
install_requires = [
    "numpy>=1.22",
    "scipy",
]
#E

#B a lots of definitions
# to parse some constants from C headers
def read_define_int(header_path, macro_name):
    pattern = re.compile(r"^\s*#\s*define\s+" + re.escape(macro_name) + r"\s+(.+?)(?:\s*/[/*].*)?$")

    with open(header_path, "r") as f:
        for line in f:
            m = pattern.match(line)
            if not m:
                continue

            value = m.group(1).strip()
            value = value.split("//", 1)[0].strip()
            value = value.split("/*", 1)[0].strip()

            try:
                return int(value)
            except ValueError:
                raise RuntimeError(
                    f"Macro {macro_name} in {header_path} is not a plain integer: {value}"
                )

    raise RuntimeError(f"Macro {macro_name} not found in {header_path}")

#B if needed fftw3
def first_existing_fftw_path(kind):
    candidates = [
        "/opt/homebrew",
        "/usr/local",
        os.path.join(os.path.expanduser("~"), "local", "fftw3"),
        os.path.join(os.path.expanduser("~"), "NagBody_pkg", "local", "fftw3"),
    ]

    for base in candidates:
        if kind == "include":
            path = os.path.join(base, "include")
            if os.path.exists(os.path.join(path, "fftw3.h")):
                return path
        elif kind == "lib":
            path = os.path.join(base, "lib")
            if os.path.isdir(path):
                for name in os.listdir(path):
                    if name.startswith("libfftw3."):
                        return path
    return None
#E

#def read_makefile_variable(path, key, default=None):
#    try:
#        with open(path, "r") as f:
#            for line in f:
#                line = line.strip()
#                if line.startswith(key) and "=" in line:
#                    return line.split("=", 1)[1].strip()
#    except FileNotFoundError:
#        pass
#    return default


def get_gsl_config(flag):
    try:
        res = subprocess.check_output(["gsl-config", flag], text=True)
        return res.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return ""

def read_makefile_setting(path, key, default="0"):
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(key):
                    return line.split("=", 1)[1].strip()
    except FileNotFoundError:
        pass
    return os.environ.get(key, default)

# gsl definition
def parse_gsl_config():
    cflags = get_gsl_config("--cflags").split()
    libs_flags = get_gsl_config("--libs").split()

    include_dirs = [
        flag[2:] for flag in cflags
        if flag.startswith("-I") and len(flag) > 2
    ]

    library_dirs = [
        flag[2:] for flag in libs_flags
        if flag.startswith("-L") and len(flag) > 2
    ]

    libraries = [
        flag[2:] for flag in libs_flags
        if flag.startswith("-l") and len(flag) > 2
    ]

    return include_dirs, library_dirs, libraries
#

#B to get any pkg-config pkg
def parse_pkg_config(package):
    cflags = get_pkg_config(package, "--cflags").split()
    libs_flags = get_pkg_config(package, "--libs").split()

    include_dirs = [x[2:] for x in cflags if x.startswith("-I")]
    library_dirs = [x[2:] for x in libs_flags if x.startswith("-L")]
    libraries = [x[2:] for x in libs_flags if x.startswith("-l")]

    return include_dirs, library_dirs, libraries

def get_pkg_config(package, flag):
    try:
        return subprocess.check_output(
            ["pkg-config", flag, package],
            text=True
        ).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return ""
#E

# generator function
def indent_block(text, spaces=8):
    prefix = " " * spaces
    return "\n".join(prefix + line if line.strip() else line
                     for line in text.strip("\n").splitlines())

#E a lots of definitions

#B
# Recover the gcc compiler runtime path
gcc_lib = subprocess.check_output(
    ["gcc", "-print-libgcc-file-name"],
    text=True,
).strip()
GCCPATH = os.path.normpath(os.path.dirname(gcc_lib))
#E

# define absolute paths
root_folder = os.path.dirname(os.path.abspath(__file__))
include_folder = os.path.join(root_folder, "include")
general_libs_folder = os.path.join(root_folder, "general_libs")
getparam_folder = os.path.join(root_folder, "getparam")
source_folder = os.path.join(root_folder, "source")
class_lib_folder = os.path.join(os.path.join(root_folder, "addons"),"class_lib")
octree_ggg_omp_folder = os.path.join(os.path.join(root_folder, "addons"),"octree_ggg_omp")
pxd_folder = os.path.join(os.path.join(root_folder, "addons"),"pxd")
cyballs_folder = os.path.join(root_folder, "python")
#
addons_folder = os.path.join(root_folder, "addons")
addons_include_folder = os.path.join(addons_folder, "addons_include")
addons_include_include_folder = os.path.join(addons_include_folder, "include")
addons_include_addons_folder = os.path.join(addons_include_folder, "addons")
addons_include_startrun_folder = os.path.join(addons_include_folder, "source", "startrun")
addons_include_cballs_folder = os.path.join(addons_include_folder, "source", "cballs")
addons_include_cballsio_folder = os.path.join(addons_include_folder, "source", "cballsio")
#

makefile_settings = os.path.join(root_folder, "Makefile_settings")
addons_settings = os.path.join(root_folder, "addons", "Makefile_addons_settings")
#makefile_machine = os.path.join(root_folder, "Makefile_machine")

#if "CC" not in os.environ:
#    make_cc = read_makefile_variable(makefile_machine, "CC", None)
#    if make_cc:
#        os.environ["CC"] = make_cc

#B GSL...
gsl_include_dirs, gsl_library_dirs, gsl_libraries = parse_gsl_config()

if os.environ.get("GSL_INCLUDE"):
    gsl_include_dirs = [os.environ["GSL_INCLUDE"]]

if os.environ.get("GSL_LIB"):
    gsl_library_dirs = [os.environ["GSL_LIB"]]

if not gsl_include_dirs:
    gsl_include_dirs = ["/usr/local/include"]

if not gsl_library_dirs:
    gsl_library_dirs = ["/usr/local/lib"]

if not gsl_libraries:
    gsl_libraries = ["gsl", "gslcblas", "m"]

#B to fix NDIM problem, the long term version...
def read_make_cyballs_env():
    try:
        output = subprocess.check_output(
            ["make", "--no-print-directory", "-s", "print-cyballs-build-env"],
            cwd=root_folder,
            text=True,
            stderr=subprocess.STDOUT,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "could not query Makefile build flags; build with `make cyballs` "
            "or fix the `print-cyballs-build-env` target"
        ) from exc

    values = {}
    for line in output.splitlines():
        key, sep, value = line.partition("=")
        if sep and key.startswith("__CBALLS_"):
            values[key] = value

    required = ["__CBALLS_LIB__", "__CBALLS_CPPFLAGS__"]
    missing = [key for key in required if key not in values]
    if missing:
        raise RuntimeError(
            "Makefile print-cyballs-build-env did not provide: "
            + ", ".join(missing)
        )

    return values

make_env = {}

#B==============================================
make_env = read_make_cyballs_env()

if "CC" not in os.environ and "__CBALLS_CC__" in make_env:
    os.environ["CC"] = make_env["__CBALLS_CC__"]

libname = os.environ.get("CBALLS_LIB") or make_env["__CBALLS_LIB__"]
make_cppflags = shlex.split(
    os.environ.get("CBALLS_CPPFLAGS") or make_env["__CBALLS_CPPFLAGS__"]
)

#B
def cppflag_has_macro(flags, name):
    prefix = f"-D{name}="
    return any(flag == f"-D{name}" or flag.startswith(prefix) for flag in flags)

ADDONSON = "1" if cppflag_has_macro(make_cppflags, "ADDONS") else "0"
CLASSLIBON = "1" if cppflag_has_macro(make_cppflags, "CLASSLIB") else "0"
PXDON = "1" if cppflag_has_macro(make_cppflags, "PXD") else "0"
CFITSIOON = "1" if cppflag_has_macro(make_cppflags, "CFITSIO") else "0"
CFITSIOLIBON = "1" if cppflag_has_macro(make_cppflags, "CFITSIOLIB") else "0"
USEGSL = "1" if cppflag_has_macro(make_cppflags, "USEGSL") else "0"
BALLTREEMPION = "1" if cppflag_has_macro(make_cppflags, "BALLTREEMPI") else "0"
BALLTREE2BALLSMPION = (
    "1" if cppflag_has_macro(make_cppflags, "BALLTREE2BALLSMPI") else "0"
)
BALLTREE2BALLSOMP3PCFON = (
    "1" if cppflag_has_macro(make_cppflags, "BALLTREE2BALLSOMP3PCF") else "0"
)
BALLTREE2BALLSMPI3PCFON = (
    "1" if cppflag_has_macro(make_cppflags, "BALLTREE2BALLSMPI3PCF") else "0"
)
OCTREE2BALLSOMPON = (
    "1" if cppflag_has_macro(make_cppflags, "OCTREE2BALLSOMP") else "0"
)
OCTREE2BALLSMPION = (
    "1" if cppflag_has_macro(make_cppflags, "OCTREE2BALLSMPI") else "0"
)
OCTREEGGGMPION = "1" if cppflag_has_macro(make_cppflags, "OCTREEGGGMPI") else "0"
OCTREEBALLS4MPION = "1" if cppflag_has_macro(make_cppflags, "OCTREEBALLS4MPI") else "0"
LYAFORESTOMPON = "1" if cppflag_has_macro(make_cppflags, "LYAFORESTOMP") else "0"
LYAFORESTMPION = "1" if cppflag_has_macro(make_cppflags, "LYAFORESTMPI") else "0"
OCTREE3PCF3DOMPON = "1" if cppflag_has_macro(make_cppflags, "OCTREE3PCF3DOMP") else "0"
OCTREE3PCF3DMPION = "1" if cppflag_has_macro(make_cppflags, "OCTREE3PCF3DMPI") else "0"
OCTREESHEAROMPON = "1" if cppflag_has_macro(make_cppflags, "OCTREESHEAROMP") else "0"

GSLINTERNAL = make_env.get("__CBALLS_GSLINTERNAL__", "1")
OPENMPMACHINE = make_env.get("__CBALLS_OPENMPMACHINE__", "0")
SINGLEPON = make_env.get("__CBALLS_SINGLEPON__", "0")
macos_deployment_target = make_env.get(
    "__CBALLS_MACOSX_DEPLOYMENT_TARGET__",
    os.environ.get("MACOSX_DEPLOYMENT_TARGET", ""),
)

if sys.platform == "darwin" and macos_deployment_target:
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = macos_deployment_target

if USEGSL != "1":
    raise RuntimeError("cyballs requires USEGSL=1")

define_macros = [
    ("__CTREEBALLSDIR__", f'"{root_folder}"'),
    ("USEGSL", None),
    ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
]

if ADDONSON == "1":
    define_macros.append(("ADDONS", None))
if CLASSLIBON == "1":
    define_macros.append(("CLASSLIB", None))
if PXDON == "1":
    define_macros.append(("PXD", None))
#E

# additional
openmp_compile_args = []
openmp_link_args = []
macos_compile_args = []
macos_link_args = []

if sys.platform == "darwin" and macos_deployment_target:
    deployment_flag = f"-mmacosx-version-min={macos_deployment_target}"
    macos_compile_args.append(deployment_flag)
    macos_link_args.append(deployment_flag)

if OPENMPMACHINE == "1":
    openmp_compile_args += ["-fopenmp"]
    openmp_link_args += ["-fopenmp", "-lgomp"]

#B
if GSLINTERNAL == "1":
    gsl_include_dirs = [
        os.path.join(root_folder, "addons", "gsl"),
        os.path.join(root_folder, "addons", "gsl", "gsl"),
    ]
    gsl_library_dirs = []
    gsl_libraries = []
else:
    define_macros.append(("NOINTERNALGSL", None))

liblist = [libname] + gsl_libraries

if sys.platform.startswith("linux"):
    liblist += ["mvec"]
#E

#E==============================================

#B==============================================
#
# Parse C ABI constants used to generate the Cython declarations.
#
common_defs_h = os.path.join(include_folder, "common_defs.h")
global_data_h = os.path.join(include_folder, "global_data.h")
class_common_h = os.path.join(class_lib_folder, "common.h")
parser_h = os.path.join(class_lib_folder, "parser.h")

filearg_size = read_define_int(parser_h, "_ARGUMENT_LENGTH_MAX_")
errormsg_size = read_define_int(class_common_h, "_ERRORMSGSIZE_")
filename_size = (
    read_define_int(class_common_h, "_FILENAMESIZE_")
    + read_define_int(class_common_h, "_BASEPATHSIZE_")
)

def cpp_macro_defined(name):
    return (
        any(flag == f"-D{name}" or flag.startswith(f"-D{name}=") for flag in make_cppflags)
        or any(macro[0] == name for macro in define_macros)
    )

def cpp_macro_value(name):
    prefix = f"-D{name}="
    for flag in make_cppflags:
        if flag.startswith(prefix):
            return flag[len(prefix):]
    return None

def active_ndim():
    value = cpp_macro_value("NDIM")
    if value is not None:
        return int(value)
    if cpp_macro_defined("ONEDIMCODE") or cpp_macro_defined("ONEDIM"):
        return 1
    if cpp_macro_defined("THREEDIMCODE") or cpp_macro_defined("THREEDIM"):
        return 3
    return 2

PXD_CONSTANTS = {
    "FILEARG_SIZE": str(filearg_size),
    "ERRORMSG_SIZE": str(errormsg_size),
    "FILENAME_SIZE": str(filename_size),
    "MAXITEMS": str(read_define_int(common_defs_h, "MAXITEMS")),
    "MAXLENGTHOFFILES": str(read_define_int(common_defs_h, "MAXLENGTHOFFILES")),
    "MAXLENGTHOFFMTFILES": str(read_define_int(common_defs_h, "MAXLENGTHOFFMTFILES")),
    "MAXLENGTHOFINDIVIDUALFILES": str(read_define_int(common_defs_h, "MAXLENGTHOFINDIVIDUALFILES")),
    "MAXLEVEL": str(read_define_int(global_data_h, "MAXLEVEL")),
    "NDIM": str(active_ndim()),
    "INTEGER_CTYPE": "long" if cpp_macro_defined("LONGINT") else "int",
}
#E==============================================

#B==============================================
def generate_ccyballs_pxd():
    template_path = os.path.join(cyballs_folder, "ccyballs.pxd.in")
    output_path = os.path.join(cyballs_folder, "ccyballs.pxd")

#B
    smoothpivot_hist_fields = ""
    smoothpivot_n2pcf_fields = ""
    balls4scanlev_fields = ""
    pivotcount_fields = ""

    if cpp_macro_defined("SMOOTHPIVOT"):
        smoothpivot_hist_fields = "double * histNNSubXi2pcftotal"
        smoothpivot_n2pcf_fields = "double * histNNSubN2pcftotal"

    if (cpp_macro_defined("BALLS4SCANLEV")
            or cpp_macro_defined("OCTREEBALLS4OMP")
            or cpp_macro_defined("OCTREEBALLS4MPI")):
        balls4scanlev_fields = "unsigned char flagBalls4Scanlevel"

    if cpp_macro_defined("NMultipoles") and cpp_macro_defined("NONORMHIST"):
        pivotcount_fields = "long pivotCount"
#E

    cmdline_classlib_fields = ""
    cmdline_cosmolib_fields = ""
    global_classlib_fields = ""
    global_cosmolib_fields = ""

    basepath_size = read_define_int(common_defs_h, "_BASEPATHSIZE_")

    if ADDONSON == "1" and CLASSLIBON == "1":
        cmdline_classlib_fields = indent_block(f"""
char base_path[{basepath_size}]
""")
        global_classlib_fields = ""


    if ADDONSON == "1":
        cmdline_cosmolib_fields = indent_block("""
""")

        global_cosmolib_fields = indent_block("""
""")

    with open(template_path, "r") as f:
        text = f.read()

    text = text.replace("{{CMDLINE_CLASSLIB_FIELDS}}", cmdline_classlib_fields)
    text = text.replace("{{CMDLINE_COSMOLIB_FIELDS}}", cmdline_cosmolib_fields)
    text = text.replace("{{GLOBAL_CLASSLIB_FIELDS}}", global_classlib_fields)
    text = text.replace("{{GLOBAL_COSMOLIB_FIELDS}}", global_cosmolib_fields)
#
    text = text.replace("{{GLOBAL_SMOOTHPIVOT_HIST_FIELDS}}", smoothpivot_hist_fields)
    text = text.replace("{{GLOBAL_SMOOTHPIVOT_N2PCF_FIELDS}}", smoothpivot_n2pcf_fields)
    text = text.replace("{{GLOBAL_BALLS4SCANLEV_FIELDS}}", balls4scanlev_fields)
    text = text.replace("{{GLOBAL_PIVOTCOUNT_FIELDS}}", pivotcount_fields)
#
    for key, value in PXD_CONSTANTS.items():
        text = text.replace("{{" + key + "}}", value)

    with open(output_path, "w") as f:
        f.write(text)
#E==============================================

#B==============================================

# This makes the .pxd match the active C macros.
if ADDONSON != "1" or CLASSLIBON != "1" or PXDON != "1":
    raise RuntimeError("cyballs requires ADDONSON=1, CLASSLIBON=1, PXDON=1")
#

# It can be necessary to export path to pkgconfig
#export PKG_CONFIG_PATH=/Users/mar/NagBody_pkg/local/cfitsio/lib/pkgconfig:$PKG_CONFIG_PATH
#cfitsio_include_dirs, cfitsio_library_dirs, cfitsio_libraries = parse_pkg_config("cfitsio")
#
#cfitsio_rpath_args = [
#    f"-Wl,-rpath,{libdir}"
#    for libdir in cfitsio_library_dirs
#]

#B CFITSIO configuration
cfitsio_include_dirs = []
cfitsio_library_dirs = []
cfitsio_libraries = []
cfitsio_rpath_args = []

if CFITSIOON == "1":
    define_macros.append(("CFITSIO", None))

    if CFITSIOLIBON == "1":
        define_macros.append(("CFITSIOLIB", None))
        cfitsio_include_dirs = [
            os.path.join(root_folder, "addons", "cfitsiolib", "cfitsio_ver_4.6.3")
        ]
        # No external libcfitsio: objects should be inside libcballs.a
    else:
        cfitsio_include_dirs, cfitsio_library_dirs, cfitsio_libraries = parse_pkg_config("cfitsio")
        cfitsio_rpath_args = [
            f"-Wl,-rpath,{libdir}"
            for libdir in cfitsio_library_dirs
        ]
#E


with open(os.path.join(class_lib_folder, 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            VERSION = line.split()[-1][2:-1]
            break

cyballs_ext = Extension(
    "cyballs",
        [
            os.path.join(cyballs_folder, "cyballs.pyx"),
            os.path.join(source_folder, "abi_check.c"),
        ],
    depends=[os.path.join(root_folder, f"lib{libname}.a")],
    include_dirs=[
        nm.get_include(),
        include_folder,
        class_lib_folder,
        general_libs_folder,
        getparam_folder,
        pxd_folder,
        source_folder,
        #
        addons_folder,
        addons_include_folder,
        addons_include_include_folder,
        addons_include_addons_folder,
        addons_include_startrun_folder,
        addons_include_cballs_folder,
        addons_include_cballsio_folder,
        #
        *gsl_include_dirs,
        *cfitsio_include_dirs,
    ],
    define_macros=define_macros,
    libraries=liblist + cfitsio_libraries,
    library_dirs=[
        root_folder,
        GCCPATH,
        *gsl_library_dirs,
        *cfitsio_library_dirs,
    ],
    extra_link_args=(openmp_link_args + macos_link_args
                     + ['-lz'] + cfitsio_rpath_args),
    extra_compile_args=(openmp_compile_args + macos_compile_args
                        + make_cppflags),
)

cyballs_ext.cython_directives = {
    'language_level': "3" if sys.version_info.major >= 3 else "2"
}
#E==============================================

#B==============================================
class build_ext(cython_build_ext):
    def run(self):
        generate_ccyballs_pxd()
        if os.environ.get("CBALLS_STATIC_LIBRARY_READY") != "1":
            make_profile = {
                "ADDONSON": ADDONSON,
                "CLASSLIBON": CLASSLIBON,
                "PXDON": PXDON,
                "CFITSIOON": CFITSIOON,
                "CFITSIOLIBON": CFITSIOLIBON,
                "USEGSL": USEGSL,
                "GSLINTERNAL": GSLINTERNAL,
                "OPENMPMACHINE": OPENMPMACHINE,
                "SINGLEPON": SINGLEPON,
                "BALLTREE2BALLSMPION": BALLTREE2BALLSMPION,
                "BALLTREE2BALLSOMP3PCFON": BALLTREE2BALLSOMP3PCFON,
                "BALLTREE2BALLSMPI3PCFON": BALLTREE2BALLSMPI3PCFON,
                "OCTREE2BALLSOMPON": OCTREE2BALLSOMPON,
                "OCTREE2BALLSMPION": OCTREE2BALLSMPION,
                "BALLTREEMPION": BALLTREEMPION,
                "OCTREEGGGMPION": OCTREEGGGMPION,
                "OCTREEBALLS4MPION": OCTREEBALLS4MPION,
                "LYAFORESTOMPON": LYAFORESTOMPON,
                "LYAFORESTMPION": LYAFORESTMPION,
                "OCTREE3PCF3DOMPON": OCTREE3PCF3DOMPON,
                "OCTREE3PCF3DMPION": OCTREE3PCF3DMPION,
                "OCTREESHEAROMPON": OCTREESHEAROMPON,
            }
            make_command = ["make"] + [
                f"{key}={value}" for key, value in make_profile.items()
            ] + ["cyballs-static-lib"]
            subprocess.check_call(make_command, cwd=root_folder)
        super().run()

#B for using with pip -m install...
setup(
    name='cyballs',
    version=VERSION,
    description='Python interface to the nPCF code cballs',
    url='http://github.com/rodriguezmeza/cTreeBalls.git',
    cmdclass={'build_ext': build_ext},
    ext_modules=[cyballs_ext],
    install_requires=install_requires,
)
#E
#E==============================================
