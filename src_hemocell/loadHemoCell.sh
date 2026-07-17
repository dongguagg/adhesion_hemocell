#!/usr/bin/env bash

# ============================================================
# HemoCell environment for jxh
#
# Assumption:
#   This file is sourced in a fresh SSH shell before compiling
#   or running HemoCell.
#
# Target:
#   System GCC/G++ + OpenMPI + system HDF5 + bundled Palabos
#
# HemoCell:
#   /home/jxh/adhesion_rbc/src_hemocell
#
# Palabos:
#   /home/jxh/adhesion_rbc/src_hemocell/palabos
#
# Usage:
#   source ~/loadHemoCell.sh
#
# Then:
#   cd $HEMOCELL_DIR
#   rm -rf build
#   mkdir build
#   cd build
#   cmake ..
#   make -j
#
# Do NOT run:
#   bash ~/loadHemoCell.sh
# ============================================================


# ------------------------------------------------------------
# Must be sourced
# ------------------------------------------------------------

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: Please source this file instead of executing it:"
    echo "  source ~/loadHemoCell.sh"
    exit 1
fi


# ------------------------------------------------------------
# User paths
# ------------------------------------------------------------

export HEMOCELL_DIR=/home/jxh/adhesion_rbc/src_hemocell
export HEMOCELL_ROOT="$HEMOCELL_DIR"

export PALABOS_DIR=/home/jxh/adhesion_rbc/src_hemocell/palabos
export PALABOS_ROOT="$PALABOS_DIR"
export PLB_ROOT="$PALABOS_DIR"
export PALABOS="$PALABOS_DIR"


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

path_prepend() {
    local var="$1"
    local dir="$2"
    local old_value
    local new_value

    [ -d "$dir" ] || return 0

    old_value="${!var-}"

    new_value="$(
        printf '%s' "$old_value" |
        awk -v RS=: -v ORS=: -v dir="$dir" '
            length($0) > 0 && $0 != dir {print}
        ' |
        sed 's/:$//'
    )"

    if [ -n "$new_value" ]; then
        printf -v "$var" '%s:%s' "$dir" "$new_value"
    else
        printf -v "$var" '%s' "$dir"
    fi

    export "$var"
}

check_cmd() {
    local cmd="$1"

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: cannot find required command: $cmd"
        _hemocell_missing=1
    fi
}

check_mpi_cxx_cmd() {
    if command -v mpicxx >/dev/null 2>&1; then
        return 0
    fi

    if command -v mpic++ >/dev/null 2>&1; then
        return 0
    fi

    echo "ERROR: cannot find required MPI C++ wrapper: mpicxx or mpic++"
    _hemocell_missing=1
    return 1
}


# ------------------------------------------------------------
# Protect caller shell options
# ------------------------------------------------------------

_hemocell_old_opts="$-"

set +e
set +u


# ------------------------------------------------------------
# Basic system toolchain
# ------------------------------------------------------------

# Keep the existing PATH, but make sure system tools have priority.
path_prepend PATH /usr/local/sbin
path_prepend PATH /usr/local/bin
path_prepend PATH /usr/sbin
path_prepend PATH /usr/bin
path_prepend PATH /sbin
path_prepend PATH /bin

# Use system GCC/G++ explicitly.
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

# Tell OpenMPI wrappers to use system GCC/G++.
export OMPI_CC="$CC"
export OMPI_CXX="$CXX"

# Very important for your OpenMPI:
# /opt/software/openmpi-4.0.2 does not provide libmpi_cxx.
# Therefore, disable deprecated OpenMPI C++ bindings.
case " ${CXXFLAGS:-} " in
    *" -DOMPI_SKIP_MPICXX "*) ;;
    *) export CXXFLAGS="${CXXFLAGS:+$CXXFLAGS }-DOMPI_SKIP_MPICXX" ;;
esac

# MPI wrapper compilers.
if command -v mpicc >/dev/null 2>&1; then
    export MPICC="$(command -v mpicc)"
fi

if command -v mpicxx >/dev/null 2>&1; then
    export MPICXX="$(command -v mpicxx)"
elif command -v mpic++ >/dev/null 2>&1; then
    export MPICXX="$(command -v mpic++)"
fi

# Common Ubuntu/Debian OpenMPI library/include paths.
if [ -d /usr/lib/x86_64-linux-gnu/openmpi/lib ]; then
    path_prepend LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/openmpi/lib
    path_prepend LIBRARY_PATH /usr/lib/x86_64-linux-gnu/openmpi/lib
fi

if [ -d /usr/lib/x86_64-linux-gnu/openmpi/include ]; then
    path_prepend CPATH /usr/lib/x86_64-linux-gnu/openmpi/include
fi


# ------------------------------------------------------------
# Python packages
# ------------------------------------------------------------

# Use Ubuntu apt Python packages first.
# This avoids mixing apt h5py 2.10.0 with pip numpy 1.24.x from /usr/local.
if [ -d /usr/lib/python3/dist-packages ]; then
    path_prepend PYTHONPATH /usr/lib/python3/dist-packages
fi


# ------------------------------------------------------------
# HDF5
# ------------------------------------------------------------

# HemoCell normally should use parallel HDF5 with OpenMPI.
# On Ubuntu this usually comes from:
#   libhdf5-openmpi-dev

export HDF5_ROOT=/usr
export HDF5_PREFER_PARALLEL=TRUE
export HDF5_USE_STATIC_LIBRARIES=OFF

if [ -d /usr/include/hdf5/openmpi ]; then
    export HDF5_FLAVOR=openmpi
    path_prepend CPATH /usr/include/hdf5/openmpi
fi

if [ -d /usr/lib/x86_64-linux-gnu/hdf5/openmpi ]; then
    export HDF5_FLAVOR=openmpi
    path_prepend LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/openmpi
    path_prepend LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/openmpi
    path_prepend PKG_CONFIG_PATH /usr/lib/x86_64-linux-gnu/hdf5/openmpi/pkgconfig
    path_prepend CMAKE_PREFIX_PATH /usr/lib/x86_64-linux-gnu/hdf5/openmpi
elif [ -d /usr/lib/x86_64-linux-gnu/hdf5/serial ]; then
    export HDF5_FLAVOR=serial
    path_prepend CPATH /usr/include/hdf5/serial
    path_prepend LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial
    path_prepend LIBRARY_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial
    path_prepend PKG_CONFIG_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial/pkgconfig
    path_prepend CMAKE_PREFIX_PATH /usr/lib/x86_64-linux-gnu/hdf5/serial
fi

if command -v h5pcc >/dev/null 2>&1; then
    export HDF5_C_COMPILER="$(command -v h5pcc)"
elif command -v h5cc >/dev/null 2>&1; then
    export HDF5_C_COMPILER="$(command -v h5cc)"
fi


# ------------------------------------------------------------
# Palabos
# ------------------------------------------------------------

path_prepend CMAKE_PREFIX_PATH "$PALABOS_DIR"

if [ -d "$PALABOS_DIR/src" ]; then
    path_prepend CPATH "$PALABOS_DIR/src"
fi

if [ -d "$PALABOS_DIR/lib" ]; then
    path_prepend LD_LIBRARY_PATH "$PALABOS_DIR/lib"
    path_prepend LIBRARY_PATH "$PALABOS_DIR/lib"
fi

if [ -d "$PALABOS_DIR/build" ]; then
    path_prepend CMAKE_PREFIX_PATH "$PALABOS_DIR/build"
fi


# ------------------------------------------------------------
# ParMETIS / METIS
# ------------------------------------------------------------

# You said ParMETIS is not needed.
# In a fresh SSH shell these variables usually do not exist,
# but unsetting them is harmless and prevents accidental detection.
unset PARMETIS_ROOT
unset ParMETIS_ROOT
unset PARMETIS_DIR
unset ParMETIS_DIR
unset METIS_ROOT
unset METIS_DIR


# ------------------------------------------------------------
# Final checks
# ------------------------------------------------------------

hash -r

_hemocell_missing=0

if [ ! -d "$HEMOCELL_DIR" ]; then
    echo "ERROR: Cannot find HemoCell directory:"
    echo "  $HEMOCELL_DIR"
    _hemocell_missing=1
fi

if [ ! -d "$PALABOS_DIR" ]; then
    echo "ERROR: Cannot find Palabos directory:"
    echo "  $PALABOS_DIR"
    _hemocell_missing=1
fi

check_cmd gcc
check_cmd g++
check_cmd mpicc
check_mpi_cxx_cmd
check_cmd mpirun
check_cmd cmake
check_cmd patch
check_cmd python3

if ! command -v h5pcc >/dev/null 2>&1 && ! command -v h5cc >/dev/null 2>&1; then
    echo "ERROR: Cannot find h5pcc or h5cc."
    echo "On Ubuntu, you may need:"
    echo "  sudo apt install libhdf5-openmpi-dev hdf5-tools"
    _hemocell_missing=1
fi

_h5py_version="$(python3 -c 'import h5py; print(h5py.__version__)' 2>/dev/null)"

if [ -z "$_h5py_version" ]; then
    echo "ERROR: Python package h5py is not available for python3."
    echo "On Ubuntu, you may need:"
    echo "  sudo apt install python3-h5py"
    _hemocell_missing=1
fi

_hdf5_version=""
_hdf5_parallel=""

if [ -n "${HDF5_C_COMPILER:-}" ] && [ -x "$HDF5_C_COMPILER" ]; then
    _hdf5_version="$(
        "$HDF5_C_COMPILER" -showconfig 2>/dev/null |
        awk -F: '/HDF5 Version/ {gsub(/^[ \t]+/, "", $2); print $2; exit}'
    )"

    _hdf5_parallel="$(
        "$HDF5_C_COMPILER" -showconfig 2>/dev/null |
        awk -F: '/Parallel HDF5/ {gsub(/^[ \t]+/, "", $2); print $2; exit}'
    )"
fi

_mpi_cxx_show=""
_mpi_cxx_link=""

if command -v mpicxx >/dev/null 2>&1; then
    _mpi_cxx_show="$(mpicxx --showme:command 2>/dev/null || true)"
    _mpi_cxx_link="$(mpicxx --showme:link 2>/dev/null || true)"
elif command -v mpic++ >/dev/null 2>&1; then
    _mpi_cxx_show="$(mpic++ --showme:command 2>/dev/null || true)"
    _mpi_cxx_link="$(mpic++ --showme:link 2>/dev/null || true)"
fi

export HEMOCELL_CMAKE_BIN="$(command -v cmake || true)"


# ------------------------------------------------------------
# Report
# ------------------------------------------------------------

echo "============================================================"
echo "HemoCell environment loaded"
echo "============================================================"
echo "HEMOCELL_DIR=$HEMOCELL_DIR"
echo "PALABOS_DIR=$PALABOS_DIR"
echo "PALABOS_ROOT=$PALABOS_ROOT"
echo "PLB_ROOT=$PLB_ROOT"
echo ""
echo "CC=$CC"
echo "CXX=$CXX"
echo "CXXFLAGS=$CXXFLAGS"
echo "OMPI_CC=$OMPI_CC"
echo "OMPI_CXX=$OMPI_CXX"
echo "MPICC=${MPICC:-}"
echo "MPICXX=${MPICXX:-}"
echo ""
echo "gcc      = $(command -v gcc || true)"
echo "g++      = $(command -v g++ || true)"
echo "mpicc    = $(command -v mpicc || true)"
echo "mpicxx   = $(command -v mpicxx || true)"
echo "mpic++   = $(command -v mpic++ || true)"
echo "mpirun   = $(command -v mpirun || true)"
echo "cmake    = ${HEMOCELL_CMAKE_BIN:-}"
echo "patch    = $(command -v patch || true)"
echo "python3  = $(command -v python3 || true)"
echo "h5pcc    = $(command -v h5pcc || true)"
echo "h5cc     = $(command -v h5cc || true)"
echo ""
echo "gcc version       = $(gcc -dumpfullversion -dumpversion 2>/dev/null || true)"
echo "g++ version       = $(g++ -dumpfullversion -dumpversion 2>/dev/null || true)"
echo "mpirun version    = $(mpirun --version 2>/dev/null | head -n 1 || true)"
echo "MPI C++ compiler  = ${_mpi_cxx_show:-unknown}"
echo "MPI C++ link      = ${_mpi_cxx_link:-unknown}"
echo "cmake version     = $(cmake --version 2>/dev/null | head -n 1 || true)"
echo "HDF5 flavor       = ${HDF5_FLAVOR:-unknown}"
echo "HDF5 compiler     = ${HDF5_C_COMPILER:-unknown}"
echo "HDF5 version      = ${_hdf5_version:-unknown}"
echo "Parallel HDF5     = ${_hdf5_parallel:-unknown}"
echo "h5py version      = ${_h5py_version:-missing}"
echo "============================================================"


# ------------------------------------------------------------
# Make plain `cmake ..` work for this HemoCell environment
# ------------------------------------------------------------

# Your OpenMPI wrapper does not link libmpi_cxx, and your OpenMPI prefix
# does not contain libmpi_cxx. Therefore, CMake's old MPI C++ binding test
# must be skipped. This function only injects options for configure-style
# cmake calls, such as:
#
#   cmake ..
#
# It does not inject options into:
#
#   cmake --build .
#   cmake --install .
#   cmake --version
#   cmake -E ...
#   cmake -P ...
#
cmake() {
    local _hc_arg
    local _hc_is_configure

    if [ -z "${HEMOCELL_CMAKE_BIN:-}" ]; then
        echo "ERROR: HEMOCELL_CMAKE_BIN is not set."
        return 1
    fi

    if [ "$#" -eq 0 ]; then
        "$HEMOCELL_CMAKE_BIN"
        return $?
    fi

    _hc_is_configure=1

    for _hc_arg in "$@"; do
        case "$_hc_arg" in
            --build|--install|-E|-P|--version|-version|--help|--help-*|--system-information)
                _hc_is_configure=0
                break
                ;;
        esac
    done

    if [ "$_hc_is_configure" -eq 1 ]; then
        "$HEMOCELL_CMAKE_BIN" \
            "-DCMAKE_C_COMPILER=${CC}" \
            "-DCMAKE_CXX_COMPILER=${CXX}" \
            "-DMPI_C_COMPILER=${MPICC:-mpicc}" \
            "-DMPI_CXX_COMPILER=${MPICXX:-mpicxx}" \
            -DMPI_CXX_SKIP_MPICXX=TRUE \
            "-DCMAKE_CXX_FLAGS=${CXXFLAGS:-}" \
            "$@"
    else
        "$HEMOCELL_CMAKE_BIN" "$@"
    fi
}


# ------------------------------------------------------------
# Restore caller shell options
# ------------------------------------------------------------

case "$_hemocell_old_opts" in
    *e*) set -e ;;
    *)   set +e ;;
esac

case "$_hemocell_old_opts" in
    *u*) set -u ;;
    *)   set +u ;;
esac

unset _hemocell_old_opts
unset _h5py_version
unset _hdf5_version
unset _hdf5_parallel
unset _mpi_cxx_show
unset _mpi_cxx_link

if [ "$_hemocell_missing" -ne 0 ]; then
    echo "ERROR: HemoCell environment has missing required components."
    unset _hemocell_missing
    return 1 2>/dev/null || exit 1
fi

unset _hemocell_missing
