#!/bin/bash

# ---------------------------------------
# Environment setup for Ursa (Intelllvm)
# Run as: source env_ursa_intel.sh
# ---------------------------------------

# Should work with UFS modulefiles
module purge
module use /scratch4/BMC/ufs-artic/Kristin.Barton/repos/kristinbarton/ufs-arctic-workflow/main/ufs-weather-model/modulefiles/
module load ufs_ursa.intelllvm.lua

# Set the compiler
export FC=mpiifx

# Use build directory
export MOD_FLAG="-module build"

export INCS="-I. \
-I${fms_ROOT}/include_r8 \
-I${esmf_ROOT}/include \
-I${netcdf_fortran_ROOT}/include \
-I${netcdf_c_ROOT}/include \
-I${hdf5_ROOT}/include \
-I${parallelio_ROOT}/include"

export LIBS="-L${fms_ROOT}/lib -lfms_r8 \
-L${esmf_ROOT}/lib -Wl,-rpath,${esmf_ROOT}/lib -lesmf \
-L${parallelio_ROOT}/lib -lpiof -lpioc \
-L${netcdf_fortran_ROOT}/lib -lnetcdff \
-L${netcdf_c_ROOT}/lib -lnetcdf \
-L${hdf5_ROOT}/lib -lhdf5_hl -lhdf5 -lz"

export FFLAGS="-traceback -real-size 64 -qopenmp"
export DEBUG_FFLAGS="-O0 -g -check all -link_mpi=dbg_mt -traceback -real-size 64 -qopenmp"

echo "Environment configured for Ursa (Intelllvm). You can now run 'make ggca'."
