# PYTHIA configuration file.
# Generated on Mon Mar 1 20:15:47 EST 2021 with the user supplied options:
# --prefix=/cvmfs/eic.opensciencegrid.org/packages/pythia8/8303/linux-ubuntu20.10-x86_64-gcc-10.2.0-pmdnmh3nhcjfggrsqhdnpcub3zypu364

# Install directory prefixes.
PREFIX_BIN=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/bin
PREFIX_INCLUDE=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/include
PREFIX_LIB=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/lib
PREFIX_SHARE=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=/cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/gcc/8.3.0.1-0a5ad/x86_64-centos7/bin/g++
#CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread -DGZIP
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=


EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/root6//bin/root-config
ROOT_BIN=/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/root6//bin/
ROOT_INCLUDE=-I/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/root6/include
ROOT_LIB=-L/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/root6/lib -Wl,-rpath,/cvmfs/eic.opensciencegrid.org/gcc-8.3/MCEG/releases/env/EIC2022a/root6/lib -lCore

YODA_USE=false
YODA_CONFIG=
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=
