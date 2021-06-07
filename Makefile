# makefile
# author: Tom Ross
# email: tom.ross@protonmail.com
# github: dgsaf
#
# Derived from a makefile provided by Dr Elahi Pascal, in 20201, for
# phys4004 (high performance computing) - assignment 1, at Curtin university.


# makefile configuration settings

# possible flags for makefile configuration
OPT_COMPILERTYPE := GCC CLANG AOMP CRAY INTEL
OPT_GPUTYPE := NVIDIA AMD
OPT_OPTLEVEL := 0 1 2 3 fast s g
OPT_PROFILING := OFF ON

# select flags for make configuration
COMPILERTYPE ?= GCC
GPUTYPE ?= NVIDIA
OPTLEVEL ?= 2
PROFILING ?= OFF

# set optmisation flags
OPTFLAGS = -O$(OPTLEVEL)

# set profiling flags
ifeq ($(PROFILING), ON)
	OPTFLAGS += -pg -g
endif


# compiler selection

# define c compilers
GCC = gcc
CLANG = clang
AOMP = aompcc
CRAYCC = cc
INTELCC = icc

# define c++ compilers
GCCCXX = g++
CLANGCXX = clang++
AOMPCXX = aompcc
CRAYCXX = CC
INTELCXX = icpc

# define fortran compilers
GCCFORT = gfortran
CLANGFORT = flang
CRAYFORT = ftn
INTELFORT = ifort

# c compilation flags
GCCCFLAGS = -std=c11

# fortran compilation flags
GCCFFLAGS = -cpp  -dM -ffixed-line-length-none \
	-Wall -Wextra -Wconversion -pedantic -fimplicit-none -fcheck=all
INTELFFLAGS = -cpp -extend-source -D_INTELFTN

# define cuda compilers
NCC = nvcc
NCXX = nvcc
CUDA_FLAGS = -DUSECUDA

# define hip compilers
HCC = hipcc
HCXX = hipcc
HIP_FLAGS = -DUSEHIP
OCL_FLAGS = -lOpenCL -DUSEOPENCL

# define mpi compilers
MPICC = mpicc
MPICXX = mpicxx
MPIFORT = mpif90

# openmp compilation flags
GCCOMP_FLAGS = -fopenmp -DUSEOPENMP
GCCOMPTARGET_FLAGS = -fopenmp -DUSEOPENMPTARGET
INTELOMP_FLAGS = -qopenmp -DUSEOPENMP
INTELMPTARGET_FLAGS = -qopenmp -DUSEOPENMPTARGET
GCCOACC_FLAGS = -fopenacc -fopt-info-optimized-omp -DUSEOPENACC
INTELOACC_FLAGS = -qopenacc -fopt-info-optimized-omp -DUSEOPENACC

# select default compilers, but change compilers if specified
CC = $(GCC)
CXX = $(GCCCXX)
FORT = $(GCCFORT)
GPUCC = $(NCC)
GPUCXX = $(NCXX)
FFLAGS = $(GCCFFLAGS)
OMP_FLAGS = $(GCCOMP_FLAGS)
OMPTARGET_FLAGS = $(GCCOMPTARGET_FLAGS)
OACC_FLAGS = $(GCCOACC_FLAGS)
CFLAGS = $(GCCCFLAGS)

ifeq ($(COMPILERTYPE), CLANG)
	CC = $(CLANG)
	CXX = $(CLANGCXX)
	FORT = $(CLANGFORT)
endif

ifeq ($(COMPILERTYPE), AOMP)
	CC = $(AOMP)
	CXX = $(AOMPCXX)
endif

ifeq ($(COMPILERTYPE), CRAY)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	FORT = $(CRAYFORT)
	FFLAGS = -eZ -ffree
endif

ifeq ($(COMPILERTYPE), INTEL)
	CC = $(INTELCC)
	CXX = $(INTELCXX)
	FORT = $(INTELFORT)
	FFLAGS = $(INTELFFLAGS)
	OMP_FLAGS = $(INTELOMP_FLAGS)
	OMPTARGET_FLAGS = $(INTELOMPTARGET_FLAGS)
endif

ifeq ($(GPUTYPE), AMD)
	GPUCC = $(HCC)
	GPUCXX = $(HCXX)
endif

# select compiler used for openmp
OMPCC = $(CC)
OMPCXX = $(CXX)
OMPFORT = $(FORT)

# select compilers used for openacc
OACCCC = $(CC)
OACCCXX = $(CXX)
OACCFORT = $(FORT)

# select compilers used for opencl
OCLC = $(CCGPU)
OCLCXX = $(CXXGPU)

# flags common to all c, c++, fortran compilers
COMMONFLAGS = $(OPTFLAGS)


# command output formatting
NULL :=
TAB := $(NULL)- $(NULL)
PFX := $(NULL)> $(NULL)


# file expansions and related variables

# define source files
SRCS_C := $(wildcard src/*.c)
SRCS_H := $(wildcard src/*.h)
SRCS_CXX := $(wildcard src/*.cpp)
SRCS_HXX := $(wildcard src/*.hpp)
SRCS_F := $(wildcard src/*.f*) $(wildcard src/*.F*)
SRCS_F_MODS := $(wildcard src/m_*.f*) $(wildcard src/m_*.F*)
SRCS_F_PRGS := $(wildcard src/p_*.f*) $(wildcard src/p_*.F*)

# select fortran files (and c header files for preprocessing)
HDRS := $(SRCS_H)
SRCS := $(SRCS_F)
SRCS_MODS := $(SRCS_F_MODS)
SRCS_PRGS := $(SRCS_F_PRGS)

# define the base-names and file-suffixes of the source files
EXTS := $(suffix $(SRCS))
BSNS := $(notdir $(basename $(SRCS)))
BSNS_MODS := $(notdir $(basename $(SRCS_MODS)))
BSNS_PRGS := $(notdir $(basename $(SRCS_PRGS)))

# define the object and module files for each source file
OBJS := $(addprefix obj/,$(addsuffix .o,$(BSNS)))
MODS := $(addprefix mod/,$(addsuffix .mod,$(BSNS_MODS)))

# define binary targets to be made
BINS := $(addprefix bin/,$(BSNS_PRGS))


# directory commands

# make obj/, mod/, bin/ if they dont already exist
.PHONY : dirs
dirs :
	[ -d obj ] || mkdir obj
	[ -d mod ] || mkdir mod
	[ -d bin ] || mkdir bin


# clean commands

# all clean commands (obj/, mod/, bin/)
.PHONY : clean
clean : clean_obj clean_mod clean_bin

# clean obj/
.PHONY : clean_obj
clean_obj :
	rm -f obj/*

# clean mod/
.PHONY : clean_mod
clean_mod :
	rm -f mod/*

# clean bin/
.PHONY : clean_bin
clean_bin :
	rm -f bin/*


# make commands

# explicit target dependencies for objects
obj/m_random.o : obj/m_parameters.o
obj/m_io.o : obj/m_parameters.o
obj/m_diffeq.o : obj/m_parameters.o
obj/m_schrodinger.o : obj/m_parameters.o obj/m_diffeq.o obj/m_integrate.o

# implicit rule for arbitrary fortan targets
obj/%.o : $(firstword $(addprefix src/%,$(EXTS)))
	@echo "$(PFX)$@ : $^"
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $< -o $@ -J mod/


# explicit target dependencies for binaries
bin/p_test_io : obj/m_random.o obj/m_io.o

# implicit rule for binary targets
bin/% : $(firstword $(addprefix src/%,$(EXTS)))
	@echo "$(PFX)$@ : $^"
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -o $@ $^ -J mod/


# info commands

# all information targets
.PHONY : info
info : info_config info_build info_commands info_files

# information about makefile configuration
.PHONY : info_config
info_config :
	@echo "Makefile configuration options (with [selected]):"
	@echo "$(PFX)COMPILERTYPE ?=" \
		"$(patsubst $(COMPILERTYPE),[$(COMPILERTYPE)],$(OPT_COMPILERTYPE))"
	@echo "$(PFX)GPUTYPE      ?=" \
		"$(patsubst $(GPUTYPE),[$(GPUTYPE)],$(OPT_GPUTYPE))"
	@echo "$(PFX)OPTLEVEL     ?=" \
		"$(patsubst $(OPTLEVEL),[$(OPTLEVEL)],$(OPT_OPTLEVEL))"
	@echo "$(PFX)PROFILING    ?=" \
		"$(patsubst $(PROFILING),[$(PROFILING)],$(OPT_PROFILING))"

# information about current build given the commands make was passed
.PHONY : info_build
info_build :
	@echo "Compilers currently selected for (C, C++, FORTRAN) codes:"
	@echo "$(PFX)CPU         := ($(CC), $(CXX), $(FORT))"
	@echo "$(PFX)MPI         := ($(MPICC), $(MPICXX), $(MPIFORT))"
	@echo "$(PFX)GPU         := ($(GPUCC), $(GPUCXX),)"
	@echo "$(PFX)GPU-OpenMP  := ($(OMPCC), $(OMPCXX), $(OMPFORT))"
	@echo "$(PFX)GPU-OpenACC := ($(OACCCC), $(OACCCXX), $(OACCFORT))"

# information about source, base name, and object files
.PHONY : info_files
info_files :
	@echo "Source files (and their derived files) and binaries:"
	@echo "$(PFX)HDRS := $(HDRS)"
	@echo "$(PFX)SRCS := $(SRCS)"
	@echo "$(PFX)$(TAB)BSNS := $(BSNS)"
	@echo "$(PFX)$(TAB)EXTS := $(EXTS)"
	@echo "$(PFX)OBJS := $(OBJS)"
	@echo "$(PFX)MODS := $(MODS)"
	@echo "$(PFX)BINS := $(BINS)"

# information about current make commands available
.PHONY : info_commands
info_commands :
	@echo "Make commands:"
	@echo "$(PFX)info  : All info_* commands"
	@echo "$(PFX)$(TAB)info_config   : Configuration options for makefile"
	@echo "$(PFX)$(TAB)info_build    : The current compilers selected in makefile"
	@echo "$(PFX)$(TAB)info_dirs     : The contents of relevant directories"
	@echo "$(PFX)$(TAB)info_files    : The file expansions used in makefile"
	@echo "$(PFX)$(TAB)info_commands : The make commands available"
	@echo "$(PFX)dirs  : Creates obj/, mod/, bin/ if they dont exist"
	@echo "$(PFX)clean : All clean_* commands"
	@echo "$(PFX)$(TAB)clean_obj : Removes contents of obj/"
	@echo "$(PFX)$(TAB)clean_mod : Removes contents of mod/"
	@echo "$(PFX)$(TAB)clean_bin : Removes contents of bin/"
