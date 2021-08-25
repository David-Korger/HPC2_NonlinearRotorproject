# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Csparse/Include -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 cs_sed.o \
 hpc_buildS_byInterpolation.o \
 hpc_cumsum.o \
 hpc_fmg.o  \
 hpc_malloc.o \
 hpc_mg.o  \
 hpc_mg_cycle.o  \
 hpc_print_vector.o \
 hpc_prol_quad.o  \
 hpc_prol.o  \
 hpc_rest.o  \
 hpc_solve_nonlin.o \
 hpc_solve_nonlin_damped.o \
 hpc_solve_nonlin_broyden.o \
 hpc_solve_nonlin_broyden_2.o \
 mesh_buildRhs_nonlin.o \
 mesh_buildRhs.o \
 mesh_getEdge2no.o \
 mesh_getFixed.o  \
 mesh_load.o  \
 mesh_print.o  \
 mesh_refine.o  \
 mesh_util.o  \
 sed_buildS_nonlin.o \
 sed_buildS.o \
 sed_cholsol_res.o \
 sed_dupl.o \
 sed_gs_constr.o \
 sed_gaxpy.o \
 sed_print.o \
 sed_util.o  \
 sky_cholsol.o  \
 sky_cholesky.o  \
 sky_print.o  \
 sky_sed.o  \
 sky_util.o  \
 ../Csparse/Lib/*.o

$(HPC): ../Include/hpc.h Makefile

%.o: ../Source/%.c ../Include/hpc.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

