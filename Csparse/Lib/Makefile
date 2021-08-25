# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
#

# --- System configuration file for Linux with gcc

LIBRARY = libcsparse

CC = gcc
ARCHIVE = ar
RANLIB = ranlib

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3 -fPIC

I = -I../Include

AR_TARGET = $(LIBRARY).a

all: install

#-------------------------------------

CS = cs_add.o cs_amd.o cs_chol.o cs_cholsol.o cs_counts.o cs_cumsum.o \
	cs_droptol.o cs_dropzeros.o cs_dupl.o cs_entry.o \
	cs_etree.o cs_fkeep.o cs_gaxpy.o cs_happly.o cs_house.o cs_ipvec.o \
	cs_lsolve.o cs_ltsolve.o cs_lu.o cs_lusol.o cs_util.o cs_multiply.o \
	cs_permute.o cs_pinv.o cs_post.o cs_pvec.o cs_qr.o cs_qrsol.o \
	cs_scatter.o cs_schol.o cs_sqr.o cs_symperm.o cs_tdfs.o cs_malloc.o \
	cs_transpose.o cs_compress.o cs_usolve.o cs_utsolve.o cs_scc.o \
	cs_maxtrans.o cs_dmperm.o cs_updown.o cs_print.o cs_norm.o cs_load.o \
	cs_dfs.o cs_reach.o cs_spsolve.o cs_ereach.o cs_leaf.o cs_randperm.o

$(CS): ../Include/cs.h Makefile

%.o: ../Source/%.c ../Include/cs.h
	$(CC) $(CF) $(I) -c $<

#-------------------------------------

$(AR_TARGET): $(CS)
	$(ARCHIVE) r $@ $^
	- $(RANLIB) $@

# install CSparse in this directory
install: $(AR_TARGET)

#-------------------------------------

purge: distclean

clean:
	- $(RM) *.o

distclean: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

