
ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Csparse/Include -I../Include 

LDLIBS += -lm
CS = $(LDFLAGS) ../Csparse/Lib/libcsparse.a ../Lib/libhpc.a $(LDLIBS)

all: lib demo
	- ./demo motor 0
	- ./demo motor 2
	- ./demo motor 4

lib:
	-( cd ../Csparse ; $(MAKE) )
	-( cd ../Lib ; $(MAKE) )

demo: lib demo.c Makefile
	$(CC) $(CF) $(I) -o demo demo.c $(CS)


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r demo \
    *.a *.dSYM *.obj *.dll
