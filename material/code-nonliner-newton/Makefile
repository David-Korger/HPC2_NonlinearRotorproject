#------------------------------------------------------------------------------
# CSparse Makefile
#------------------------------------------------------------------------------

C:
	( cd Csparse/Lib ; $(MAKE) )
	( cd Lib ; $(MAKE) )
	( cd Demo ; $(MAKE) )

all: C 

library:
	( cd Csparse/Lib ; $(MAKE) )
	( cd Lib ; $(MAKE) )

clean:
	( cd Csparse/Lib ; $(MAKE) clean )
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )

purge:
	( cd Csparse/Lib ; $(MAKE) purge )
	( cd Lib ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )


