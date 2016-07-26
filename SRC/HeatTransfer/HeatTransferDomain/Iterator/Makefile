include ../../../../Makefile.def

OBJS       = AllTempBCIter.o  BoundaryPatternIter.o HeatFluxBCIter.o HT_ElementIter.o HT_NodeIter.o HTDomain_Iter.o \
			MP_TempBCIter.o TempBCIter.o

all:         $(OBJS)
	
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
