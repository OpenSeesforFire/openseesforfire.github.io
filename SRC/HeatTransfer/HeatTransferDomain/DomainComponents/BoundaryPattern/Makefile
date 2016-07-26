include ../../../../../Makefile.def

OBJS       = BoundaryPattern.o  FireImposedPattern.o HeatFluxPattern.o TemperaturePattern.o TimeSeriesPattern.o

all:         $(OBJS)
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
