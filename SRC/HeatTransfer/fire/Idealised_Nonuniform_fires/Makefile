include ../../../Makefile.def

OBJS       = AlpertCeilingJetModel.o LocalizedFireEC1.o  Idealised_Local_Fire.o LocalizedFireSFPE.o MovingLocalizedFireEC1.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
