include ../../Makefile.def

OBJS       = FireModel.o  UserDefinedFire.o 

all:         $(OBJS)
	@$(CD) $(FE)/fire/CFD_Interface; $(MAKE);
	@$(CD) $(FE)/fire/EmpericalCompartmentFire; $(MAKE);
	@$(CD) $(FE)/fire/EmpericalLocalizedFire; $(MAKE);
	@$(CD) $(FE)/fire/NorminalFire; $(MAKE);
	@$(CD) $(FE)/fire/ZoneModel; $(MAKE);
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/fire/CFD_Interface; $(MAKE) wipe;
	@$(CD) $(FE)/fire/EmpericalCompartmentFire; $(MAKE) wipe;
	@$(CD) $(FE)/fire/EmpericalLocalizedFire; $(MAKE) wipe;
	@$(CD) $(FE)/fire/NorminalFire; $(MAKE) wipe;
	@$(CD) $(FE)/fire/ZoneModel; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
