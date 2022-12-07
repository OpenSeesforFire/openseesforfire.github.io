include ../../Makefile.def

OBJS       = FireModel.o  UserDefinedFire.o 

all:         $(OBJS)
	@$(CD) $(FE)/fire/CFD_Interface; $(MAKE);
	@$(CD) $(FE)/fire/Idealised_Uniform_Fires; $(MAKE);
	@$(CD) $(FE)/fire/Idealised_Nonuniform_fires; $(MAKE);
	@$(CD) $(FE)/fire/ZoneModel; $(MAKE);
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/fire/CFD_Interface; $(MAKE) wipe;
	@$(CD) $(FE)/fire/Idealised_Uniform_Fires; $(MAKE) wipe;
	@$(CD) $(FE)/fire/Idealised_Nonuniform_fires; $(MAKE) wipe;
	@$(CD) $(FE)/fire/ZoneModel; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
