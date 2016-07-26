include ../../Makefile.def

all:
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HTModelBuilder; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HTRecorder; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/SimpleMesh; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/Tcl; $(MAKE);
# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) Makefile.bak $(OBJS) *.o *~ #*# core

spotless: clean
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HTModelBuilder; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HTRecorder; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/SimpleMesh; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/Tcl; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(HEATTRANSFER_LIBRARY);

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
