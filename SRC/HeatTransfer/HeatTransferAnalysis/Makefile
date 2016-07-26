include ../../../Makefile.def

OBJS       = HeatTransferAnalysis.o  HT_TransientAnalysis.o 

all:         $(OBJS)
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/ConvergenceTest; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferAnalysisModel; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferDOFNumber; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferIntegrator; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferSolutionAlgorithm; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/Iterator; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/TemperatureBCHandler; $(MAKE);
	
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/ConvergenceTest; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferAnalysisModel; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferDOFNumber; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferIntegrator; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/HeatTransferSolutionAlgorithm; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/Iterator; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferAnalysis/TemperatureBCHandler; $(MAKE) wipe;
	

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
