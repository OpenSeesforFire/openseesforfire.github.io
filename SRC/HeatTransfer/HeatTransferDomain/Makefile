include ../../../Makefile.def

OBJS       = 

all:         $(OBJS)
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/Iterator; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/TimeSeries; $(MAKE);
		

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/Iterator; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/TimeSeries; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
