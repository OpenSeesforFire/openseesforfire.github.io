include ../../../../Makefile.def

OBJS       = HeatTransferDomain.o  HeatTransferDomainComponent.o 

all:         $(OBJS)
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/BoundaryPattern; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatFluxBC; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatTransferElement; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatTransferNode; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/MP_TemperatureBC; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/TemperatureBC; $(MAKE);
	

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/BoundaryPattern; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatFluxBC; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatTransferElement; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/HeatTransferNode; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/MP_TemperatureBC; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferDomain/DomainComponents/TemperatureBC; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
