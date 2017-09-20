include ../../../Makefile.def

OBJS       =  CarbonSteelEC3.o ConcreteEC2.o HeatTransferMaterial.o SimpleMaterial.o \
             SteelASCE.o TestMaterial.o TestMaterial2.o SFRMCoating.o StainlessSteelEC.o

all:         $(OBJS)
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/AnisotropicMaterial; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/IsotropicMaterial; $(MAKE);
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/OrthotropicMaterial; $(MAKE);
		

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/AnisotropicMaterial; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/IsotropicMaterial; $(MAKE) wipe;
	@$(CD) $(FE)/HeatTransfer/HeatTransferMaterial/OrthotropicMaterial; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
