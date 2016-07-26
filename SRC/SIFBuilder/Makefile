include ../../Makefile.def

OBJS       = SIFBuilderDomain.o SIFCompartment.o SIFCompartmentIter.o \
			SIFfireAction.o SIFfireActionIter.o \
			SIFHTforMember.o SIFJoint.o SIFJointIter.o SIFMaterial.o  SIFSection.o \

# Compilation control
all:         $(OBJS)
	@$(CD) $(FE)/SIFBuilder/tcl; $(MAKE);
	@$(CD) $(FE)/SIFBuilder/SIFMembers; $(MAKE);


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
