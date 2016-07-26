include ../../../Makefile.def

OBJS       = SIFColumn.o SIFMember.o SIFSlab.o \
SIFXBeam.o SIFXWall.o SIFYBeam.o  SIFYWall.o SIFXBeamSec.o\

# Compilation control
all:         $(OBJS)
	@$(CD) $(FE)/SIFBuilder/SIFMembers/SIFMemberIter; $(MAKE);


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
