include ../../../../Makefile.def

OBJS       = HT_AnalysisModel.o HT_DOF_Group.o HT_FE_Element.o Penalty_FE.o HT_PenaltyMP_FE.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
