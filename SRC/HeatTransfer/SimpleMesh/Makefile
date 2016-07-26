include ../../../Makefile.def

OBJS       =  Simple_Entity.o Simple_Line.o Simple_Block.o \
Simple_Brick.o Simple_Isection.o  Simple_Isection3D.o Simple_IsecProtected.o\
Simple_Boundary.o Simple_Mesh.o Simple_Composite3D.o \
Simple_Composite2D.o  HTEleSet.o  HTConstants.o  \
Simple_MeshIter.o HTNodeSet.o



all:         $(OBJS)
	
# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean


wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

