#ifndef SIFCompartment3D_h
#define SIFCompartment3D_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>

class ID;
class Vector;

class SIFCompartment3D: public TaggedObject
{
public:
	SIFCompartment3D(int tag, double cxx1, double cxx2, double cxx3, double cxx4, double rbxx1, double rbxx2, double rbxx3, double rbxx4, double rsx1, int wll, int wlr, double wlf, double wlb);
	
	//Added by Liming 
	const ID& getConnectedXBeams(void);
	const ID& getConnectedYBeams(void);
	const ID& getConnectedColumns(void);
	const ID& getConnectedSlabs(void);

	
	virtual ~SIFCompartment3D ();
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
    Vector* Mem;
};
#endif