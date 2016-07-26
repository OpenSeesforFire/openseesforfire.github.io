#ifndef SIFCompartment_h
#define SIFCompartment_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>
#include <SIFfireAction.h>

class ID;
class Vector;
class SIFfireAction;

class SIFCompartment: public TaggedObject
{
public:
	SIFCompartment(int tag, const ID& Compinfo, const Vector& CompOrigin);
	
	//Added by Liming 
	int AddXBeam(int MemberTag);
	const ID& getConnectedXBeams(void);
	int AddYBeam(int MemberTag);
	const ID& getConnectedYBeams(void);
	
	int AddXBeamSec(int MemberTag);
	const ID& getConnectedSecXBeams(void);
	//int AddYBeamSec(int MemberTag);
	//const ID& getConnectedSecYBeams(void);

	int AddColumn(int MemberTag);
	const ID& getConnectedColumns(void);
	int AddSlab(int MemberTag);
	const ID& getConnectedSlabs(void);
	
	const ID& getCompartmentInfo(void);

	const Vector& getCompartmentOrigin();

	virtual ~SIFCompartment ();
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
    ID theCompInfo;
	Vector theCompOrigin;
	
	
	ID* theXBeamTags;
	ID* theYBeamTags;
	ID* theSecXBeamTags;
	ID* theSecYBeamTags;
	ID* theColumnTags;
	ID* theSlabTags;
	
};
#endif