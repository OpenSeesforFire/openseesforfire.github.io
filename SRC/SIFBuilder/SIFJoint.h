
#ifndef SIFJoint_h
#define SIFJoint_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>

class ID;
class Vector;

class SIFJoint: public TaggedObject
{
public:
	SIFJoint(int tag, double crd1, double crd2, double crd3);
	
	
	~SIFJoint ();
    //virtual int GetEntityTypeTag();
	//int cTag;

	// public methods for obtaining SIFJoints coordinates 
    const Vector &getCrds(void) const;
	int setNodeTag(const ID& nodeTag);
	void setBCtype(int BC_Type);
	int getBCtype();

	void addLoad(const Vector& theLoadVec);
	const Vector& getLoad();

	const ID& getNodeTag(void);

	void  Print(OPS_Stream&, int = 0) {return;};

private:
    Vector* Crd;
	ID NodeTag;      //For nodes meshed at this sifjoint
	int BCtypeTag;
	Vector theLoad;  //For nodal load;
};
#endif