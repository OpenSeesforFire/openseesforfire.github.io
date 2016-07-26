
#ifndef HTNodeSet_h
#define HTNodeSet_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>

class ID;

class HTNodeSet: public TaggedObject
{
public:
	HTNodeSet(int tag);
	
  ~HTNodeSet ();
	
  int addNodeID(ID& newNodeID);
	const ID& getNodeID(void);
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
  ID NodeID;
};
#endif