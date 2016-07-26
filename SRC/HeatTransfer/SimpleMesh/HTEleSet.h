
#ifndef HTEleSet_h
#define HTEleSet_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>

class ID;

class HTEleSet: public TaggedObject
{
public:
	HTEleSet(int tag);
	
  ~HTEleSet ();
	
  int addEleID(ID& newEleID);
	const ID& getEleID(void);
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
  ID EleID;
};
#endif