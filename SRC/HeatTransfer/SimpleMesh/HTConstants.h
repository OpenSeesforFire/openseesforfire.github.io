
#ifndef HTConstants_h
#define HTConstants_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>

class ID;
class Vector;

class HTConstants: public TaggedObject
{
public:
	HTConstants(int tag, Vector& constants);
	
  ~HTConstants ();
	
	const Vector& getConstants(void);
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
  Vector Constants;
};
#endif