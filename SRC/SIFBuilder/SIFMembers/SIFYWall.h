#ifndef SIFYWall_h
#define SIFYWall_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFYWall: public SIFMember
{
public:
  SIFYWall(int tag, int typeTag, int jt1, int jt2, int jt3, int jt4,const ID& theMemInfo);
  
  
  const char* getClassType(void) const {return "SIFYWall";};
  
  SIFMember* getConnectedBeam(void);
  
  
  
  virtual ~SIFYWall ();
  //virtual int GetEntityTypeTag();
  //int cTag;
  
  
  virtual void  Print(OPS_Stream&, int = 0) {return;};
  
private:
  
  SIFYWall* ConnectedSlab;
  
};
#endif