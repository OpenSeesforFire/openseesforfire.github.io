#ifndef SIFXWall_h
#define SIFXWall_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFXWall: public SIFMember
{
public:
  SIFXWall(int tag, int typeTag, int jt1, int jt2, int jt3, int jt4,const ID& theMemInfo);
  
  
  const char* getClassType(void) const {return "SIFXWall";};
  
  SIFMember* getConnectedBeam(void);
  
  
  
  virtual ~SIFXWall ();
  //virtual int GetEntityTypeTag();
  //int cTag;
  
  
  virtual void  Print(OPS_Stream&, int = 0) {return;};
  
private:
 
  SIFXWall* ConnectedSlab;
  
};
#endif