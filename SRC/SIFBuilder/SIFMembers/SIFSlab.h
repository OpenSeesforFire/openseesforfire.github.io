#ifndef SIFSlab_h
#define SIFSlab_h

#include <SIFMember.h>

class ID;
class Vector;
class SIFMember;

class SIFSlab: public SIFMember
{
public:
  SIFSlab(int tag, int typeTag, int jt1, int jt2, int jt3, int jt4,const ID& theMemInfo);
  SIFSlab(int tag, int typeTag, int jt1, int jt2, const ID& theMemInfo);
  
  
  const char* getClassType(void) const {return "SIFSlab";};
  
  //SIFMember* getConnectedBeam(void);
  
   ~SIFSlab ();
  //virtual int GetEntityTypeTag();
  //int cTag;
  
   void addLoad(const Vector& theLoadVec);
	const Vector& getLoad();

  int getConnectedSlab(void){return 0;};

  double getGamma(void){return 0;};
  void  Print(OPS_Stream&, int = 0) {return;};
  
private:
  double Gamma;
  Vector theLoad;
  //SIFSlab* ConnectedSlab;
  
};
#endif