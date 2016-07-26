#ifndef SIFYBeam_h
#define SIFYBeam_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFYBeam: public SIFMember
{
public:
  SIFYBeam(int tag, int typeTag, int jt1, int jt2,const ID& theMemInfo, double gamma=0);
  
  
  const char* getClassType(void) const {return "SIFYBeam";};
  
  int setConnectedSlab(int slabID);
  int getConnectedSlab(void);
  
  double getGamma(void);
   void addLoad(const Vector& theLoadVec);
	const Vector& getLoad();

  
  
  virtual ~SIFYBeam ();
  //virtual int GetEntityTypeTag();
  //int cTag;
  
  
  virtual void  Print(OPS_Stream&, int = 0) {return;};
  
private:
  Vector theLoad;
  double Gamma;
   int ConnectedSlab;
  
};
#endif