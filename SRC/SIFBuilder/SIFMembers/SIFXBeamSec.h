#ifndef SIFXBeamSec_h
#define SIFXBeamSec_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFXBeamSec: public SIFMember
{
public:
	SIFXBeamSec(int tag, int typeTag,  int jt1, int jt2, const ID& theMemInfo, double gamma=0);
  
	
  const char* getClassType(void) const {return "SIFXBeamSec";};
   int setConnectedSlab(int slabID);
  int getConnectedSlab(void);
  
  double getGamma(void);
    void addLoad(const Vector& theLoadVec);
	const Vector& getLoad();

	~SIFXBeamSec ();
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	void  Print(OPS_Stream&, int = 0) {return;};

private:
	Vector theLoad;
  double Gamma;
  int ConnectedSlab;
  int SecBeamTag;
  
};
#endif