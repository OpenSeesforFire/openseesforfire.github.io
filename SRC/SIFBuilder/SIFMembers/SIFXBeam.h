#ifndef SIFXBeam_h
#define SIFXBeam_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFXBeam: public SIFMember
{
public:
	SIFXBeam(int tag, int typeTag, int jt1, int jt2, const ID& theMemInfo, double gamma=0);
  
	
  const char* getClassType(void) const {return "SIFXBeam";};
   int setConnectedSlab(int slabID);
  int getConnectedSlab(void);
  
  double getGamma(void);
    void addLoad(const Vector& theLoadVec);
	const Vector& getLoad();

	
	~SIFXBeam ();
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	void  Print(OPS_Stream&, int = 0) {return;};

private:
	Vector theLoad;
  double Gamma;
  int ConnectedSlab;
  
};
#endif