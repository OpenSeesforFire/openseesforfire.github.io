#ifndef SIFColumn_h
#define SIFColumn_h

#include <SIFMember.h>

class ID;
class Vector;

class SIFColumn: public SIFMember
{
public:
	SIFColumn(int tag, int typeTag, int jt1, int jt2,const ID& theMemInfo, double gamma=0);
  
	
  const char* getClassType(void) const {return "SIFColumn";};

  int getConnectedSlab(void){return 0;};
  
  double getGamma(void){return 0;};

	
	virtual ~SIFColumn ();
    //virtual int GetEntityTypeTag();
	//int cTag;
	
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
  //double Gamma;

};
#endif