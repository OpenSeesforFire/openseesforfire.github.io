
#ifndef SIFMaterial_h
#define SIFMaterial_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>
#include <ID.h>
#include <HeatTransferMaterial.h>
#include <UniaxialMaterial.h>
class ID;

class SIFMaterial: public TaggedObject
{
public:
	SIFMaterial(int tag, int MaterialTypeTag, 
	double fy, double E0);

	SIFMaterial(int tag, int MaterialTypeTag,
		double fy, double fu, double E0);
	
	SIFMaterial(int tag, int MaterialTypeTag, 
	double fc, double epsc0, double fcu, double epscu, double rat, double ft,double Ets, double moisture);
	
	
	~SIFMaterial();
	
	int getMaterialTypeTag(void);
	const Vector& getMaterialPars(void);
  
  //define HeatTransfer Material for HT analysis;
  HeatTransferMaterial* getHeatTransferMaterial( );
  
  UniaxialMaterial* getUniaxialMaterial(bool isElastic = false);

  double getInitialModulus();
  
	virtual void  Print(OPS_Stream&, int = 0) {return;};
   
private:
  int MaterialTypeTag;
	Vector MaterialPars;
  HeatTransferMaterial* theHTMaterial;
	
};
#endif