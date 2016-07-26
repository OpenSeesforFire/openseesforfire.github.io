
#ifndef SIFfireAction_h
#define SIFfireAction_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <FireModel.h>
#include <HeatTransferDomain.h>
#include <SIFMember.h>
#include <SIFHTforMember.h>

class ID;
class Vector;
class Matrix;
class SIFBuilderDomain;
class FireModel;
class HeatTransferDomain;
class TaggedObjectStorage;
class TaggedObjectIter;

class SIFfireAction: public TaggedObject
{
public:
	SIFfireAction(int tag, int FireModelType, int CompartmentID);
	
	
	~SIFfireAction();
	
	int SetFireOrigin(const Vector& fireOrigin);
  int SetFireHRR(double );
  int SetFireDia(double );
  

	void setSIFDomain(SIFBuilderDomain *theSifDomain);
	
	int SetFirePath(const Matrix& FirePath);
	int SetFireDuration(double fireDuration);

	int SetStartTime(double startTime);
	
	int UpdateFireModel(int MemberType);
	
	int Apply( int LoadPatternTag, double timeStep=30,double fireDuration=0);
	int RunHTforMember(SIFMember* theMember,int LoadPatternTag);
	
	virtual void  Print(OPS_Stream&, int = 0) {return;};
   
private:
  int FireModelType;
	int CompartmentID;
  double FireHRR;
  double FireDia;
  double StartTime;
	double FireDuration;
	double TimeStep;
	Vector* FireOrigin;
	Matrix* FirePath;
	FireModel* theFireModel;
	SIFBuilderDomain* theSIFDomain;
    HeatTransferDomain* theHTDomain;
		
};
#endif
