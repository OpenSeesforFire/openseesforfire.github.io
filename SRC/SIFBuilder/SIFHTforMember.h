
#ifndef SIFHTforMember_h
#define SIFHTforMember_h

#include <TaggedObject.h>
#include <MovableObject.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <FireModel.h>
#include <SIFMember.h>
#include <SIFBuilderDomain.h>
#include <PathTimeSeriesThermal.h>
#include <FireModel.h>
#include <TaggedObject.h>
#include <SIFMaterial.h>
#include <Simple_Mesh.h>



class SIFBuilderDomain;
class ID;
class SIFMember;
class HT_AnalysisModel;
class CTestNormTempIncr;
class HT_SolutionAlgorithm;
class HT_TransientIntegrator;
class TemperatureBCHandler;
class RCM;
class HT_DOF_Numberer;
class BandGenLinLapackSolver;
class BandGenLinSOE;
class HT_TransientAnalysis;

class SIFHTforMember: public TaggedObject
{
public:
	SIFHTforMember(int tag, SIFMember* theMember,SIFBuilderDomain* theSifDomain, double fireDuration, double timeStep );
	
	~SIFHTforMember();

	int SetHTDomain(HeatTransferDomain* theHTdomain);

	int SetFireModel(FireModel* thefireModel);

	int SetFireExposureCons(FireModel* fireModel,ID* theFireExposedfaces);

	int SetAmbExposureCons(const ID& theAmbExposedfaces);
  
    int applyFire(const Vector& fireOrigin =0, Matrix* crdMat=0);

	int BuildHTModel2D(const Vector& SectionLocs, int SecTag=0 );
	int BuildHTModel1D(const Vector& SectionLocs, int SecTag=0 );

	int RunHTAnalysis(double fireDuration, double timeStep);
    int SetPathTimeSeries(PathTimeSeriesThermal** thePathTimeSeries);
	PathTimeSeriesThermal* getHTResults(void);
    
	const Vector& getRecLocations();

	virtual void  Print(OPS_Stream&, int = 0) {return;};

private:
	int Tag;
  double FireDuration, TimeStep;
  double theSlabT;
 
  SIFBuilderDomain* theSIFDomain;
  HeatTransferDomain* theHTDomain;
  
  //SIF Member, Section, Materials
  SIFMember* theMember;
  SIFSection* theSection;
  SIFSection* theSection1;
  
  Simple_Entity* theEntity;
  Simple_Entity* theEntity1;
  
  Simple_Mesh* theMesh;
  Simple_Mesh* theMesh1;

  HeatTransferMaterial* theHTMaterial;
  HeatTransferMaterial* theHTMaterial1;
  
  FireModel** theFireModel;
  ID** theFireExpFacesID;

  ID* theAmbExpFacesID;

  Vector RecLocations;
  PathTimeSeriesThermal** thePathSeries;
  int FireActionIndex;

};
#endif