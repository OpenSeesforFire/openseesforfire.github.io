#ifndef SimpleBoundary_h
#define SimpleBoundary_h


#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>
#include <stdlib.h>

#include <Domain.h>
#include <HeatTransferDomain.h>
#include <Node.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <FourNodeQuad.h>
#include <HeatTransferMaterial.h>
#include <OPS_Globals.h>

#include <SP_Constraint.h>
#include <TemperatureBC.h>
#include <MP_Constraint.h>
#include <MP_TemperatureBC.h>
#include <BoundaryPattern.h>
#include <HeatFluxBC.h>
#include <Convection.h>
#include <Radiation.h>
#include <PrescribedSurfFlux.h>

class HeatTransferDomain;
class HeatTransferNode;
class QuadFour;
class HeatTransferMaterial;
using std::vector;

//typedef vector<double> VECTOR;
	
class Simple_Boundary
{
public:
	Simple_Boundary(int tag, HeatTransferDomain* theDomain);
	Simple_Boundary(int tag, Domain* theDomain);

	~Simple_Boundary();
   


	virtual void GeneratingSP_BC(const ID& NodesRange,int nDoF, double value, bool ISconstant,int LoadPatternTag=-1);
	virtual void GeneratingMP_BC(const ID& NodesRange,const ID& NodesRange1,int DoFTag0,int DoFTag1=-1);

	virtual void GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType=0);
	virtual int getTag();
   

private:	
	int BoundaryTag;
	Domain* theDomain;
	HeatTransferDomain* theHTDomain;
	bool isHTDomain;

};
#endif
