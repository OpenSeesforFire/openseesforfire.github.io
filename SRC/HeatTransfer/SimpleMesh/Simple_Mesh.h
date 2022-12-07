#ifndef SimpleMesh_h
#define SimpleMesh_h


#include <TaggedObject.h>
#include <MovableObject.h>

#include <Simple_Entity.h>
#include <ID.h>
#include <vector>
#include <Simple_Line.h>
#include <Simple_Block.h>
#include <Simple_Isection.h>
#include <Simple_Brick.h>
#include <Simple_Isection3D.h>
#include <Simple_Composite3D.h>
#include <Simple_Composite2D.h>
#include <Simple_IsecProtected.h>
#include <Matrix.h>
#include <Vector.h>
#include <stdlib.h>

//#include <Domain.h>
#include <HeatTransferDomain.h>
//#include <Node.h>
#include <HeatTransferNode.h>
#include <NDMaterial.h>
#include <HeatTransferMaterial.h>
#include <QuadFour.h>
#include <BrickEight.h>
#include <LineTwo.h>
//#include <FourNodeQuad.h>
#include <OPS_Globals.h>
#include <TemperatureBC.h>
#include <BoundaryPattern.h>

class HeatTransferMaterial;
using std::vector;

//typedef vector<double> VECTOR;
	
class Simple_Mesh: public TaggedObject
{

public:
	//Simple_Mesh(int tag, Simple_Block* Block, Domain* theDomain, NDMaterial* theNDMaterial,int NumCtrX, int NumCtrY);
	//Simple_Mesh(int tag, Simple_Isection* Isection, Domain* theDomain,NDMaterial* theNDMaterial, double EleX, double EleY);
	//Simple_Mesh(int tag, Simple_Isection* Isection , Domain* theDomain,NDMaterial* theNDMaterial, double EleX, double EleY,double EleX_Web,double EleY_Web);
	
	Simple_Mesh(int tag, Simple_Entity* Entity, HeatTransferDomain* theDomain, HeatTransferMaterial* theHTMaterial,Vector& MeshCtrls,HeatTransferMaterial* theHTMateria1=0,bool numCtrl=false);
	//Simple_Mesh(int tag, Simple_Brick* Brick, HeatTransferDomain* theDomain, HeatTransferMaterial* theHTMaterial,double EleX, double EleY, double EleZ);
	//Simple_Mesh(int tag, Simple_Isection* Isection, HeatTransferDomain* theDomain, HeatTransferMaterial* theHTMaterial, double EleX, double EleY);
	//Simple_Mesh(int tag, Simple_Isection* Isection , HeatTransferDomain* theDomain, HeatTransferMaterial* theHTMaterial, double EleX, double EleY,double EleX_Web,double EleY_Web);
	//Simple_Mesh(int tag, Simple_Isection3D* Isection3D, HeatTransferDomain* theDomain, HeatTransferMaterial* theHTMaterial, double EleX, double EleY,double EleX_Web,double EleY_Web,double EleZ);
	Simple_Mesh(int tag=0);
	~Simple_Mesh();
    int GetNumOfNodes();
	int GetNumOfEles();
	int SetOriginLocs(const Vector& originLocs);
	int SetEleParameters(const ID& eleParameters);
  
	int GeneratingNodes(const Vector& originLocs=0);
	int GeneratingEles(const ID& EleParameters=0);
	
    int SelectingNodes(ID& NodesRange, int crdTag,double MinValue, double MaxValue, double Tolerance=0.00001 );
	int SelectingNodesbyFace(ID& NodesRange, int FaceTag=1); 
	int SelectingEles(ID&ElesRange, const ID& NodesRange, int eleFaceTag=0);
	int SelectingElesbyFace(ID& ElesRange, int FaceTag, int &eleFaceID);
	
	
    int GetNodesForRecorder(ID& NodesRange, int dimTag, double PlaneLoc=0.0) ;
	int GetNodesForRecorder(ID& NodesRange, const Vector& MonitorLocX, double PlaneLoc=0.0);
	int GetNodesForRecorder(ID& NodesRange, double xLoc=0.0, double yLoc=0.0);
	
	//virtual int getMeshTag();
    const ID& getNumCtrlID();
	void  Print(OPS_Stream&, int = 0) {return;};

private:	
    NDMaterial* theNDMaterial;
	HeatTransferMaterial* theHTMaterial;
	HeatTransferMaterial* theHTMaterial1;
	int OriginNodeTag,OriginEleTag;
	Simple_Entity* theEntity;
	bool isHTDomain;
  
  Vector OriginLocs;
  ID EleParameters;
	
	//Domain* theDomain;
	HeatTransferDomain* theHTDomain;
	


};
#endif
