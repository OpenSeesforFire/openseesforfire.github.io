/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                       


/**********************************************************************
** This project is aiming to provide a Tcl interface to define models  **
** for simulating structural behaviours under fire action.           **
** Developed by:                  `                                   **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $


#include <SIFBuilderDomain.h>

#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>
#include <RigidRod.h>
#include <RigidBeam.h>
//Iterator for MapofTaggedObjects


#include <BeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>

#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <DispBeamColumn3d.h>
#include <DispBeamColumn3dThermal.h>
//#include <DispBeamColumn3dThermal.h>
#include <ShellMITC4Thermal.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <NodalLoad.h>
#include <Beam3dUniformLoad.h>
// analysis model
#include <AnalysisModel.h>
#include <RCM.h>
// convergence tests
#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>
#include <CTestRelativeNormUnbalance.h>
#include <CTestRelativeNormDispIncr.h>
#include <CTestRelativeEnergyIncr.h>
#include <CTestRelativeTotalNormDispIncr.h>
#include <CTestFixedNumIter.h>
#include <NormDispAndUnbalance.h>
#include <NormDispOrUnbalance.h>

// soln algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <NewtonLineSearch.h>
#include <ModifiedNewton.h>

// constraint handlers
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
//#include <PenaltyHandlerNoHomoSPMultipliers.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <DOF_Numberer.h>
#include <PlainNumberer.h>
#include <DOF_Numberer.h>

// integrators
#include <LoadControl.h>
#include <ArcLength.h>

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

#include <SymSparseLinSolver.h>
#include <SymSparseLinSOE.h>

#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ElasticIsotropic3DThermal.h>
#include <ElasticIsotropicMaterial.h>
#include <PlateFiberMaterialThermal.h>
#include <MembranePlateFiberSectionThermal.h>
#include <SuperLU.h>
#include <SparseGenColLinSOE.h>

static AnalysisModel *theAnalysisModel =0;
static EquiSolnAlgo *theAlgorithm =0;
static ConstraintHandler *theHandler =0;
static DOF_Numberer *theNumberer =0;
static LinearSOE *theSOE =0;
static StaticIntegrator *theStaticIntegrator =0;
static ConvergenceTest *theTest =0;
static StaticAnalysis *theStaticAnalysis = 0;
static DirectIntegrationAnalysis *theTransientAnalysis = 0;
static VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis = 0;

SIFBuilderDomain::SIFBuilderDomain()
{
	
	theSIFMaterials = new ArrayOfTaggedObjects(32);
	theSIFSections = new ArrayOfTaggedObjects(32);
	//theSIFMembers = new ArrayOfTaggedObjects(500);
	
    
  theSIFJoints = new MapOfTaggedObjects();
  theSIFJoint_Iter  = new SIFJointIter(theSIFJoints);

  theSIFCompartments = new MapOfTaggedObjects();
  theSIFCompartment_Iter= new SIFCompartmentIter(theSIFCompartments);

  theSIFSecXBeams = new MapOfTaggedObjects();
  theSIFSecXBeam_Iter  = new SIFSecXBeamIter(theSIFSecXBeams);

  theSIFXBeams = new MapOfTaggedObjects();
  theSIFXBeam_Iter  = new SIFXBeamIter(theSIFXBeams);
  theSIFYBeams = new MapOfTaggedObjects();
  theSIFYBeam_Iter  = new SIFYBeamIter(theSIFYBeams);
  theSIFColumns = new MapOfTaggedObjects();
  theSIFColumn_Iter  = new SIFColumnIter(theSIFColumns);
  theSIFSlabs = new MapOfTaggedObjects();
  theSIFSlab_Iter  = new SIFSlabIter(theSIFSlabs);

  theSIFXWalls = new ArrayOfTaggedObjects(1000);
  theSIFYWalls = new ArrayOfTaggedObjects(1000);
  
  theSIFfireActions = new MapOfTaggedObjects();
  theSIFfireAction_Iter  = new SIFfireActionIter(theSIFfireActions);
	
  //check whether the space is properly allocated
  if (theSIFMaterials == 0 || theSIFSections == 0 || theSIFXBeams == 0 || theSIFXWalls == 0
      || theSIFCompartments == 0 || theSIFfireActions == 0
      || theSIFfireAction_Iter == 0 )
  {
    opserr << "HeatTransferDomain::HeatTransferDomain() - out of memory\n";
    exit(-1);
  }

  isElastic = false;
  LoadApplied = 0;
  
  theNodeTag = 0;                     //static tag for node; 
  theEleTag = 0;

  InterpPWD = 0;
  theHTdir = 0;
  SIFBuilderInfo = 0;
  offset =0;

  fireFloorTag = 0;

  	theNodalLoadTag=0;
	theEleLoadTag=0;
	theDomain=0;
	nDoF =6;
}


SIFBuilderDomain::~SIFBuilderDomain()
{
  this->clearAll();
  if(theSIFMaterials != 0)
	  delete theSIFMaterials;

    if(theSIFSections != 0)
	  delete theSIFSections;
	
	if(theSIFCompartments != 0)
	  delete theSIFCompartments;

	if(theSIFJoints != 0)
	  delete theSIFJoints;

	if(theSIFXBeams != 0)
	  delete theSIFXBeams;

	if(theSIFXBeam_Iter != 0)
	  delete theSIFXBeam_Iter;

	if(theSIFYBeams != 0)
	  delete theSIFYBeams;

	if(theSIFYBeam_Iter != 0)
	  delete theSIFYBeam_Iter;

	
	if(theSIFColumns != 0)
	  delete theSIFColumns;
	
	if(theSIFColumn_Iter != 0)
	  delete theSIFColumn_Iter;
	
	if(theSIFSlabs != 0)
	  delete theSIFSlabs;
	
	if(theSIFSlab_Iter != 0)
	  delete theSIFSlab_Iter;

	if(theSIFXWalls != 0)
	  delete theSIFXWalls;
	
	if(theSIFfireActions != 0)
	  delete theSIFfireActions;

	if(theSIFfireAction_Iter != 0)
	  delete theSIFfireAction_Iter;

	if(theNodeTag !=0||theEleTag !=0){
		theNodeTag =0;
		theEleTag = 0;
	}



}
int 
SIFBuilderDomain::SetStructureDomain(Domain* theStructDomain)
{
		
	theDomain= theStructDomain;
		
	return 0;

}

Domain* 
SIFBuilderDomain::getStructureDomain()
{
	return theDomain;
}

void 
SIFBuilderDomain::setInterpPWD(const char *pwd)
{
	InterpPWD = pwd;
}

const char* 
SIFBuilderDomain::getInterpPWD()
{
	return InterpPWD;
}

void 
SIFBuilderDomain::setHTdir(const char *HTdir)
{
	theHTdir = HTdir;
}

const char* 
SIFBuilderDomain::getHTdir()
{
	return theHTdir;
}

int 
SIFBuilderDomain::setSIFBuilderInfo(const ID& theBuilderInfo)
{
	SIFBuilderInfo=ID();
	SIFBuilderInfo = theBuilderInfo;
	return 0;
}

const ID& 
SIFBuilderDomain::getSIFBuilderInfo()
{

	return SIFBuilderInfo;
}


void
SIFBuilderDomain::clearAll(){

	if(theDomain!=0){
		theDomain->clearAll();
	}
	
	theDomain = 0 ;                 //static pointer to Domain;
	nDoF =6;                            //static tag for number of degree of freedom
	theNodalLoadTag =0;
	theEleLoadTag =0;

	
  if (theStaticAnalysis != 0) {
      theStaticAnalysis->clearAll();//all the analysis facilities are released
      theStaticAnalysis =0;
  }

  theAnalysisModel= 0;
  theStaticIntegrator = 0;
  theTest = 0;
  theAlgorithm=0;
  theHandler=0;
  theNumberer=0;
  theSOE =0;
  theStaticAnalysis = 0;
}

int 
SIFBuilderDomain::addSIFMaterial(SIFMaterial *theSIFMaterial)
{
  bool result = theSIFMaterials->addComponent(theSIFMaterial);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFMaterial() - failed to add SIF Material: " << theSIFMaterial;
    return -1;
  }
}

SIFMaterial *
SIFBuilderDomain::getSIFMaterial(int tag)
{
  TaggedObject *mc = theSIFMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFMaterial *result = (SIFMaterial *)mc;
  return result;
}


int 
SIFBuilderDomain::addSIFSection(SIFSection *theSIFSection)
{
  bool result = theSIFSections->addComponent(theSIFSection);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFSection() - failed to add SIF Material: " << theSIFSection;
    return -1;
  }
}


SIFSection *
SIFBuilderDomain::getSIFSection(int tag)
{
  TaggedObject *mc = theSIFSections->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFSection *result = (SIFSection *)mc;
  return result;
}



int 
SIFBuilderDomain::addSIFJoint(SIFJoint *theSIFJoint)
{
 bool result = theSIFJoints->addComponent(theSIFJoint);
 if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFJoint() - failed to add SIF Joint: " << theSIFJoint;
    return -1;
  }
}


SIFJoint*
SIFBuilderDomain::getSIFJoint(int tag)
{
  TaggedObject *mc = theSIFJoints->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFJoint *result = (SIFJoint *)mc;
  return result;	
}

SIFJointIter &
SIFBuilderDomain::getSIFJoints()
{
  theSIFJoint_Iter->reset();
  return *theSIFJoint_Iter;
}


//-------------------------Sec-XBeam----------------------------//
int 
SIFBuilderDomain::addSIFXBeamSec(SIFXBeamSec *theSIFXBeamSec)
{
  bool result = theSIFSecXBeams->addComponent(theSIFXBeamSec);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFXBeamSec() - failed to add SIFSecXBeam: " << theSIFXBeamSec;
    return -1;
  }
}


SIFXBeamSec*
SIFBuilderDomain::getSIFXBeamSec(int tag)
{
  TaggedObject *mc = theSIFSecXBeams->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFXBeamSec *result = (SIFXBeamSec *)mc;
  return result;
}


SIFSecXBeamIter &
SIFBuilderDomain::getSIFSecXBeams()
{
  theSIFSecXBeam_Iter->reset();
  return *theSIFSecXBeam_Iter;
}

//--------------------------XBeam----------------------------//
int 
SIFBuilderDomain::addSIFXBeam(SIFXBeam *theSIFXBeam)
{
  bool result = theSIFXBeams->addComponent(theSIFXBeam);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFXBeam() - failed to add SIF XBeam: " << theSIFXBeam;
    return -1;
  }
}


SIFXBeam*
SIFBuilderDomain::getSIFXBeam(int tag)
{
  TaggedObject *mc = theSIFXBeams->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFXBeam *result = (SIFXBeam *)mc;
  return result;
}


SIFXBeamIter &
SIFBuilderDomain::getSIFXBeams()
{
  theSIFXBeam_Iter->reset();
  return *theSIFXBeam_Iter;
}


int
SIFBuilderDomain::addSIFYBeam(SIFYBeam *theSIFYBeam)
{
  bool result = theSIFYBeams->addComponent(theSIFYBeam);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFYBeam() - failed to add SIF YBeam: " << theSIFYBeam;
    return -1;
  }
}


SIFYBeam*
SIFBuilderDomain::getSIFYBeam(int tag)
{
  TaggedObject *mc = theSIFYBeams->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  SIFYBeam *result = (SIFYBeam *)mc;
  return result;
}

SIFYBeamIter &
SIFBuilderDomain::getSIFYBeams()
{
  theSIFYBeam_Iter->reset();
  return *theSIFYBeam_Iter;
}


int
SIFBuilderDomain::addSIFColumn(SIFColumn *theSIFColumn)
{
  bool result = theSIFColumns->addComponent(theSIFColumn);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addtheSIFColumn() - failed to add SIF Column: " << theSIFColumn;
    return -1;
  }
}


SIFColumn*
SIFBuilderDomain::getSIFColumn(int tag)
{
  TaggedObject *mc = theSIFColumns->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  SIFColumn *result = (SIFColumn *)mc;
  return result;
}

SIFColumnIter &
SIFBuilderDomain::getSIFColumns()
{
  theSIFColumn_Iter->reset();
  return *theSIFColumn_Iter;
}

int 
SIFBuilderDomain::addSIFXWall(SIFXWall *theSIFXWall)
{
bool result = theSIFXWalls->addComponent(theSIFXWall);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFXWall() - failed to add SIF Wall: " << theSIFXWall;
    return -1;
  }
}

SIFXWall*
SIFBuilderDomain::getSIFXWall(int tag)
{
TaggedObject *mc = theSIFXWalls->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFXWall *result = (SIFXWall *)mc;
  return result;
}

int 
SIFBuilderDomain::addSIFCompartment(SIFCompartment *theSIFCompartment)
{
bool result = theSIFCompartments->addComponent(theSIFCompartment);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFCompartment3D() - failed to add SIF Compartment3D: " << theSIFCompartment;
    return -1;
  }
}


SIFCompartmentIter &
SIFBuilderDomain::getSIFCompartments()
{
  theSIFCompartment_Iter->reset();
  return *theSIFCompartment_Iter;
}

SIFCompartment*
SIFBuilderDomain::getSIFCompartment(int tag)
{
TaggedObject *mc = theSIFCompartments->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFCompartment *result = (SIFCompartment *)mc;
  return result;
}

SIFCompartment*
SIFBuilderDomain::getSIFCompartment(const ID& compInfo)
{
	SIFCompartment *theSIFCompartment;
	if(compInfo.Size()!=3){
		opserr<<"WARNING:: getSIFCompartment recieves an invalid compInfo"<<endln;
		return 0;
		}
    SIFCompartmentIter &theSIFCompartments  = this->getSIFCompartments();
	while((theSIFCompartment = theSIFCompartments())!=0){
		ID theCompInfo = theSIFCompartment->getCompartmentInfo();
	if(theCompInfo(0)==compInfo(0)&&theCompInfo(1)==compInfo(1)&&theCompInfo(2)==compInfo(2)){
			return theSIFCompartment;
		}
	}
}

int 
SIFBuilderDomain::addSIFSlab(SIFSlab *theSIFSlab)
{
  bool result = theSIFSlabs->addComponent(theSIFSlab);
  if (result == true)
    return 0;
  else {
    opserr << "SIFBuilderDomain::addSIFSlab() - failed to add SIF Slab: " << theSIFSlab;
    return -1;
  }
}


SIFSlab*
SIFBuilderDomain::getSIFSlab(int tag)
{
TaggedObject *mc = theSIFSlabs->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFSlab *result = (SIFSlab *)mc;
  return result;	
}


SIFSlabIter &
SIFBuilderDomain::getSIFSlabs()
{
  theSIFSlab_Iter->reset();
  return *theSIFSlab_Iter;
}



int 
SIFBuilderDomain::addSIFfireAction(SIFfireAction *theSIFfireAction)
{
  // first check if a mesh with a similar tag exists in model
  int tag = theSIFfireAction->getTag();
  TaggedObject *other = theSIFfireActions->getComponentPtr(tag);
  if (other != 0) {
    opserr << "SIFBuilderDomain::addSIFfireAction - cannot add as SIFfireACction with tag" <<
    tag << "already exists in model\n";
    
    return false;
  }

 ID CompId = (this->getSIFCompartment(tag))->getCompartmentInfo();
  if (fireFloorTag == 0)
	  fireFloorTag = Vector(10);
  fireFloorTag(0) = CompId(2);

  
  // now we add the load pattern to the container for load pattrens
  bool result = theSIFfireActions->addComponent(theSIFfireAction);
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addSIFfireAction() - failed to add SIFfireAction: " <<theSIFfireAction;
    return -1;
  }


}


SIFfireAction*
SIFBuilderDomain::getSIFfireAction (int tag)
{
TaggedObject *mc = theSIFfireActions->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  SIFfireAction *result = (SIFfireAction *)mc;
  return result;
}


SIFfireActionIter &
SIFBuilderDomain::getSIFfireActions()
{
  theSIFfireAction_Iter->reset();
  return *theSIFfireAction_Iter;
}



////////////--------------Building Strucutral Model-------------------------------//////
int
SIFBuilderDomain::GenStructuralModel(const ID& MeshCtrlPars, bool isElasticModel, bool isDynamicModel, bool isPinnedModel)
{
  isElastic = isElasticModel;
  isPinned = isPinnedModel;
  LoadApplied = 0;
  //Mesh Number control
  int NumElesX; int NumElesY; int NumElesZ;
  NumElesX =MeshCtrlPars(0); NumElesY =MeshCtrlPars(1); NumElesZ = MeshCtrlPars(2);



  //First: To generate nodes on SIFJoints
  SIFJoint *theSIFJoint;
  SIFJointIter &theSIFJoints  = this->getSIFJoints();

  while((theSIFJoint = theSIFJoints())!=0){
    Vector theSIFJointCrds = theSIFJoint->getCrds();
	if(theSIFJointCrds.Size()!=3){
		opserr<<"WARNING::theSIFJoint "<<theSIFJoint->getTag()
			  <<" doesn't have compatible coordinates"<<endln;
	}
	
	//Node::Node(int tag, int ndof, double Crd1, double Crd2, double Crd3)
	//Create a node attached to SIFJoint
	theNodeTag++;
	//double offset;
	Node* theNode = new Node(theNodeTag, nDoF,theSIFJointCrds(0),theSIFJointCrds(1),theSIFJointCrds(2) );
	theDomain->addNode(theNode);
	//-------For dynamic analysis, nodes are created with mass attached.
	

    ID theNodeID = ID(1);
	theNodeID(0)= theNodeTag;
	theSIFJoint->setNodeTag(theNodeID);

  }


  int NumElesPerXBeam = NumElesX; 
  int NumElesPerYBeam = NumElesZ; 
  int NumElesPerColumn =NumElesY;
  //Secon00d: Generate interior nodes in SIFMember and mesh;
  //Here meshCtrl can be easily transformed to element size control;
  SIFXBeam *theSIFXBeam;
  SIFXBeamIter &theSIFXBeams  = this->getSIFXBeams();
  while((theSIFXBeam = theSIFXBeams())!=0){
	this->MeshSIFMember(theSIFXBeam,NumElesPerXBeam);
  }
  
  SIFYBeam *theSIFYBeam;
  SIFYBeamIter &theSIFYBeams  = this->getSIFYBeams();
  while((theSIFYBeam = theSIFYBeams())!=0){
    this->MeshSIFMember(theSIFYBeam,NumElesPerYBeam);
  }
  
  SIFColumn *theSIFColumn;
  SIFColumnIter &theSIFColumns  = this->getSIFColumns();
  while((theSIFColumn = theSIFColumns())!=0){
    this->MeshSIFMember(theSIFColumn,NumElesPerColumn);
  }

SIFXBeamSec *theSIFXBeamSec;
  SIFSecXBeamIter &theSIFSecXBeams  = this->getSIFSecXBeams();
  while((theSIFXBeamSec = theSIFSecXBeams())!=0){
	this->MeshSIFMember(theSIFXBeamSec,NumElesPerXBeam);
  }

 
   SIFSlab *theSIFSlab;
  SIFSlabIter &theSIFSlabs  = this->getSIFSlabs();
  while((theSIFSlab = theSIFSlabs())!=0){
	  if(this->getSIFBuilderInfo()(0)==12)
		this->MeshSIFMember(theSIFSlab,NumElesPerXBeam);
	  else if(this->getSIFBuilderInfo()(0)==13)
		this->MeshSIFMember(theSIFSlab,NumElesPerYBeam);
	  else{
		  if(theSIFSlab->getMemberInfo()(2)==fireFloorTag(0))
			this->MeshSIFSlab(theSIFSlab,NumElesPerXBeam,NumElesPerYBeam);
	  }
  }
  
   opserr<<theDomain->getNumNodes()<<" Nodes, "<<theDomain->getNumElements()<< " elements being created"<<endln;
  
  
  //Third: To generate shell element for the slab
   //------------------------------------SP constraints for frame model------------------------------
   if (this->getSIFBuilderInfo()(0) == 11 || this->getSIFBuilderInfo()(0) == 12) {
	   NodeIter &theNodes = theDomain->getNodes();
	   Node *theNode;

	   int SPdofTag = 0;
	   if (this->getSIFBuilderInfo()(0) == 11)
		   SPdofTag = 0;
	   else if (this->getSIFBuilderInfo()(0) == 12)
		   SPdofTag = 2;

	   while ((theNode = theNodes()) != 0) {
		   SP_Constraint *theSP1 = new SP_Constraint(theNode->getTag(), SPdofTag, 0.0, true);
		   SP_Constraint *theSP2 = new SP_Constraint(theNode->getTag(), 3, 0.0, true);
		   SP_Constraint *theSP3 = new SP_Constraint(theNode->getTag(), 4, 0.0, true);
		   theDomain->addSP_Constraint(theSP1);
		   theDomain->addSP_Constraint(theSP2);
		   theDomain->addSP_Constraint(theSP3);
	   }
   }
  return 0;    
}

////////////////----------------------Define Nodal Masses-----------------------------//////////
	/*
int
SIFBuilderDomain::defineMass(){
	//Matrix mass(6,6);
	//theNode->setMass(mass);

	SIFJoint *theSIFJoint;
	int theJtNodeTag;
	double NodalMass;
	NodalMass = 2250; //Temp Setting;
	SIFJointIter &theSIFJoints  = this->getSIFJoints();

	while((theSIFJoint = theSIFJoints())!=0){
		if(theSIFJoint->getTag()<10000){
			int theJtNodeTag = theSIFJoint->getNodeTag()(0);


		}
	
	}




}
*/



////////////////----------------------Define Nodal Masses-----------------------------//////////





////////////////----------------------Define Boundary conditions-----------------------------//////////
int
SIFBuilderDomain::defineBC(){

//-----------SP_Constraint---------------

	SIFJoint *theSIFJoint;
    SIFJointIter &theSIFJoints  = this->getSIFJoints();
	
	//loop over every joint to define the contraints
	while((theSIFJoint = theSIFJoints())!=0){
	//always choose the first node on this SIFJoint
	//The other nodes are defined for multi-point connection
	int theJointNodeTag = (theSIFJoint->getNodeTag())(0);
	
	if(theSIFJoint->getBCtype()==1){
		//pinned connection
		for(int i =0; i<3; i++){
		   SP_Constraint *theSP = new SP_Constraint(theJointNodeTag, i, 0.0, true);
		   theDomain->addSP_Constraint(theSP);
		}
	}
	else if (theSIFJoint->getBCtype()==2)
	{	
		//fixed connection
		//SP_Constraint(int nodeTag, int ndof, double value, bool isConstant);
		for(int i =0; i<6; i++){
		   SP_Constraint *theSP = new SP_Constraint(theJointNodeTag, i, 0.0, true);
		   theDomain->addSP_Constraint(theSP);
		} 
	}
	}//end of while loop
	//--------end of SP_Constraint------------
 return 0;
}



////////////////----------------------MeshSIFMember--(FramBeamColumn)------------------//////////

int
SIFBuilderDomain::MeshSIFMember(SIFMember* theSIFMember, int NumEles)
{
    //SIFMemberType:   //1:xBeam,2:YBeam,3:Column,4:Slab, 21 xSecBeam, 11: slab modelled as plane frame
	
	int nIP =5;//Number fo integration points
	ID theJoints = theSIFMember->getConnectedJoints();

	int Jt1Tag = theJoints(0);
	int Jt2Tag = theJoints(1);
	
	SIFJoint* Jt1 = this->getSIFJoint (Jt1Tag);
    SIFJoint* Jt2 = this->getSIFJoint (Jt2Tag);

	Vector Jt1Crds = Jt1->getCrds();
	Vector Jt2Crds = Jt2->getCrds();

	double Jt1xCrd = Jt1Crds (0);
	double Jt1yCrd = Jt1Crds (1);
	double Jt1zCrd = Jt1Crds (2);

	double Jt2xCrd = Jt2Crds (0);
	double Jt2yCrd = Jt2Crds (1);
	double Jt2zCrd = Jt2Crds (2);

	//Secondary XBeam
	//if(theSIFMember->getMemberTypeTag() ==21){
		//int theCompTag = (SIFXBeamSec*)theSIFMember->getCompartmentTag();
		//SIFCompartment* theComp = this->getSIFCompartment(theCompTag);
		//theComp->getConnectedXBeams();
	//}


  ID intNodes = ID(NumEles-1);  // for storing the interior nodes
 
  double NdCrdx, NdCrdy, NdCrdz;
  //double offset= 0;

  //determine offset
  if(this->getSIFBuilderInfo()(0)!=2){
	  //Frame structure without slab
  if(theSIFMember->getMemberTypeTag()==1||theSIFMember->getMemberTypeTag()==2||theSIFMember->getMemberTypeTag()==21){
	  if(offset==0){
	  int SlabID = theSIFMember->getConnectedSlab();
	 SIFSlab* theSlab = 0;
	 if(SlabID!=0)
		 theSlab = this->getSIFSlab(SlabID);

	 SIFSection* theSlabSection = theSlab->getSIFSectionPtr();
	 double slabT = theSlabSection->getSectionPars()(0);
	 SIFSection* theBeamSection = theSIFMember->getSIFSectionPtr();
	 Vector SectionPars = theBeamSection->getSectionPars();
	 if(theBeamSection->getSectionTypeTag()==1){
			offset = SectionPars(1)/2+slabT/2;
	 }
	 else if(theBeamSection->getSectionTypeTag()==2){
			offset = SectionPars(0)/2+slabT/2;
	 }
	}
//end of if offset ==0;
   }
  }
  //end of determining offset
 // offset =0;
  ///////////////////////////////Mesh member in nodes///////////////////////////////////
  for (int i = 1; i < NumEles; i++)
  {
    theNodeTag++;
    
    NdCrdx = Jt1xCrd + (Jt2xCrd- Jt1xCrd)/NumEles*i;
	if(theSIFMember->getMemberTypeTag()==11)
		NdCrdy = Jt1yCrd + (Jt2yCrd- Jt1yCrd)/NumEles*i +offset;
	else 
		NdCrdy = Jt1yCrd + (Jt2yCrd- Jt1yCrd)/NumEles*i;

    NdCrdz = Jt1zCrd + (Jt2zCrd- Jt1zCrd)/NumEles*i;
    
    Node *theNode = new Node(theNodeTag, nDoF, NdCrdx, NdCrdy, NdCrdz);
  
    theDomain -> addNode(theNode);
    
    intNodes(i-1) = theNodeTag;
    
  }


// For compiste beam, pinned beam, extra nodes are needed
  int NodeJt1Tag = 0; int NodeJt2Tag = 0;
  NodeJt1Tag = Jt1->getNodeTag()(0);
  NodeJt2Tag = Jt2->getNodeTag()(0);

  //---------------Extra nodes and connections for Secondary Xbeam---------------------------------------//
  if (theSIFMember->getMemberTypeTag() == 21) {
	  ID theCompInfo = theSIFMember->getMemberInfo();
	  SIFCompartment* theComp = this->getSIFCompartment(theCompInfo(0) * 1000 + theCompInfo(1) * 100 + theCompInfo(2));
	  ID theZBeams = theComp->getConnectedYBeams();
	  SIFYBeam* SIFZBeam1 = this->getSIFYBeam(theZBeams(0));
	  SIFYBeam* SIFZBeam2 = this->getSIFYBeam(theZBeams(1));
	  ID IntNodes1 = SIFZBeam1->getIntNodeTags();
	  ID IntNodes2 = SIFZBeam2->getIntNodeTags();
	  double xCrd1, xCrd2, zCrd;
	  xCrd1 = Jt1->getCrds()(0);
	  xCrd2 = Jt2->getCrds()(0);
	  zCrd = Jt1->getCrds()(2);

	  for (int j = 0; j<IntNodes1.Size(); j++) {
		  int NodeTag1 = IntNodes1(j);
		  int NodeTag2 = IntNodes2(j);
		  Node* theNode1 = theDomain->getNode(NodeTag1);
		  Node* theNode2 = theDomain->getNode(NodeTag2);
		  if ((fabs(theNode1->getCrds()(0) - xCrd1)<1e-4) && (fabs(theNode1->getCrds()(2) - zCrd)<1e-4))
			  NodeJt1Tag = NodeTag1;
		  else if ((fabs(theNode2->getCrds()(0) - xCrd1)<1e-4) && (fabs(theNode2->getCrds()(2) - zCrd)<1e-4))
			  NodeJt1Tag = NodeTag2;

		  if ((fabs(theNode1->getCrds()(0) - xCrd2)<1e-4) && (fabs(theNode1->getCrds()(2) - zCrd)<1e-4))
			  NodeJt2Tag = NodeTag1;
		  else if ((fabs(theNode2->getCrds()(0) - xCrd2)<1e-4) && (fabs(theNode2->getCrds()(2) - zCrd)<1e-4))
			  NodeJt2Tag = NodeTag2;
	  }

	  int theDOFid[] = { 0,1,2,3,4,5 };
	  ID theConDOF = ID(theDOFid, 6);
	  ID theRetDOF = ID(theDOFid, 6);
	  Matrix Ccr(6, 6);
	  Ccr.Zero();
	  Ccr(0, 0) = 1;
	  Ccr(1, 1) = 1; Ccr(2, 2) = 1;
	  Ccr(3, 3) = 1;
	  Ccr(4, 4) = 1; Ccr(5, 5) = 1;

	  MP_Constraint* theMP = new MP_Constraint(Jt1->getNodeTag()(0), NodeJt1Tag, Ccr, theConDOF, theRetDOF);
	  theDomain->addMP_Constraint(theMP);

	  theMP = new MP_Constraint(Jt2->getNodeTag()(0), NodeJt2Tag, Ccr, theConDOF, theRetDOF);
	  theDomain->addMP_Constraint(theMP);
  }
  //----------------------end of secXbeam--------------------------------
  else if((theSIFMember->getMemberTypeTag()==11)|| ((theSIFMember->getMemberTypeTag()==1)&(isPinned))){
	  
	  theNodeTag++;
	  
	  NdCrdx = Jt1xCrd ;
	  if(theSIFMember->getMemberTypeTag()==11)
	  NdCrdy = Jt1yCrd+offset ;
	  else if(theSIFMember->getMemberTypeTag()==21)
	  NdCrdy = Jt1yCrd;
	  else if((theSIFMember->getMemberTypeTag()==1)&(isPinned))
	  NdCrdy = Jt1yCrd;
	  
	  NdCrdz = Jt1zCrd;
	  
	  intNodes.resize(NumEles+1);

	  Node *theNode = new Node(theNodeTag, nDoF, NdCrdx, NdCrdy, NdCrdz);	
	   theDomain -> addNode(theNode);
	
		intNodes(NumEles-1)= theNodeTag;
	   
	  theNodeTag++;
  	  NdCrdx = Jt2xCrd ;
	  if(theSIFMember->getMemberTypeTag()==11)
	  NdCrdy = Jt2yCrd+offset ;
	  else if(theSIFMember->getMemberTypeTag()==21)
	  NdCrdy = Jt2yCrd;
	  else if((theSIFMember->getMemberTypeTag()==1)&(isPinned))
	  NdCrdy = Jt2yCrd;
	  

	  NdCrdz = Jt2zCrd ;
	  theNode = new Node(theNodeTag, nDoF, NdCrdx, NdCrdy, NdCrdz);	  
	   theDomain -> addNode(theNode);

	  intNodes(NumEles)= theNodeTag;
	 
   }
   //SIFMember stores the nodes
  theSIFMember->setIntNodeTags(intNodes );

  ///-------------------------Adding rigidbeam connections for Slab in plane frame--------------------------------//
  //slab modelled as beam elements
  if (theSIFMember->getMemberTypeTag() == 11) {
	  //add muilti point constraints
	  //----------------Connecting joints------------------
	  //MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, 
	  // ID &constrainedDOF, 
	  // ID &retainedDOF, int clasTag)
	
	  ID theCompInfo = theSIFMember->getMemberInfo();
	  SIFCompartment* theSIFComp = this->getSIFCompartment(theCompInfo);
	  ID XBeams = theSIFComp->getConnectedXBeams();
	  ID YBeams = theSIFComp->getConnectedYBeams();
	  //rigidLink for connecting beam and slab;
	  SIFXBeam* theXBeam = 0;
	  SIFYBeam* theYBeam = 0;
	  ID theBIntNodes;
	  if (XBeams != 0) {
		  theXBeam = this->getSIFXBeam(XBeams(0));
		  theBIntNodes = theXBeam->getIntNodeTags();
	  }
	  if (YBeams != 0) {
		  theYBeam = this->getSIFYBeam(YBeams(0));
		  theBIntNodes = theYBeam->getIntNodeTags();
	  }

	  int theBeamNodeTag = 0; int theSlabNodeTag=0;
	  for (int i = 0; i<intNodes.Size(); i++) {

			  if (i == intNodes.Size() - 2) {
				  if (isPinned)
					  theBeamNodeTag = theBIntNodes(i);
				  else
					  theBeamNodeTag = NodeJt1Tag;
			  }
			  else if (i == intNodes.Size() - 1) {
				  if (isPinned)
					  theBeamNodeTag = theBIntNodes(i);
				  else
					  theBeamNodeTag = NodeJt2Tag;
			  }

			 else
				theBeamNodeTag = theBIntNodes(i);
		  

		  theSlabNodeTag = intNodes(i);

		  //const Vector crds1 = (theDomain->getNode(theBeamNodeTag))->getCrds();
		  //const Vector crds2 = (theDomain->getNode(theSlabNodeTag))->getCrds();
		  //if(crds1(0)!=crds2(0)||crds1(2)!=crds2(2))
		  //opserr<<"WARNING:: incompatiable nodes for defining rigidlink--Node1 "<<crds1<<"Node2  "<<crds2<<endln;
		  if(theBeamNodeTag!=0&& theSlabNodeTag!=0)
		  RigidBeam theLink(*theDomain, theBeamNodeTag, theSlabNodeTag);
	  }

	  //end of adding  rigidlink to connect XBeam 
  }
  // primary beams with pinned connections
  else if ((theSIFMember->getMemberTypeTag() == 1)&(isPinned)){
	int theDOFid[] = { 0,1,2 };
	ID theConDOF = ID(theDOFid, 3);
	ID theRetDOF = ID(theDOFid, 3);
	Matrix Ccr(3, 3);
	Ccr.Zero();
	Ccr(0, 0) = 1;
	Ccr(1, 1) = 1; Ccr(2, 2) = 1;
	
	MP_Constraint* theMP1= new MP_Constraint(Jt1->getNodeTag()(0), intNodes(NumEles-1), Ccr, theConDOF, theRetDOF);
	theDomain->addMP_Constraint(theMP1);
	MP_Constraint* theMP2 = new MP_Constraint(Jt2->getNodeTag()(0), intNodes(NumEles), Ccr, theConDOF, theRetDOF);
	theDomain->addMP_Constraint(theMP2);
	//Pinned Connections
  }

  ///////////////////////////////Mesh member in nodes///////////////////////////////////

  //---------------------Now adding eles for frame beam and columns--------------------------
  SectionForceDeformation* theSection =0;
  SIFSection* theSIFSection = 0;
  BeamIntegration* beamIntegr = 0;
  CrdTransf* theTransf3d  = 0;
  Element* theElement =0;
  
  ID eles = ID(NumEles);

 
  
 

  int theNodeTag1, theNodeTag2;
  for (int i = 0; i < NumEles; i++)
  {
	  //slab modelled as beam elements, or beams are pinned to column
	if((theSIFMember->getMemberTypeTag()==11)||((theSIFMember->getMemberTypeTag()==1)&isPinned)){
		if(i==0){
		theNodeTag1 = intNodes(NumEles-1);
		theNodeTag2 = intNodes(0);
		}
		else if(i==(NumEles-1)){
		theNodeTag1 = intNodes(NumEles-2);
		theNodeTag2 = intNodes(NumEles);
		}
		else{
		theNodeTag1 = intNodes(i-1);
		theNodeTag2 = intNodes(i);
		}

	}
	//other normal beams&columns
	else{
		if(i==0){
		theNodeTag1 = NodeJt1Tag;
		theNodeTag2 = intNodes(0);
		}
		else if(i==(NumEles-1)){
		theNodeTag1 = intNodes(NumEles-2);
		theNodeTag2 = NodeJt2Tag;
		}
		else{
		theNodeTag1 = intNodes(i-1);
		theNodeTag2 = intNodes(i);
		}
	}
    theEleTag++;
    theSIFSection =  theSIFMember->getSIFSectionPtr();
	//theSIFSection = this->getSIFSection(1);
	if(theSIFSection==0){
		opserr<<"WARNING::SIFBuilderDomain can't find the section assigned to Member: "<<theSIFMember->getTag()<<endln;
		return -1;
		}

    theSection = theSIFSection->DefineBeamSection();
	
	SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
	for (int j = 0; j < nIP; j++)
		sections[j] = theSection;
	
	beamIntegr = new LobattoBeamIntegration();   //Beam integration
	int theMemberTypeTag = theSIFMember->getMemberTypeTag();

	theTransf3d = this->getCrdTransf(theMemberTypeTag,1) ;
	//////

    //DispBeamColumn2dThermal(int tag, int nd1, int nd2, int numSections, SectionForceDeformation **s, BeamIntegration &bi, CrdTransf &coordTransf, double rho = 0.0);
    
    //have not finished yet
	double mass = 0.0;
   theElement = new DispBeamColumn3dThermal(theEleTag, theNodeTag1, theNodeTag2, nIP, sections, *beamIntegr, *theTransf3d, mass);
    
   theDomain->addElement(theElement);

    eles(i) = theEleTag;
    
  }
  //end of mesh eles
   theSIFMember->setIntEleTags( eles );
#ifdef _DEBUG
   opserr<< intNodes <<eles<<endln;
#endif

   if (theSIFMember->getMemberTypeTag() == 2 && theSIFMember->getTag() == 3302)
	   opserr << "Beam:" << theSIFMember->getTag() << eles << endln;
   if(theSIFMember->getMemberTypeTag()==3&&theSIFMember->getTag()==3301)
	   opserr<<"Column:"<< theSIFMember->getTag() <<eles<<endln;

// If it is a x-y plane frame of model type 12;
// and the slab is defined, the slab should be meshed as beamColumn elements


///////////////////////////////////////Slab in plane frame///////////////////////////////////////////


  return 0;    
}


////////////////----------------------MeshSIFSlab-----------------------------//////////
int
SIFBuilderDomain::MeshSIFSlab(SIFMember* theSIFMember, int NumEleX,int NumEleZ)
{
	int OriginNodeTag = theNodeTag;
	// obtain the connected joints
	ID theJoints = theSIFMember->getConnectedJoints();
	ID theCompInfo = theSIFMember->getMemberInfo();
	Matrix JtCrds = Matrix(4,3);

	for(int i = 0; i<4;i++){
		int jointTag = theJoints(i);
		SIFJoint* theJoint = this->getSIFJoint (jointTag);
		const Vector theCrds = theJoint->getCrds();
		for(int j=0;j<theCrds.Size();j++)
			JtCrds(i,j)= theCrds(j);
	}

	double crdXincr = 0; double crdZincr = 0;
	crdXincr = (JtCrds(1,0)- JtCrds(0,0))/NumEleX;
	crdZincr = (JtCrds(3,2)- JtCrds(0,2))/NumEleZ;

	//Generating internal nodes
	ID intNodes = ID((NumEleX+1)*(NumEleZ+1));  // for storing the interior nodes
  
    double NdCrdx, NdCrdy, NdCrdz;

	for(int i=0; i<= NumEleZ; i++){
		for(int j=0 ;j<=NumEleX ; j++){
			
			theNodeTag++;
			NdCrdx = JtCrds(0,0) + crdXincr*j;
			NdCrdy = JtCrds(0,1)+offset;
			//NdCrdy = JtCrds(0,1);
			NdCrdz = JtCrds(0,2) + crdZincr*i;

			Node *theNode = new Node(theNodeTag, nDoF, NdCrdx, NdCrdy, NdCrdz);
  
			theDomain -> addNode(theNode);
    
			intNodes((NumEleX+1)*i+j) = theNodeTag;

		}

	}
	theSIFMember->setIntNodeTags( intNodes );


	//ShellMITC4Thermal::ShellMITC4Thermal(  int tag, 
                         //int node1,int node2,int node3,int node4,
	                     //SectionForceDeformation &theMaterial )
	SectionForceDeformation* theShellSection = 0;
	//theShellSection
	//MembranePlateFiberSectionThermal::MembranePlateFiberSectionThermal( int tag,  double thickness, NDMaterial &Afiber ) 
	//ElasticIsotropic3DThermal(int tag, double e, double nu, double rho=0,double alpha=0.0);
	SIFSection* theSIFSection = theSIFMember->getSIFSectionPtr();
	//isElastic = true;
	theShellSection = theSIFSection->DefineShellSection(isElastic);
	Element* theElement =0;
	ID eles = ID(NumEleX*NumEleZ);
	int NodeTag1,NodeTag2,NodeTag3,NodeTag4;
	for(int i=0; i< NumEleZ; i++){
		for(int j=0 ;j<NumEleX ; j++){
			theEleTag++;
			NodeTag1 = intNodes((NumEleX+1)*i+j);
			NodeTag2 = intNodes((NumEleX+1)*i+j+1);
			NodeTag3 = intNodes((NumEleX+1)*(i+1)+j+1);
			NodeTag4 = intNodes((NumEleX+1)*(i+1)+j);
			
			theElement = new ShellMITC4Thermal(theEleTag, NodeTag1, NodeTag2,NodeTag3,NodeTag4,
				*theShellSection);	 
			theDomain->addElement(theElement);
			eles(NumEleX*i+j) = theEleTag;
		}
	}
	theSIFMember->setIntEleTags( eles );
#ifdef _DEBUG
   opserr<< eles<<endln;
#endif
	//add muilti point constraints
	//----------------Connecting joints------------------
   //MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, 
			    // ID &constrainedDOF, 
			    // ID &retainedDOF, int clasTag)
	int theDOFid[] = {0,1,2};
	ID theConDOF = ID(theDOFid,3);
	ID theRetDOF = ID(theDOFid,3);
	Matrix Ccr (3, 3);
	Ccr.Zero();
	Ccr(0,0)=1;
	Ccr(1,1)=1;Ccr(2,2)=1;
	//Ccr(3,3) =1;
	//Ccr(4,4) =1; Ccr(5,5) =1;

	for(int i = 0; i<4;i++){
		int jointTag = theJoints(i);
		SIFJoint* theJoint = this->getSIFJoint (jointTag);
		int theJtNodeTag = theJoint->getNodeTag()(0);
		int SlabJtNodeTag=0;

		if(i==0)
			SlabJtNodeTag =OriginNodeTag+1;
		else if(i==1)
			SlabJtNodeTag =OriginNodeTag+NumEleX+1;
		else if(i==2)
			SlabJtNodeTag =OriginNodeTag+(NumEleX+1)*(NumEleZ+1);
	    else 
			SlabJtNodeTag =OriginNodeTag+(NumEleX+1)*NumEleZ+1;

		//MP_Constraint* theMP= new MP_Constraint(theJtNodeTag, SlabJtNodeTag, theConDOF, theRetDOF);
		//MP_Constraint* theMP= new MP_Constraint(theJtNodeTag, SlabJtNodeTag, Ccr, theConDOF, theRetDOF);
		  //  theDomain->addMP_Constraint(theMP);
		 RigidBeam theLink(*theDomain, theJtNodeTag, SlabJtNodeTag);
	}
 
	SIFCompartment* theSIFComp = this->getSIFCompartment(theCompInfo);
	ID XBeams = theSIFComp->getConnectedXBeams();
	ID YBeams = theSIFComp->getConnectedYBeams();
	ID XSecBeams = theSIFComp->getConnectedSecXBeams();
	//rigidLink for connecting beam and slab;
	 ///////////////////////////////Connections betweem XBeams and slab ///////////////////////////////////
	for(int k=0; k<XBeams.Size();k++){
		int XBeamTag = XBeams(k);
		SIFXBeam* theXBeam = this->getSIFXBeam(XBeamTag);	
		ID theIntNodes = theXBeam->getIntNodeTags();

		for(int i=0; i<(theIntNodes.Size()); i++){
			int theBeamNodeTag = theIntNodes(i);
			int theSlabNodeTag;
			if(k==0)
				theSlabNodeTag = intNodes(i+1);
			else if(k==1)
				theSlabNodeTag = intNodes((NumEleX+1)*NumEleZ+i+1);

			const Vector crds1 = (theDomain->getNode(theBeamNodeTag))->getCrds();
			const Vector crds2 = (theDomain->getNode(theSlabNodeTag))->getCrds();
			if(crds1(0)!=crds2(0)||crds1(2)!=crds2(2))
				opserr<<"WARNING:: incompatiable nodes for defining rigidlink--Node1 "<<crds1<<"Node2  "<<crds2<<endln;
			RigidBeam theLink(*theDomain, theBeamNodeTag, theSlabNodeTag);
			
			//MP_Constraint* theMP= new MP_Constraint(theBeamNodeTag, theSlabNodeTag,Ccr, theConDOF, theRetDOF);
		    //theDomain->addMP_Constraint(theMP);
		}

	}
	 ///////////////////////////////Connections between YBeams and slab///////////////////////////////////
	for(int k=0; k<YBeams.Size();k++){
		int YBeamTag = YBeams(k);
		SIFYBeam* theYBeam = this->getSIFYBeam(YBeamTag);	
		ID theIntNodes = theYBeam->getIntNodeTags();

		for(int i=0; i<(theIntNodes.Size()); i++){
			int theBeamNodeTag = theIntNodes(i);
			int theSlabNodeTag;
			if(k==0)
				theSlabNodeTag = intNodes((NumEleX+1)*(i+1));
			else if(k==1)
				theSlabNodeTag = intNodes((NumEleX+1)*(i+2)-1);

			const Vector crds1 = (theDomain->getNode(theBeamNodeTag))->getCrds();
			const Vector  crds2 = (theDomain->getNode(theSlabNodeTag))->getCrds();
			if(crds1(0)!=crds2(0)||crds1(2)!=crds2(2))
				opserr<<"WARNING:: incompatiable nodes for defining rigidlink--Node1 "<<crds1<<"Node2  "<<crds2<<endln;
			RigidBeam theLink(*theDomain, theBeamNodeTag, theSlabNodeTag);
			//MP_Constraint* theMP= new MP_Constraint(theBeamNodeTag, theSlabNodeTag,Ccr, theConDOF, theRetDOF);
		    //theDomain->addMP_Constraint(theMP);
		}

	}


 ///////////////////////////////Connections between Secondary Beam and slabs///////////////////////////////////
	if(XSecBeams!=0){
		
	for(int k=0; k<XSecBeams.Size();k++){
		int SecXBeamTag = XSecBeams(k);
		SIFXBeamSec* theSecXBeam = this->getSIFXBeamSec(SecXBeamTag);	
		SIFJoint* theJT1  = this->getSIFJoint(theSecXBeam->getConnectedJoints()(0));
		SIFJoint* theJT2  = this->getSIFJoint(theSecXBeam->getConnectedJoints()(1));
		double SexBeamCrdZ = theJT1->getCrds()(2);
		int StartTaginSlb=0;
		for(int i =0; i<=NumEleZ; i++){
			int SlbNodeTag = intNodes((NumEleX+1)*i);
			Node* theNode = theDomain->getNode(SlbNodeTag);
			double crdZ = theNode->getCrds()(2);
			if(crdZ-SexBeamCrdZ<1e-4&&crdZ-SexBeamCrdZ>-1e-4)
				StartTaginSlb = SlbNodeTag;
		}
		if(StartTaginSlb==0)
			opserr<<"WARNING: Secondary beam is not identified when meshing slab "<<theSIFMember->getTag()
			      <<"! Please check the mesh control parametres"<<endln;
	
		//RigidBeam theLink1(*theDomain, BeamNode1,StartTaginSlb);
		//RigidBeam theLink2(*theDomain, BeamNode2,StartTaginSlb+NumEleX);

		ID theIntNodes = theSecXBeam->getIntNodeTags();

		for(int i=0; i<(theIntNodes.Size()); i++){
			int theBeamNodeTag = theIntNodes(i);
			int theSlabNodeTag;
			theSlabNodeTag = StartTaginSlb+i+1;

			const Vector crds1 = (theDomain->getNode(theBeamNodeTag))->getCrds();
			const Vector crds2 = (theDomain->getNode(theSlabNodeTag))->getCrds();
			if(crds1(0)!=crds2(0)||crds1(2)!=crds2(2))
				opserr<<"WARNING:: incompatiable nodes for defining rigidlink--Node1 "<<crds1<<"Node2  "<<crds2<<endln;
			//MP_Constraint* theMP= new MP_Constraint(theBeamNodeTag, theSlabNodeTag,Ccr, theConDOF, theRetDOF);
		    //theDomain->addMP_Constraint(theMP);
			RigidBeam theLink(*theDomain, theBeamNodeTag, theSlabNodeTag);
		}
		

	}
	}
	//for SecXBeam is not 0

	return 0;

}

CrdTransf* 
SIFBuilderDomain::getCrdTransf(int MemberType, int TransfType)
{
	Vector vecxzPlane(3);  
	vecxzPlane(0)=0;vecxzPlane(1)=0;vecxzPlane(2)=1;
	Vector jntOffsetI(3),jntOffsetJ(3);

	if(MemberType == 1||MemberType == 21){
		vecxzPlane(0)=0;vecxzPlane(1)=0; vecxzPlane(2)=1;
	   //001 for xbeam
	}
	else if(MemberType == 2){
	   vecxzPlane(0)=-1;vecxzPlane(1)=0; vecxzPlane(2)=0;
	   //
	}
	else if(MemberType == 3){
	   vecxzPlane(0)=1;vecxzPlane(1)=0; vecxzPlane(2)=0;
	   //
	}
	else if(MemberType == 11){
		if(this->getSIFBuilderInfo()(0) == 12){
			vecxzPlane(0)=0;vecxzPlane(1)=0; vecxzPlane(2)=1;
		}
		else if(this->getSIFBuilderInfo()(0) == 13){
			vecxzPlane(0)=-1;vecxzPlane(1)=0; vecxzPlane(2)=0;
		}
		else 
			opserr<<"WARNING::SIFBuilderDomain can not assign GeomTransf for 2D slab"<<endln;
	   //
	}

	else {
		opserr<<"WARNING:Invalid member type "<<MemberType<<" to sepcify the geometrical transformation"<<endln;
	}

	CrdTransf* theTransf3d = new  CorotCrdTransf3d(0, vecxzPlane, jntOffsetI, jntOffsetJ);

	return theTransf3d;

}




//-////------------------------------Apply Gravity Load-------------------------/////
int
SIFBuilderDomain::applySelfWeight(int thePatternTag, double dt)
{
	if(LoadApplied != 0)
		LoadApplied = 0;

	int theEleLoadTag = 0;
	int loadPatternTag = thePatternTag;
	
	//------------------Gravity Load-------------------
	LoadPattern *theLoadPattern0 = new LoadPattern(loadPatternTag);
	TimeSeries *theSeries0 = new LinearSeries();
    theLoadPattern0->setTimeSeries(theSeries0);
    theDomain->addLoadPattern(theLoadPattern0);

	Beam3dUniformLoad *theLoad = 0;

	SIFXBeam *theSIFXBeam;
	SIFXBeamIter &theSIFXBeams  = this->getSIFXBeams();
	while((theSIFXBeam = theSIFXBeams())!=0){

		SIFSection *theSIFSection = theSIFXBeam->getSIFSectionPtr();
		ID EleTags = theSIFXBeam->getIntEleTags();
		int NumEles = EleTags.Size();
		int sectionType = theSIFSection->getSectionTypeTag();

		if (sectionType == 1) {									//Rectangular section;
			Vector SectionPars(2);
			SectionPars = theSIFSection->getSectionPars();
			double Area = SectionPars(0) * SectionPars(1);
			double Weight = Area * 7850 * 9.8;					//SI unit i.e. meter, newton is assumed;
			
			for (int i = 0; i < NumEles; i++){
				theLoad = new Beam3dUniformLoad(theEleLoadTag, -Weight, 0.0, 0.0, EleTags(i));
				//opserr<<"Gravity applied at SIFXBeam "<<theSIFXBeam->getTag()<<" with Element "<<EleTags(i) <<endln;

				if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
				//opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theLoad;
				return -2;
				}

				theEleLoadTag++;
			}
		}
		else if (sectionType == 2) {							//I section;
			Vector SectionPars(4);
			SectionPars = theSIFSection->getSectionPars();
			double Area = 2 * SectionPars(2) * SectionPars(3) + (SectionPars(0) - 2 * SectionPars(3)) * SectionPars(1);
			double Weight = Area * 7850 * 9.8+5000;					//SI unit i.e. meter, newton is assumed;
			
			for (int i = 0; i < NumEles; i++){
				theLoad = new Beam3dUniformLoad(theEleLoadTag, -Weight, 0.0, 0.0, EleTags(i));
				//opserr<<"Gravity applied at SIFXBeam "<<theSIFXBeam->getTag()<<" with Element "<<EleTags(i) <<endln;

				if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theLoad;
				return -2;
				}

				theEleLoadTag++;
			}
		}
	}
  
	SIFYBeam *theSIFYBeam;
	SIFYBeamIter &theSIFYBeams  = this->getSIFYBeams();
	while((theSIFYBeam = theSIFYBeams())!=0){

		SIFSection *theSIFSection = theSIFYBeam->getSIFSectionPtr();
		ID EleTags = theSIFYBeam->getIntEleTags();
		int NumEles = EleTags.Size();
		int sectionType = theSIFSection->getSectionTypeTag();

		if (sectionType == 1) {									//Rectangular section;
			Vector SectionPars(2);
			SectionPars = theSIFSection->getSectionPars();
			double Area = SectionPars(0) * SectionPars(1);
			double Weight = Area * 7850 * 9.8;					//SI unit i.e. meter, newton is assumed;
			
			for (int i = 0; i < NumEles; i++){
				theLoad = new Beam3dUniformLoad(theEleLoadTag, -Weight, 0.0, 0.0, EleTags(i));
				//opserr<<"Gravity applied at SIFYBeam "<<theSIFYBeam->getTag()<<" with Element "<<EleTags(i) <<endln;

				if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theLoad;
				return -2;
				}

				theEleLoadTag++;
			}
		}
		else if (sectionType == 2) {							//I section;
			Vector SectionPars(4);
			SectionPars = theSIFSection->getSectionPars();
			double Area = 2 * SectionPars(2) * SectionPars(3) + (SectionPars(0) - 2 * SectionPars(3)) * SectionPars(1);
			double Weight = Area * 7850 * 9.8+5000;					//SI unit i.e. meter, newton is assumed;
			
			for (int i = 0; i < NumEles; i++){
				theLoad = new Beam3dUniformLoad(theEleLoadTag, -Weight, 0.0, 0.0, EleTags(i));
#ifdef _DEBUG
				//opserr<<"Gravity applied at SIFYBeam "<<theSIFYBeam->getTag()<<" with Element "<<EleTags(i) <<endln;
#endif
				if (theDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theLoad;
				return -2;
				}

				theEleLoadTag++;
			}
		}
	}
//end of selfweight YBeam;

	//opserr << *theDomain;
	//------------end of gravity Load------------------
    int numSteps = 1/dt;
	double TimeStep = dt;

	opserr<<"SIFBuilder::Now applying gravity load in " << numSteps <<" steps"<<endln;

    if(theAnalysisModel== 0)
	  theAnalysisModel = new AnalysisModel();

	if(theStaticIntegrator == 0){
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
		//(double dLambda, int numIncr, double min, double max)
	} 
	else {
		delete theStaticIntegrator;
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
	}
	
	if(theTest == 0)
		theTest =  new CTestNormDispIncr(1e-3, 1000, 0);

	if(theAlgorithm==0)
		theAlgorithm = new NewtonRaphson(*theTest);
    
	if(theHandler==0)
		theHandler = new PlainHandler(); 
    
	if(theNumberer==0){
		RCM *theRCM = new RCM(false);
		theNumberer = new DOF_Numberer(*theRCM);  
	}

	if(theSOE ==0){
		ProfileSPDLinDirectSolver* theSolver = new ProfileSPDLinDirectSolver();
		theSOE = new ProfileSPDLinSOE(*theSolver);   
	}

    theStaticAnalysis = new StaticAnalysis(*theDomain,
				  *theHandler,
				  *theNumberer,
				  *theAnalysisModel,
				  *theAlgorithm,
				  *theSOE,
				  *theStaticIntegrator,
				  theTest);

    // perform the analysis & print out the results for the domain
    theStaticAnalysis->analyze(numSteps);
    //opserr << *theDomain;

  LoadApplied++;
  return 0;
}



///////////-----------------------Apply Miscellaneous Load------------------------/////-
int
SIFBuilderDomain::applyMiscLoad(int thePatternTag, double dt)
{
	//LoadType = 0: All the load will be applied
	//1: Gravity Load, 2:MiscLoad (NodalLoad or elemental load), 3:FireAction
	
	//------------------Misc Load-------------------
	//LoadPattern *theLoadPattern = new LoadPattern(2);
	// theLoadPattern->setTimeSeries(theSeries);
    //theDomain->addLoadPattern(theLoadPattern);
	 if(LoadApplied != 0){
	  theDomain->setCurrentTime(0.0);
	  theDomain->setCommittedTime(0.0); 
	  theDomain->setLoadConstant();
     }

	 theNodalLoadTag =0;
	int loadPatternTag = thePatternTag;
	
	//------------------Nodal Load-------------------
	LoadPattern *theLoadPattern = new LoadPattern(loadPatternTag);
	TimeSeries *theSeries = new LinearSeries();
    theLoadPattern->setTimeSeries(theSeries);
    theDomain->addLoadPattern(theLoadPattern);

	NodalLoad *theLoad = 0;
	ElementalLoad *theEleLoad = 0;
	
	SIFJoint *theSIFJoint;
    SIFJointIter &theSIFJoints  = this->getSIFJoints();
	while((theSIFJoint = theSIFJoints())!=0){
    Vector theJointLoad = theSIFJoint->getLoad();
	 if(theJointLoad!=0){
		int theJointNodeTag = (theSIFJoint->getNodeTag())(0);
		theLoad = new NodalLoad(theNodalLoadTag, theJointNodeTag, theJointLoad, false);
		opserr<<"Load applied at SIFJoint "<<theSIFJoint->getTag()<<" with Node "<<theJointNodeTag <<endln;

		//theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
		if (theDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
			opserr << "WARNING TclModelBuilder - could not add load to domain\n";
			delete theLoad;
			return -2;
		}
		theNodalLoadTag++;
	 }
	}
	
//Applying Misc load for xBeams
	SIFXBeam *theSIFXBeam;
	SIFXBeamIter &theSIFXBeams  = this->getSIFXBeams();
	while((theSIFXBeam = theSIFXBeams())!=0){
		Vector theXBeamLoad = theSIFXBeam->getLoad();
		if(theXBeamLoad!=0){
		//ID theIntEles = theSIFXBeam->getIntEleTags();
			ID theIntNodes = theSIFXBeam->getIntNodeTags();
		int NumNodes = theIntNodes.Size(); 
		for (int i = 0; i < NumNodes; i++){
			theLoad = new NodalLoad(theNodalLoadTag, theIntNodes(i), theXBeamLoad, false);
#ifdef _DEBUG
				//opserr<<"Miscload applied at SIFYBeam "<<theSIFYBeam->getTag()<<" with Element "<<EleTags(i) <<endln;
#endif
			if (theDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add load to domain\n";
				delete theLoad;
				return -2;
			}

				theNodalLoadTag++;
		}
	 }
		
	}
//Applying Misc load for xBeams
	SIFYBeam *theSIFYBeam;
	SIFYBeamIter &theSIFYBeams  = this->getSIFYBeams();
	while((theSIFYBeam = theSIFYBeams())!=0){
		Vector theYBeamLoad = theSIFYBeam->getLoad();
		ID theIntEles = theSIFYBeam->getIntEleTags();
		int NumEles = theIntEles.Size(); 
		if(theYBeamLoad!=0){
		for (int i = 0; i < NumEles; i++){
				theEleLoad = new Beam3dUniformLoad(theEleLoadTag, theYBeamLoad(1), theYBeamLoad(0), theYBeamLoad(2), theIntEles(i));
#ifdef _DEBUG
				//opserr<<"Miscload applied at SIFYBeam "<<theSIFYBeam->getTag()<<" with Element "<<EleTags(i) <<endln;
#endif
				if (theDomain->addElementalLoad(theEleLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theEleLoad;
				return -2;
				}

				theEleLoadTag++;
		}
	  }
		
	}
//Applying Misc load for slabs
	SIFSlab *theSIFSlab;
    SIFSlabIter &theSIFSlabs  = this->getSIFSlabs();
	while((theSIFSlab = theSIFSlabs())!=0){
		Vector theSlabLoad = theSIFSlab->getLoad();
		ID theIntEles = theSIFSlab->getIntEleTags();
		Element* theEle = theDomain->getElement(theIntEles(0));
		if(theEle !=0){
		ID theConnectedNodes = theEle->getExternalNodes();
			if(theConnectedNodes.Size()==4){
			const Vector crd1 = (theDomain->getNode(theConnectedNodes(0)))->getCrds();
			const Vector crd2 = (theDomain->getNode(theConnectedNodes(2)))->getCrds();
			//opserr<<crd1<<crd2<<endln;
			double EleArea = fabs((crd2(2)-crd1(2))*(crd2(0)-crd1(0)));
			theSlabLoad*=EleArea;
			}
		}
#ifdef _DEBUG
	opserr<<theSlabLoad;
#endif
	
	 if(theSlabLoad!=0){
		 if(theSIFSlab->getMemberTypeTag()==10){
		 ID theSlabNodeTags = theSIFSlab->getIntNodeTags();
		 for(int i=0; i<theSlabNodeTags.Size(); i++){
			int theSlabNodeTag = theSlabNodeTags(i);
			theLoad = new NodalLoad(theNodalLoadTag, theSlabNodeTag, theSlabLoad, false);
			//opserr<<"Load applied at SIFJoint "<<theSIFJoint->getTag()<<" with Node "<<theJointNodeTag <<endln;
			
			if (theDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
			opserr << "WARNING TclModelBuilder - could not add load to domain\n";
			delete theLoad;
			return -2;
			}
			theNodalLoadTag++;
		}
		//theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst)
	  }
		 //if slab is 3D;
		 else if(theSIFSlab->getMemberTypeTag()==11){
			ID theIntEles = theSIFSlab->getIntEleTags();
			int NumEles = theIntEles.Size(); 
			for (int i = 0; i < NumEles; i++){
				theEleLoad = new Beam3dUniformLoad(theEleLoadTag, theSlabLoad(1), theSlabLoad(0), theSlabLoad(2), theIntEles(i));
#ifdef _DEBUG
				//opserr<<"Miscload applied at SIFYBeam "<<theSIFYBeam->getTag()<<" with Element "<<EleTags(i) <<endln;
#endif
				if (theDomain->addElementalLoad(theEleLoad, loadPatternTag) == false) {
				opserr << "WARNING TclModelBuilder - could not add gravity to domain\n";
				delete theEleLoad;
				return -2;
				}

				theEleLoadTag++;
			}
		  }
		}
		 //if the SlabLoad!=0
	}
	 //while SIFSlab exists


	//opserr << *theDomain;
	//------------end of Nodal Load------------------
    int numSteps = 1/dt;
	double TimeStep = dt;

	opserr<<"SIFBuilder;:Now applying miscellaneous load in " << numSteps <<" steps"<<endln;
    if(theAnalysisModel== 0)
	  theAnalysisModel = new AnalysisModel();

	if(theStaticIntegrator == 0){
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
		//(double dLambda, int numIncr, double min, double max)
	} 
	else {
		delete theStaticIntegrator;
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
	}
	
	if(theTest == 0)
		theTest =  new CTestNormDispIncr(1e-3, 500, 1);

	if(theAlgorithm==0)
		theAlgorithm = new NewtonRaphson(*theTest);
    
	if(theHandler==0)
		//theHandler = new PlainHandler(); 
		theHandler  = new PenaltyConstraintHandler( 1.0e10, 1.0e10);
    
	if(theNumberer==0){
		RCM *theRCM = new RCM(false);
		theNumberer = new DOF_Numberer(*theRCM);  
	}

	//if(theSOE ==0){
	//ProfileSPDLinDirectSolver* theSolver = new ProfileSPDLinDirectSolver();
	//theSOE = new ProfileSPDLinSOE(*theSolver);   
	//}
	double thresh = 0.0;
	int npRow = 1;
	int npCol = 1;
	int np = 1;
	int permSpec = 0;
	int panelSize = 6;
	int relax = 6;
	char symmetric = 'N';
	double drop_tol = 0.0;
	
	SparseGenColLinSolver* theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
	theSOE = new SparseGenColLinSOE(*theSolver);
	



    theStaticAnalysis = new StaticAnalysis(*theDomain,
				  *theHandler,
				  *theNumberer,
				  *theAnalysisModel,
				  *theAlgorithm,
				  *theSOE,
				  *theStaticIntegrator,
				  theTest);

    // perform the analysis & print out the results for the domain
    theStaticAnalysis->analyze(numSteps);
    //opserr << *theDomain;

  LoadApplied++;
  return 0;
}

//---------------------------Apply Fire Action--------------------
int
SIFBuilderDomain::applyFireAction(int thePatternTag, double timeStep, double fireDuration)
{
	//if there was load applied
	if (LoadApplied != 0) {
		theDomain->setCurrentTime(0.0);
		theDomain->setCommittedTime(0.0);
		theDomain->setLoadConstant();
	}

	int loadPatternTag = thePatternTag;
	double FireDuration = fireDuration;
	double TimeStep = timeStep;

	LoadPattern *theLoadPattern = new LoadPattern(loadPatternTag);
	TimeSeries *theSeries = new LinearSeries();
	theLoadPattern->setTimeSeries(theSeries);
	theDomain->addLoadPattern(theLoadPattern);

	//Apply SIFfireAction
	opserr << "SIFBuilder;:Now analyzing the heat transfer to the structural members..." << endln;
	SIFfireAction *theSIFfireAction;
	SIFfireActionIter &theSIFfireActions = this->getSIFfireActions();
	while ((theSIFfireAction = theSIFfireActions()) != 0) {
		theSIFfireAction->Apply(loadPatternTag, timeStep, fireDuration);
	}



	//---------------------------------------Analysis for Fire Action----------------------
	int numSteps = FireDuration / TimeStep;

	opserr << "SIFBuilder;:Now applying fire action in " << numSteps << " steps" << endln;
	if (theAnalysisModel == 0)
		theAnalysisModel = new AnalysisModel();

	if (theStaticIntegrator == 0) {
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
		//(double dLambda, int numIncr, double min, double max)
	}
	else {
		delete theStaticIntegrator;
		theStaticIntegrator = new LoadControl(TimeStep, numSteps, TimeStep, TimeStep);
	}

	//if(theTest == 0)
	theTest = new CTestNormDispIncr(5e-3, 2000, 1);
	//theTest = new CTestNormUnbalance(1e-3, 200, 1, 2);


	if (theAlgorithm == 0)
		theAlgorithm = new NewtonRaphson(*theTest);

	if (theHandler == 0)
		//theHandler  = new PlainHandler();
		theHandler = new PenaltyConstraintHandler(1.0e10, 1.0e10);

	if (theNumberer == 0) {
		RCM *theRCM = new RCM(false);
		theNumberer = new DOF_Numberer(*theRCM);
	}

	//if(theSOE ==0){
		//ProfileSPDLinDirectSolver* theSolver = new ProfileSPDLinDirectSolver();
		//theSOE = new ProfileSPDLinSOE(*theSolver);   
	//}
	double thresh = 0.0;
	int npRow = 1;
	int npCol = 1;
	int np = 1;
	int permSpec = 0;
	int panelSize = 6;
	int relax = 6;
	char symmetric = 'N';
	double drop_tol = 0.0;
	if (theSOE == 0) {
		SparseGenColLinSolver* theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
		theSOE = new SparseGenColLinSOE(*theSolver);
	}
	else {
		delete theSOE;
		SparseGenColLinSolver* theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric);
		theSOE = new SparseGenColLinSOE(*theSolver);
	}


    theStaticAnalysis = new StaticAnalysis(*theDomain,
				  *theHandler,
				  *theNumberer,
				  *theAnalysisModel,
				  *theAlgorithm,
				  *theSOE,
				  *theStaticIntegrator,
				  theTest);

    // perform the analysis & print out the results for the domain
    int result= theStaticAnalysis->analyze(numSteps);
	
    //opserr << *theDomain;

  LoadApplied++;
  return 0;

  
}


int 
SIFBuilderDomain::getNodalLoadTag(){
	return theNodalLoadTag ;
}
	
void 
SIFBuilderDomain::incrNodalLoadTag(){
	theNodalLoadTag++;
}

int
SIFBuilderDomain::getEleLoadTag(){
	return theEleLoadTag;
}
void
SIFBuilderDomain::incrEleLoadTag(){
	theEleLoadTag++;
}