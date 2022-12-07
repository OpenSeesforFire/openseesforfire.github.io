/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)								  **
**   Liming Jiang (Liming.Jiang@ed.ac.uk)                             **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 2.4.2 $


// This file is added by Liming Jiang for creating Tcl interface to Heat Transfer analysis
// Written by  Liming Jiang (Liming.Jiang@polyu.edu.hk)

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>	
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
// includes for the simple mesh classes
#include <Simple_Line.h>
#include <Simple_Brick.h>
#include <Simple_Block.h>
#include <Simple_Isection.h>
#include <Simple_IsecProtected.h>
#include <Simple_Isection3D.h>
#include <Simple_Boundary.h>
#include <HTConstants.h>

#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <QuadFour.h>
#include <BrickEight.h>
#include <TemperatureBC.h>
#include <MP_TemperatureBC.h>
#include <Convection.h>


#include <CarbonSteelEC3.h>
#include <ConcreteEC2.h>
#include <SimpleMaterial.h>
#include <SFRMCoating.h>
#include <TimberHTMaterial.h>

// includes for the analysis classes
#include <HT_TransientAnalysis.h>
#include <BackwardDifference.h>
#include <HT_AnalysisModel.h> 
#include <HT_SolutionAlgorithm.h>
#include <LinearAlgorithm.h>
#include <NewtonMethod.h>
//#include <ModifiedNewtonMethod.h>
#include <CTestNormTempIncr.h>
#include <PenaltyBC_Handler.h>
#include <RCM.h>
#include <HT_DOF_Numberer.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <SuperLU.h>
#include <SparseGenColLinSOE.h>

//include fire models
#include <ParametricFireEC1.h>
#include <NorminalFireEC1.h>
#include <LocalizedFireEC1.h>
#include <Idealised_Local_Fire.h>
#include <LocalizedFireSFPE.h>
#include <AlpertCeilingJetModel.h>
#include <NaturalFire.h>
#include <UserDefinedFire.h>

#include <BoundaryPattern.h>
#include <FireImposedPattern.h>

//include recorders
#include <HTNodeRecorder.h>
#include <HTRecorderToStru.h>

#include <elementAPI.h>
#include <HeatTransferModule.h>
#include <PythonWrapper.h>
 #include <StandardStream.h>
 #include <DataFileStream.h>

 #include <SimulationInformation.h>
#include < PathTimeSeriesThermal.h>
//#include <OpenSeesCommands.h>

static PythonWrapper* theWrapper = 0;
 extern SimulationInformation simulationInfo;
 static HeatTransferModule* theHTModule = 0;
static HeatTransferDomain* theHTDomain = 0;
static BoundaryPattern* theHTPattern=0;

//Analysis classes
static HT_AnalysisModel *theAnalysisModel =0;
static HT_SolutionAlgorithm *theAlgorithm =0;
static TemperatureBCHandler *theHandler =0;
static HT_DOF_Numberer *theNumberer =0;
static LinearSOE *theSOE =0;
static HT_TransientAnalysis* theHTAnalysis =0;
static HT_TransientIntegrator* theTransientIntegrator=0;
static HT_ConvergenceTest* theTest =0;
static int HTReorderTag =0;

static double lasttime = 0;

//Declaration of Heat transfer Tcl commands
//void* OPS_addHTMaterial();


HeatTransferModule::HeatTransferModule(int ndm)
{
  
 // Tcl_CreateCommand(interp, "HTNode",	TclHeatTransferModule_addHTNode,(ClientData)NULL, NULL);
  opserr<<"               ----HeatTransfer Module is developed by UoE Group----       "<<endln;


  theHTModule = this;

  theHTMaterials = new ArrayOfTaggedObjects(10);
  theHTEntities = new ArrayOfTaggedObjects(10);
  theHTCons     = new ArrayOfTaggedObjects(5);
  theHTEleSets = new ArrayOfTaggedObjects(10);
  theHTNodeSets = new ArrayOfTaggedObjects(10);
  theFireModels = new ArrayOfTaggedObjects(10);
  
  theHTMeshes = new MapOfTaggedObjects();
  theHTMeshIter = new Simple_MeshIter(theHTMeshes);
  
  if(theHTDomain==0) {
	  theHTDomain= new HeatTransferDomain();
  }
  //if(theTclHTBoudary==0) {
	  //theTclHTBoudary= new Simple_Boundary();
  //}




}

HeatTransferModule::~HeatTransferModule()
{
	theHTDomain->clearAll();
  
  if (theHTMaterials!=0)
    delete theHTMaterials;
  
  if (theHTEntities !=0)
    delete theHTEntities;
  
  if (theHTCons != 0){
	  theHTCons->clearAll();
	  delete theHTCons;
  }
  
  if (theHTEleSets != 0){
	  theHTEleSets->clearAll();
	  delete theHTEleSets;
  }

  if  (theHTNodeSets != 0){
	  theHTNodeSets->clearAll();
	  delete theHTNodeSets;
  }

  if(theFireModels != 0){
	  theFireModels->clearAll();
	delete theFireModels;
  }
  
if (theHTMeshes != 0){
	theHTMeshes->clearAll();
    delete theHTMeshes;
}
  
  if (theHTMeshIter != 0)
    delete theHTMeshIter;
  

 theHTDomain = 0;
 theHTModule =0;
 theHTPattern=0;
 theHTAnalysis->clearAll();

 //Analysis classes
  theAnalysisModel =0;
  theAlgorithm =0;
  theHandler =0;
theNumberer =0;
theSOE =0;
theHTAnalysis =0;
 theTransientIntegrator=0;
theTest =0;
 HTReorderTag =0;
  
}


int 
HeatTransferModule::addHTMaterial(HeatTransferMaterial &theHTMaterial)
{
  bool result = theHTMaterials->addComponent(&theHTMaterial);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addNDMaterial() - failed to add material: " << theHTMaterial;
    return -1;
  }
}

HeatTransferMaterial *
HeatTransferModule::getHTMaterial(int tag)
{
  TaggedObject *mc = theHTMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  HeatTransferMaterial *result = (HeatTransferMaterial *)mc;
  return result;
}

//Proc for adding entity into OpenSees HeatTransfer module
int 
HeatTransferModule::addHTEntity(Simple_Entity &theHTEntity)
{
  bool result = theHTEntities->addComponent(&theHTEntity);
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addHTEntity() - failed to add Entity: " << theHTEntity;
    return -1;
  }
}

Simple_Entity *
HeatTransferModule::getHTEntity(int tag)
{
  TaggedObject *mc = theHTEntities->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  Simple_Entity *result = (Simple_Entity *)mc;
  return result;
}

//for adding Meshes into OpenSees HeatTransfer module
int 
HeatTransferModule::addHTMesh(Simple_Mesh *theHTMesh)
{
  
  // first check if a mesh with a similar tag exists in model
  int tag = theHTMesh->getTag();
  TaggedObject *other = theHTMeshes->getComponentPtr(tag);
  if (other != 0) {
    opserr << "HeatTransferModule::addHTMesh - cannot add as HTMesh with tag" <<
    tag << "already exists in model\n";
    
    return false;
  }
  
  // now we add the load pattern to the container for load pattrens
  bool result = theHTMeshes->addComponent(theHTMesh);
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addHTMesh() - failed to add Mesh: " << theHTMesh;
    return -1;
  }
}

Simple_Mesh *
HeatTransferModule::getHTMesh(int tag)
{
  TaggedObject *mc = theHTMeshes->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  Simple_Mesh *result = (Simple_Mesh *)mc;
  return result;
}

//To add the pointer to HTConstants into the arraryofTaggedObjects
//for adding HTConstants into OpenSees HeatTransfer module
int
HeatTransferModule::addHTConstants(HTConstants* theHTConstants)
{
  
  // now we add the load pattern to the container for load pattrens
  bool result = theHTCons->addComponent(theHTConstants);
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addHTConstants() - failed to add HTConstants: " << theHTConstants;
    return -1;
  }
}

// To get the pointer back
HTConstants*
HeatTransferModule::getHTConstants(int tag)
{
  TaggedObject *mc = theHTCons->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  HTConstants *result = (HTConstants *)mc;
  return result;
}

//for adding --HTElesSet--- into OpenSees HeatTransfer module
int
HeatTransferModule::addHTEleSet(HTEleSet* theHTEleSet)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theHTEleSets->addComponent(theHTEleSet);
  
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addHTEleSet() - failed to add HTEleSet: " << theHTEleSet;
    return -1;
  }
}

HTEleSet*
HeatTransferModule::getHTEleSet(int tag)
{
  TaggedObject *mc = theHTEleSets->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  HTEleSet *result = (HTEleSet *)mc;
  return result;
}

//for adding --HTNodesSet--- into OpenSees HeatTransfer module
int
HeatTransferModule::addHTNodeSet(HTNodeSet* theHTNodeSet)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theHTNodeSets->addComponent(theHTNodeSet);
  
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addHTNodeSet() - failed to add HTNodeSet: " << theHTNodeSet;
    return -1;
  }
}

HTNodeSet*
HeatTransferModule::getHTNodeSet(int tag)
{
  TaggedObject *mc = theHTNodeSets->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  HTNodeSet *result = (HTNodeSet *)mc;
  return result;
}

//for adding --FireModel--- into OpenSees HeatTransfer module
int
HeatTransferModule::addFireModel(FireModel* theFireModel)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theFireModels->addComponent(theFireModel);
  
  if (result == true)
    return 0;
  else {
    opserr << "HeatTransferModule::addFireModel() - failed to add FireModel: " << theFireModel;
    return -1;
  }
}

FireModel*
HeatTransferModule::getFireModel(int tag)
{
  TaggedObject *mc = theFireModels->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  FireModel *result = (FireModel *)mc;
  return result;
}




Simple_MeshIter &
HeatTransferModule::getHTMeshs()
{
  theHTMeshIter->reset();
  return *theHTMeshIter;
}


HeatTransferDomain*
HeatTransferModule::getHeatTransferDomain(){
  
  return theHTDomain;
}
//**********************************************************************************************//
///-----------------------------------Heat Transfer Functions----------------------------------//


int
OPS_HeatTransfer()
{

	int ndm = 0;
	if (OPS_GetNumRemainingInputArgs() < 1) {
		ndm = 2;
	}
	const char* type = OPS_GetString();

	if ((strcmp(type, "1d") == 0) || (strcmp(type, "1D") == 0)) {
		ndm = 1;
	}
	else if ((strcmp(type, "2d") == 0) || (strcmp(type, "2D") == 0)) {
		ndm = 2;
	}
	else if ((strcmp(type, "3d") == 0) || (strcmp(type, "3D") == 0)) {
		ndm = 3;
	}
	else
		opserr << "WARNING: HeatTransfer recieves invalid ndm tag" << endln;

	theHTModule = new HeatTransferModule(ndm);



	return 0;
}

int OPS_addHTMaterial()
{
	// num args
	
	if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "WARNING insufficient args: HTMaterial type tag\n";
		return -1;
	}

	// model type
	const char* HTmatType = OPS_GetString();

	int HTMaterialTag;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &HTMaterialTag) < 0) {
		opserr << "WARNING:HeatTransferMaterial failed to identify material tag\n";
		return -1;
	}
		

	//opserr << "material type"<<HTmatType<<"material tag:" << matTag;
	
	
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	HeatTransferMaterial* theHTMaterial = 0;

	//Adding CarbonSteelEC3
	if (strcmp(HTmatType, "SteelEC3") == 0 || strcmp(HTmatType, "CarbonSteelEC3") == 0) {

		theHTMaterial = new CarbonSteelEC3(HTMaterialTag);

	}
	//adding test material
	else if (strcmp(HTmatType, "Test") == 0 || strcmp(HTmatType, "test") == 0) {
		int size = 1;
		//opserr << "python test: " << size << endln;
		int* sizeList = &size;
		double* data = new double [50];
		OPS_GetListInput(sizeList, data);
		opserr <<"size "<< *sizeList<< ", python data: "<<data[1] << endln;
		theHTMaterial = new CarbonSteelEC3(HTMaterialTag);
	}
	//Adding ConcreteEC2
	else if(strcmp(HTmatType, "ConcreteEC2") == 0) {

		double moisture = 0;
		bool isLower = false;
		int numArgs = OPS_GetNumRemainingInputArgs();
		if (numArgs <1) {
			moisture = 0;
		}
		else {
			if (OPS_GetDoubleInput(&numdata, &moisture) < 0) {
				opserr << "WARNING:HeatTransferMaterial failed to identify moisture for material"<<HTMaterialTag<<endln;
				return -1;
			}

			if (OPS_GetNumRemainingInputArgs()>0) {
				const char* LowerInd = OPS_GetString();
				if (strcmp(LowerInd, "Lower") == 0 || strcmp(LowerInd, "lower") == 0) {
					isLower = true;
				}
			}
			
		}
	

		theHTMaterial = new ConcreteEC2(HTMaterialTag, moisture, isLower);
	}
	else if (strcmp(HTmatType, "SFRMCoating") == 0 || strcmp(HTmatType, "SFRM") == 0) {

		int typeTag = 0;

		if (OPS_GetNumRemainingInputArgs()<2) {
			typeTag = 1;
		}
		else {
			
			if (OPS_GetIntInput(&numdata,&typeTag) <0) {
				opserr << "WARNING invalid typeTag" << endln;
				opserr << " for HeatTransfer material: " << HTMaterialTag << endln;
				return -1;
			}
		
		
		}

		theHTMaterial = new SFRMCoating(HTMaterialTag, typeTag);

	}
	else if (strcmp(HTmatType, "Timber") == 0 || strcmp(HTmatType, "timber") == 0) {

		int typeTag = 0;
		double* paras = new double[6];
		if (OPS_GetNumRemainingInputArgs() < 2) {
			typeTag = 1;
		}
		else {

			if (OPS_GetIntInput(&numdata, &typeTag) < 0) {
				opserr << "WARNING invalid typeTag" << endln;
				opserr << " for HeatTransfer material: " << HTMaterialTag << endln;
				return -1;
			}
			if (OPS_GetNumRemainingInputArgs() > 5) {
				numdata = 6;
				if (OPS_GetDoubleInput(&numdata, paras) < 0) {
					opserr << "WARNING:: HTMaterial failed to identify parameters for material: " << HTmatType << endln;
					return -1;
				}
			}


		}

		Vector Pars = Vector(paras, 6);
		theHTMaterial = new TimberHTMaterial(HTMaterialTag, typeTag, theHTDomain, Pars);

	}
	else if (strcmp(HTmatType, "GenericMaterial") == 0) {

		double density = 0;
		double cp = 0;
		double conduct = 0;

		if (OPS_GetNumRemainingInputArgs() > 2) {
			if (OPS_GetDoubleInput(&numdata, &density) < 0) {
				opserr << "WARNING:: HTMaterial failed to identify density for material: " << HTmatType << endln;
				return -1;
			}

			if (OPS_GetDoubleInput(&numdata, &cp) < 0) {
				opserr << "WARNING:: HTMaterial failed to identify specific heat for material: " << HTmatType << endln;
				return -1;
			}
			
			if (OPS_GetDoubleInput(&numdata, &conduct) < 0){
				opserr << "WARNING:: HTMaterial failed to identify conductivity for material: " << HTmatType << endln;
				return -1;
			}
		}
		else
			opserr << "WARNING:: Defining HeatTransfer material should specify density, specific heat, and conductivity" << "\n";
		//SimpleMaterial(int tag, double rho, double cp, double kc); 
		theHTMaterial = new SimpleMaterial(HTMaterialTag, density, cp, conduct);

	}

	if (theHTMaterial != 0) {
		theHTModule->addHTMaterial(*theHTMaterial);
	}
	else {
		opserr << "WARNING: TclHTModule fail to add HeatTransfer Material: " << HTmatType << endln;
		return -1;
	}
		

return 0;
}


//Adding simpleEntity
int
OPS_addHTEntity()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 4) {
		opserr << "WARNING:: Creating an entity requires at least 7 arguments." << "\n";
		return -1;
	}

	// model type
	const char* EntityType = OPS_GetString();


	Simple_Entity* theHTEntity = 0;
	int HTEntityTag = 0;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &HTEntityTag) < 0) return -1;


	//Adding 2D entity:Block
	if (strcmp(EntityType, "Line") == 0 || strcmp(EntityType, "Line1D") == 0 || strcmp(EntityType, "Line1d") == 0)
	{
		double centerX = 0.0; double lengthX = 0.0;
		if (OPS_GetDoubleInput(&numdata, &centerX) < 0) {
			opserr << "WARNING invalid center" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &lengthX) < 0) {
			opserr << "WARNING invalid Length" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		//Simple_Line(int tag, double centerX, double lengthX);
		theHTEntity = new Simple_Line(HTEntityTag, centerX, lengthX);

	}
	else if (strcmp(EntityType, "Block") == 0 || strcmp(EntityType, "Block2D") == 0 || strcmp(EntityType, "Block2d") == 0)
	{
		double centerX = 0.0; double centerY = 0.0; double BreadthX = 0.0; double HeightY = 0.0;

		if (OPS_GetDoubleInput(&numdata, &centerX) < 0) {
			opserr << "WARNING invalid centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &centerY) <0) {
			opserr << "WARNING invalid centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata,&BreadthX) <0) {
			opserr << "WARNING invalid centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HeightY) <0) {
			opserr << "WARNING invalid centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Block(HTEntityTag, centerX, centerY, BreadthX, HeightY);

	}
	//Adding 2D entity:Isection
	else if (strcmp(EntityType, "Isection") == 0 || strcmp(EntityType, "Isection2D") == 0 || strcmp(EntityType, "Isection2d") == 0)
	{

		double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB;

		if (OPS_GetDoubleInput(&numdata, &HTI_centerX)<0 ) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_centerY) <0) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_BF) <0) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_HB)<0) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tw)<0) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tf)<0) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Isection(HTEntityTag, HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB - 2 * HTI_Tf);

	}

	else if (strcmp(EntityType, "ProtectedIsection") == 0 || strcmp(EntityType, "InsIsection2D") == 0 || strcmp(EntityType, "ProtectedIsection2d") == 0)
	{

		double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB, HTI_coat;

		if (OPS_GetDoubleInput(&numdata, &HTI_centerX) <0) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_centerY) <0) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_BF) <0) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_HB) <0) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tw) <0) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tf) <0) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_coat) <0) {
			opserr << "WARNING invalid HTI_coat" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		//Simple_IsecProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw, double HTI_coat)
		theHTEntity = new Simple_IsecProtected(HTEntityTag, HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB - 2 * HTI_Tf, HTI_coat);

	}

	//adding BRICK3D
	else if (strcmp(EntityType, "Brick") == 0 || strcmp(EntityType, "Brick3D") == 0 || strcmp(EntityType, "Brick3d") == 0)
	{

		//Simple_Brick(int tag, double CenterX, double CenterY, double CenterZ, double Breadth_X, double Height_Y,double Length_Z)
		double CenterX, CenterY, CenterZ, Breadth_X, Height_Y, Length_Z;

		if (OPS_GetDoubleInput(&numdata, &CenterX) <0) {
			opserr << "WARNING invalid centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &CenterY) <0) {
			opserr << "WARNING invalid centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &CenterZ) <0) {
			opserr << "WARNING invalid centerZ" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &Breadth_X) <0) {
			opserr << "WARNING invalid Breadth_X" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &Height_Y) <0) {
			opserr << "WARNING invalid Height_Y" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &Length_Z) <0) {
			opserr << "WARNING invalid Length_Z" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Brick(HTEntityTag, CenterX, CenterY, CenterZ, Breadth_X, Height_Y, Length_Z);

	}
	else if (strcmp(EntityType, "Isection3D") == 0 || strcmp(EntityType, "Isection3d") == 0)
	{

		//Simple_Isection3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
							  //double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw, double HTI_Len);
		double HTI_centerX, HTI_centerY, HTI_centerZ, HTI_BF, HTI_Tf, HTI_HB, HTI_Tw, HTI_Len;


		if (OPS_GetDoubleInput(&numdata, &HTI_centerX) <0) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_centerY) <0) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_centerZ) <0) {
			opserr << "WARNING invalid HTI_centerZ" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_BF) <0) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_HB) <0) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		if (OPS_GetDoubleInput(&numdata, &HTI_Tw) <0) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tf) <0) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}


		if (OPS_GetDoubleInput(&numdata, &HTI_Len) <0) {
			opserr << "WARNING invalid HTI_Len" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Isection3D(HTEntityTag, HTI_centerX, HTI_centerY, HTI_centerZ,
			HTI_BF, HTI_Tf, HTI_HB - 2 * HTI_Tf, HTI_Tw, HTI_Len);

	}
	//Adding 2D entity:Composite Section
	else if (strcmp(EntityType, "Composite2D") == 0 || strcmp(EntityType, "composite2D") == 0 || strcmp(EntityType, "Composite2d") == 0)
	{

		double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB, Slab_H, Slab_W;

		if (OPS_GetDoubleInput(&numdata, &HTI_centerX) <0) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_centerY) <0) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_BF) <0) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_HB) <0) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if (OPS_GetDoubleInput(&numdata, &HTI_Tw) <0) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		if (OPS_GetDoubleInput(&numdata, &HTI_Tf) <0) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		if (OPS_GetDoubleInput(&numdata, &Slab_W) <0) {
			opserr << "WARNING invalid Slab_w" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}
		if(OPS_GetDoubleInput(&numdata, &Slab_H) <0) {
			opserr << "WARNING invalid Slab_H" << endln;
			opserr << " for HeatTransfer entity: " << EntityType << endln;
			return -1;
		}

		//Simple_Composite2D(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf,
				 //double HTI_Tf, double HTI_Tw, double HTI_Hw,double slabW, double slabH)
		theHTEntity = new Simple_Composite2D(HTEntityTag, HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB - 2 * HTI_Tf, Slab_W, Slab_H);

	}


	if (theHTEntity != 0) {
		theHTModule->addHTEntity(*theHTEntity);
	}
	else
		opserr << "WARNING: TclHTModule fail to add HeatTransfer Material: " << EntityType << endln;

#ifdef _DEBUG
	opserr << "Entity successfully added" << endln;
#endif

	return 0;

}


//Adding Simple Mesh
int
OPS_addHTMesh()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	int HTMeshTag = 0;
	int HTEntityTag = 0;

	int HTMaterialTag = 0;
	int HTMaterialTag1 = 0;
	int PhaseChangeTag = 0;
	int PhaseChangeTag1 = 0;
	Vector* SectionLocs = new Vector(0);

	Simple_Entity* theHTEntity = 0;
	Simple_Mesh* theHTMesh = 0;
	HeatTransferMaterial* theHTMaterial = 0;
	HeatTransferMaterial* theHTMaterial1 = 0;
	Vector MeshCtrls = 0;


	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) {
		opserr << "WARNING::HTMesh requires at least 3 arguments." << "\n";
		return -1;
	}


	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &HTMeshTag) < 0) {
		opserr << "WARNING:: invalid mesh tag for defining simple mesh\n";
		return -1;
	}
	if (OPS_GetIntInput(&numdata, &HTEntityTag) < 0) {
		opserr << "WARNING:: invalid mesh tag for defining simple mesh: " << HTMeshTag << "\n";
		return -1;
	}
	if (OPS_GetIntInput(&numdata, &HTMaterialTag) < 0) {
		opserr << "WARNING:: invalid mesh tag for defining simple mesh: " << HTMeshTag << "\n";
		return -1;
	}

	const char* option = 0;
	if (OPS_GetNumRemainingInputArgs() > 0) {
		
		option = OPS_GetString();

		//if second material tag is detected
		if (strcmp(option, "-SecondMat") == 0 || strcmp(option, "-secondmaterial") == 0) {
			if (OPS_GetIntInput(&numdata, &HTMaterialTag1) < 0) {
				opserr << "WARNING:: invalid second material tag when defining simple mesh: " << HTMeshTag << "\n";
				return -1;
			}
		}
		else 
			OPS_ResetCurrentInputArg(-1);
		//if phasechange is defined
		option = OPS_GetString();
		if (strcmp(option, "-PhaseChange") == 0 || strcmp(option, "-phasechange") == 0) {
			if (OPS_GetIntInput(&numdata, &PhaseChangeTag) < 0) {
				opserr << "WARNING:: invalid phase change tag when defining simple mesh: " << HTMeshTag << "\n";
				return -1;
			}
			//if second material is defined
			if (HTMaterialTag1 != 0) {
				if (OPS_GetIntInput(&numdata, &PhaseChangeTag1) < 0) {
					opserr << "WARNING:: invalid Phasechange tag for second material when defining simple mesh\n";
					return -1;
				}
			}
		}
		else
			OPS_ResetCurrentInputArg(-1);

		//if section loc identified
		option = OPS_GetString();
		if (strcmp(option, "-SectionLoc") == 0 || strcmp(option, "-sectionLoc") == 0) {
			double Loc1, Loc2;
			if (OPS_GetDoubleInput(&numdata, &Loc1) < 0) {
				opserr << "WARNING:: invalid section location for defining simple mesh: " << HTMeshTag << "\n";
				return -1;
			}

			SectionLocs = new Vector(1);
			(*SectionLocs)(0) = Loc1;
			//if second material is defined
			
			if (OPS_GetDoubleInput(&numdata, &Loc2) < 0) {
				//opserr << "WARNING:: invalid section location for defining simple mesh: " << HTMeshTag << "\n";
				//return -1;
				OPS_ResetCurrentInputArg(-1);
			}
			else {
				SectionLocs->resize(2);
				(*SectionLocs)(0) = Loc1;
				(*SectionLocs)(1) = Loc2;
			}
			
			
		}
		else
			OPS_ResetCurrentInputArg(-1);
	}

	//if meshctrls tag is detected
	option = OPS_GetString();
	if (strcmp(option, "-MeshCtrls") == 0 ||strcmp(option, "-meshctrls") == 0) {
		int numRArgs = 0;
		if (OPS_GetNumRemainingInputArgs() > 0) {
			numRArgs = OPS_GetNumRemainingInputArgs();
			MeshCtrls.resize(numRArgs);
		}
		else {
			opserr << "WARNING:: no meshCtrl defined for simple mesh: " << HTMeshTag << "\n";
			return -1;
		}
		//-----for geting uncertain number of double data.
		double data;
		
		for (int i = 0; i < numRArgs; i++) {
				OPS_GetDoubleInput(&numdata,&data);
				MeshCtrls(i) = data;
		}
		
#ifdef _DEBUG
		opserr << "OriginLocs " << *SectionLocs << " MeshCtrls " << MeshCtrls << endln;
#endif
		// for geting uncertain number of doubel values  
	}
	else {
		opserr << "WARNING:TclHTModule- MeshCtrls not found for Mesh " << HTMeshTag << endln;
		return -1;
	}

	theHTEntity = theHTModule->getHTEntity(HTEntityTag);
	if (theHTEntity == 0) {
		opserr << "WARNING:: TclHTModule failed to get the requested entity when defining simple mesh: " << HTMeshTag << "\n";
		return -1;
	}

	theHTMaterial = theHTModule->getHTMaterial(HTMaterialTag);
	if (theHTMaterial == 0) {
		opserr << "WARNING:: TclHTModule failed to get the requested HT material when defining Simple Mesh: " << HTMeshTag << "\n";
		return -1;
	}

	if (HTMaterialTag1 != 0) {
		theHTMaterial1 = theHTModule->getHTMaterial(HTMaterialTag1);
		if (theHTEntity == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested HT material when defining Simple Mesh: " << HTMeshTag << "\n";
			return -1;
		}
	}

	ID PHaseIDs = 0;
	if (PhaseChangeTag1 != 0) {
		PHaseIDs.resize(2);
		PHaseIDs(0) = PhaseChangeTag;
		PHaseIDs(1) = PhaseChangeTag1;
	}
	else
	{
		if (PhaseChangeTag != 0) {
			PHaseIDs.resize(1);
			PHaseIDs(0) = PhaseChangeTag;
		}

	}

	if (HTMaterialTag1 != 0) {
		theHTMesh = new Simple_Mesh(HTMeshTag, theHTEntity, theHTDomain, theHTMaterial, MeshCtrls, theHTMaterial1);
	}
	else {
		theHTMesh = new Simple_Mesh(HTMeshTag, theHTEntity, theHTDomain, theHTMaterial, MeshCtrls);
	}

	if (theHTMesh != 0) {
		if (PHaseIDs != 0) {
			theHTMesh->SetEleParameters(PHaseIDs);
		}
		if (SectionLocs->Size() != 0)
			theHTMesh->SetOriginLocs(*SectionLocs);

		theHTModule->addHTMesh(theHTMesh);
	}
	else
		opserr << "WARNING: TclHTModule fail to add Simple Mesh: " << HTMeshTag << endln;


	return 0;

}



//Opensees Command_HTRefineMesh
int
OPS_HTRefineMesh()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	int HTEntityTag = 0;
	Simple_Entity* theHTEntity = 0;
	int SeedTag = 0;
	int SeedTag1 = 0;
	Vector MeshVec = 0;
	const char* option = OPS_GetString();
	int numdata = 1;

	//if HTEntity tag is detected
	if (strcmp(option, "-HTEntity") == 0 || strcmp(option, "-Entity") == 0 || strcmp(option, "HTEntity") == 0) {

		if (OPS_GetIntInput(&numdata, &HTEntityTag) <0) {
			opserr << "WARNING:: invalid HT Entity tag for refining the mesh\n";
			return -1;
		}

	}
	option = OPS_GetString();
	if (strcmp(option, "-SeedTag") == 0 || strcmp(option, "-seedTag") == 0 || strcmp(option, "SeedTag") == 0) {

		if (OPS_GetIntInput(&numdata, &SeedTag)<0) {
			opserr << "WARNING:: invalid seed tag for refining the mesh: " << HTEntityTag << "\n";
			return -1;
		}
	

		if (OPS_GetIntInput(&numdata, &SeedTag1)<0) {
			OPS_ResetCurrentInputArg(-1);
			SeedTag1 = 0;
		}
	

		//obtaining second seed tag if necessary

	}
	//Seed tag
	option = OPS_GetString();
	if (strcmp(option, "-space") == 0 || strcmp(option, "-Space") == 0 || strcmp(option, "Space") == 0) {
		int numRargs = 0;
		if (OPS_GetNumRemainingInputArgs() > 0) {
			numRargs = OPS_GetNumRemainingInputArgs();
			MeshVec.resize(numRargs);
		}
		else {
			opserr << "WARNING:: no Seed info defined for refining mesh: " << HTEntityTag << "\n";
			return -1;
		}

		//-----for geting uncertain number of double data.

		double data;

		for (int i = 0; i < numRargs; i++) {
			OPS_GetDoubleInput(&numdata, &data);
			MeshVec(i) = data;
		}

		// for geting uncertain number of doubel values

#ifdef _DEBUG
		opserr << "HTEntity " << HTEntityTag << " has a refined mesh as" << MeshVec << endln;
#endif
	}

	//space info obtained
	if (SeedTag != 0) {
		if (theHTEntity->RefineSeeds(SeedTag, MeshVec) != 0) {
			opserr << "WARNING::the Entity " << HTEntityTag << "failed to refine Seeds" << HTEntityTag << endln;
		}
	}
	if (SeedTag1 != 0) {
		if (theHTEntity->RefineSeeds(SeedTag1, MeshVec) != 0) {
			opserr << "WARNING::the Entity " << HTEntityTag << "failed to refine Seeds " << SeedTag1 << endln;
		}
	}

	////////////////////////////////////////////
}

//MeshAll
int
OPS_HTMeshAll()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	Simple_Mesh* theHTMesh;
	Simple_MeshIter& theHTMeshes = theHTModule->getHTMeshs();
	while ((theHTMesh = theHTMeshes()) != 0) {
		theHTMesh->GeneratingNodes();
		theHTMesh->GeneratingEles();
	}
#ifdef _DEBUG
	opserr << "Heat Transfer module has completed the mesh"<< endln;
#endif
	return 0;

}


//SetInitialT
int
OPS_SetInitialT()
{

	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	double initialT;
	int numdata = 1;
	if (OPS_GetDoubleInput(&numdata, &initialT) <0) {
		opserr << "WARNING invalid initial Temperature" << endln;
		opserr << " for TclHeatTransferModule: setInitialT " << endln;
		return -1;
	}

	theHTDomain->setInitial(initialT);

	return 0;

}


//HTConstants
int
OPS_addHTConstants()
{

	//Recieving heat transfer coefficients: convective coefficient, ambient temperature, absorption radiative coefficient,
	 //                                     steffan-B constant ,fire radiative cofficient, irradiation
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
		return -1;
	}

	HTConstants* theHTConstants = 0;
	int HTConstantTag = 0;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &HTConstantTag) <0) {
		opserr << "WARNING:: invalid constants tag for defining HT constants\n";
		return -1;
	}

	int numArgs = OPS_GetNumRemainingInputArgs();
	Vector constants(5);
	if ((numArgs != 4) && (numArgs!= 5)) {
		opserr << "WARNING:: Non-sufficent or oversized data for defining heat transfer contants: " << HTConstantTag << "\n";
		return -1;
	}
	//-----for geting uncertain number of double data.
	double data;

	for (int i = 0; i < numArgs; i++) {
		if (OPS_GetDoubleInput(&numData, &data) < 0) {
			opserr << "HTConstants command is unable to identify an input for constants: " << HTConstantTag << "/n";
			return -1;
		}
			
		constants(i) = data;
	}


	// for geting uncertain number of doubel values 
	if (numArgs == 4) {
		constants(4) = 5.67e-8 * constants(1) * constants(1) * constants(1) * constants(1);
#ifdef _DEBUG
		opserr << "TclHeatTransfer:: AddHTconstants: " << constants << endln;
#endif
		//calculating radiation from the ambient air;
	}

	theHTConstants = new HTConstants(HTConstantTag, constants);

	if (theHTConstants != 0) {
		theHTModule->addHTConstants(theHTConstants);
	}
	else
		opserr << "WARNING: TclHTModule fail to add HTConstants: " << HTConstantTag << endln;

	return 0;

}



//HTEleSet
int
OPS_HTEleSet()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTEleSet\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTNodeSet\n";
		return -1;
	}

	HTEleSet* theHTEleSet = 0;
	int HTEleSetTag = 0;
	int HTEntityTag = 0;
	int FaceID = 0;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &HTEleSetTag) <0) {
		opserr << "WARNING:: invalid EleSet tag for slelecting elements\n";
		return -1;
	}


	//if HTEntity tag is detected
	const char* option = OPS_GetString();
	if (strcmp(option, "-HTEntity") == 0 || strcmp(option, "-Entity") == 0 || strcmp(option, "HTEntity") == 0) {

		if (OPS_GetIntInput(&numData, &HTEntityTag) <0) {
			opserr << "WARNING:: invalid HT Entity tag for defining HTEleSet: " << HTEleSetTag << "\n";
			return -1;
		}
	}
	else {
		opserr << "WARNING TclHTCommand HTEleSet wants tag -HTEntity, HTEntity, -Entity " << endln;
	}


	Simple_Entity* theHTEntity = theHTModule->getHTEntity(HTEntityTag);
	if (theHTEntity == 0) {
		opserr << "WARNING:: TclHTModule failed to get the requested entity when defining simple mesh: " << HTEleSetTag << "\n";
		return -1;
	}

	int EntityMeshTag = theHTEntity->getMeshTag();
	Simple_Mesh* theHTMesh = theHTModule->getHTMesh(EntityMeshTag);


	ID EleRange(0);
	if (OPS_GetNumRemainingInputArgs() > 0) {
		const char* option = OPS_GetString();
		//if face tag is detected
		if (strcmp(option, "-face") == 0 || strcmp(option, "-Face") == 0 || strcmp(option, "face") == 0)
		{
			if (OPS_GetIntInput(&numData, &FaceID)<0) {
				opserr << "WARNING:: invalid face tag for defining HTEleSet: " << HTEleSetTag << "\n";
				return -1;
			}
			int eleFaceID;
			theHTMesh->SelectingElesbyFace(EleRange, FaceID, eleFaceID);
		}
		else if (strcmp(option, "-NodeSet") == 0 || strcmp(option, "-nodest") == 0 || strcmp(option, "NodeSet") == 0)
		{
			int NodeSetID;
			if (OPS_GetIntInput(&numData, &NodeSetID) <0) {
				opserr << "WARNING:: invalid NodesetID for defining HTEleSet: " << HTEleSetTag << "\n";
				return -1;
			}

			HTNodeSet* theHTNodeSet = theHTModule->getHTNodeSet(NodeSetID);

			if (theHTNodeSet == 0) {
				opserr << "WARNING:: invalid NodesetID for defining HTEleSet: " << HTEleSetTag << "\n";
				return -1;
			}

			ID NodesRange = theHTNodeSet->getNodeID();
			if (OPS_GetNumRemainingInputArgs() > 0) {
				option = OPS_GetString();
				if (strcmp(option, "-face") == 0 || strcmp(option, "-Face") == 0 || strcmp(option, "face") == 0) {

					if (OPS_GetIntInput(&numData, &FaceID) <0) {
						opserr << "WARNING:: invalid face tag for defining HTEleSet: " << HTEleSetTag << "\n";
						return -1;
					}
				}
			}
			
			theHTMesh->SelectingEles(EleRange, NodesRange, FaceID);

		}

	}
	

	theHTEleSet = new HTEleSet(HTEleSetTag);

	if (EleRange != 0) {
		theHTEleSet->addEleID(EleRange);
	}

	if (theHTEleSet != 0) {
		theHTModule->addHTEleSet(theHTEleSet);
	}
	else
		opserr << "WARNING: TclHTModule failed to add HTEleSet: " << HTEleSetTag << endln;

	return 0;

}

//--------------------
//HTNodeSet
int
OPS_HTNodeSet()
{
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTNodeSet\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTNodeSet\n";
		return -1;
	}

	HTNodeSet* theHTNodeSet = 0;
	Simple_Entity* theHTEntity = 0;
	int HTNodeSetTag = 0;
	int HTEntityTag = 0;
	int FaceID = 0;
	int numData = 1;
	ID NodeRange(0);

	if (OPS_GetNumRemainingInputArgs() < 3)
		opserr << "WARNING:no arguments for HTNodeSet" << endln;

	if (OPS_GetIntInput(&numData, &HTNodeSetTag) <0) {
		opserr << "WARNING:: insufficient information for defining HTNodeSet\n";
		return -1;
	}

	const char* option = OPS_GetString();
	//if HTEntity tag is detected, then we definitely expect faceID
	if (strcmp(option, "-HTEntity") == 0 || strcmp(option, "-Entity") == 0 || strcmp(option, "HTEntity") == 0) {

		if (OPS_GetIntInput(&numData, &HTEntityTag) <0) {
			opserr << "WARNING:: invalid HT Entity tag for defining HTNodeSet: " << HTNodeSetTag << "\n";
			return -1;
		}

		theHTEntity = theHTModule->getHTEntity(HTEntityTag);
		if (theHTEntity == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested entity when defining simple mesh: " << HTNodeSetTag << "\n";
			return -1;
		}

		int EntityMeshTag = theHTEntity->getMeshTag();
		Simple_Mesh* theHTMesh = theHTModule->getHTMesh(EntityMeshTag);

		//to obtain faceTag
		if (OPS_GetNumRemainingInputArgs() > 0) {
			option = OPS_GetString();
		}
		else
			opserr << "WARNING:: insufficient information for defining HTNodeSet: " << HTNodeSetTag << "\n";

		if (strcmp(option, "-face") == 0 || strcmp(option, "-Face") == 0 || strcmp(option, "face") == 0) {
			if (OPS_GetIntInput(&numData, &FaceID) < 0) {
				opserr << "WARNING:: invalid face tag for defining HTEleSet: " << HTNodeSetTag << "\n";
				return -1;
			}
			theHTMesh->SelectingNodesbyFace(NodeRange, FaceID);
		}
		else {
			OPS_ResetCurrentInputArg(-1);
		}

	}
	else {
		OPS_ResetCurrentInputArg(-1);
	}
	//HTEntityTag obtained


	
	//end of HTEntity tag 

  //search node in the range of xloc
  //check whether it is the end of arguments
	if (OPS_GetNumRemainingInputArgs() > 0) {
		option = OPS_GetString();
		if (strcmp(option, "-Locx") == 0 || strcmp(option, "Locx") == 0 || strcmp(option, "locx") == 0 || strcmp(option, "-locx") == 0) {
			
			double xlocLB, xlocUB;

			if (OPS_GetDoubleInput(&numData, &xlocLB) <0) {
				opserr << "WARNING invalid xloc" << endln;
				opserr << " for adding node to HTNodeSet: " << HTNodeSetTag << endln;
				return -1;
			}

			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &xlocUB) < 0) {
					xlocUB = xlocLB;
					OPS_ResetCurrentInputArg(-1);
				}
			}
			else
				xlocUB = xlocLB;


			if (xlocUB < xlocLB) {
				opserr << "WARNING::HTNodeSet "<< HTNodeSetTag<<" should receive xlocUB " << xlocUB << " greater than xlocLb " << xlocLB;
				return -1;
			}
			//opserr << xlocUB << " " << xlocLB;
			theHTDomain->SelectingNodes(NodeRange, 0, xlocLB, xlocUB);
			// for geting uncertain number of doubel values

		}
		else
			OPS_ResetCurrentInputArg(-1);
	}
	

	//search node in the range of yloc
	if (OPS_GetNumRemainingInputArgs() > 0) {
		option = OPS_GetString();
		if (strcmp(option, "-Locy") == 0 || strcmp(option, "Locy") == 0 || strcmp(option, "locy") == 0) {
			double ylocLB, ylocUB;

			if (OPS_GetDoubleInput(&numData, &ylocLB) < 0) {
				opserr << "WARNING invalid yloc" << endln;
				opserr << " for adding node to HTNodeSet: " << HTNodeSetTag << endln;
				return -1;
			}

			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &ylocUB) < 0) {
					ylocUB = ylocLB;
					OPS_ResetCurrentInputArg(-1);
				}
					
			}
			else
				ylocUB = ylocLB;


			if (ylocUB < ylocLB) {
				opserr << "WARNING::HTNodeSet " << HTNodeSetTag << " should receive ylocUB " << ylocUB << " greater than ylocLb " << ylocLB;
				return -1;
			}

			theHTDomain->SelectingNodes(NodeRange, 1, ylocLB, ylocUB);
			// for geting uncertain number of doubel values

		}
		else
			OPS_ResetCurrentInputArg(-1);
	}
	//search node in the range of zloc
	if (OPS_GetNumRemainingInputArgs() > 0) {
		option = OPS_GetString();
		if (strcmp(option, "-Locz") == 0 || strcmp(option, "Locz") == 0 || strcmp(option, "locz") == 0) {
			double zlocLB, zlocUB;

			if (OPS_GetDoubleInput(&numData, &zlocLB) < 0) {
				opserr << "WARNING invalid zloc" << endln;
				opserr << " for adding node to HTNodeSet: " << HTNodeSetTag << endln;
				return -1;
			}

			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDoubleInput(&numData, &zlocUB) < 0) {
					zlocUB = zlocLB;
					OPS_ResetCurrentInputArg(-1);
				}
				
			}
			else
				zlocUB = zlocLB;


			if (zlocUB < zlocLB) {
				opserr << "WARNING::HTNodeSet " << HTNodeSetTag << " should receive zlocUB " << zlocUB << " greater than zlocLb " << zlocLB;
				return -1;
			}

			theHTDomain->SelectingNodes(NodeRange, 2, zlocLB, zlocUB);
			// for geting uncertain number of doubel values

		}
		else
			OPS_ResetCurrentInputArg(-1);
	}


	theHTNodeSet = new HTNodeSet(HTNodeSetTag);

	if (NodeRange != 0) {
		theHTNodeSet->addNodeID(NodeRange);
	}
	else {
		opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << HTNodeSetTag << endln;
	}

	if (theHTNodeSet != 0)
		theHTModule->addHTNodeSet(theHTNodeSet);
	else
		opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << HTNodeSetTag << endln;

#ifdef _DEBUG
	opserr << "NodeSet " << HTNodeSetTag << " grouped nodes:" << NodeRange << endln;
#endif

	int numNodes = NodeRange.Size();
	int* data = new int[numNodes];
	for (int i = 0; i < numNodes; i++) {
		data[i] = (int)(NodeRange(i));
	}
	if (OPS_SetIntOutput(&numNodes, data) < 0) {
		opserr << "WARNING HeatTransferModule failed to set outputs for HTNodeSet\n";
		delete[] data;
		return -1;
	}


	return 0;

}

//HTPattern
int
OPS_addHTPattern()
{
	if (OPS_GetNumRemainingInputArgs()< 2) {
		opserr << "WARNING invalid command-want:HTPattern type\n";
		opserr << " valid types: AmbientBC, fireAction(fire) \n";
		return -1;
	}

	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTPattern\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - HTPattern\n";
		return -1;
	}

	BoundaryPattern* HTPattern = 0;
	int HTPatternTag = 0;
	int numData = 1;
	const char* option = OPS_GetString();
	if (OPS_GetIntInput(&numData, &HTPatternTag) <0) {
		opserr << "WARNING:: invalid pattern tag for defining HTPattern: " << option<< "\n";
		return -1;
	}

	int  commandEndMarker = 0;
	//if face tag is detected
	if (strcmp(option, "AmbientBC") == 0 || strcmp(option, "ambientBC") == 0 || strcmp(option, "ambient") == 0) {

		HTPattern = new BoundaryPattern(HTPatternTag);

	}
	//if fireAction tag is detected here
	else if (strcmp(option, "fireAction") == 0 || strcmp(option, "FireAction") == 0 || strcmp(option, "fire") == 0) {

		int FireModelID = 0;
		HTPattern = new FireImposedPattern(HTPatternTag);
		if(OPS_GetNumRemainingInputArgs()<0)
			opserr<< "WARNING:: insufficient information for defining HTPattern: " << option << "\n";
		else 
			option = OPS_GetString();

		if (strcmp(option, "model") == 0 || strcmp(option, "Model") == 0 || strcmp(option, "-model") == 0) {
			if (OPS_GetIntInput(&numData, &FireModelID) <0) {
				opserr << "WARNING:: invalid fire model tag for defining HTPattern: " << HTPatternTag << "\n";
				return -1;
			}
		}

		FireModel* theFireModel = theHTModule->getFireModel(FireModelID);
		if (theFireModel == 0) {
			opserr << "WARNING:: no fire model is defined for heat flux BC: " << "\n";
			return -1;
		}
		((FireImposedPattern*)HTPattern)->setFireModel(theFireModel);

	}
	//END OF FIRE ACTION PATTERN
	else {
		opserr << "WARNING unknown pattern type" << option;
		opserr << "- want : HTPattern patternType " << HTPatternTag;
		opserr << " valid types: AmbientBC, fireAction(fire) \n";
		return -1;
	}

	//Add boundaryPattern into the heat transfer domain
	if (theHTDomain->addBoundaryPattern(HTPattern) == false) {
		opserr << "WARNING could not add boundary pattern to the heat transfer domain" << HTPatternTag;
		delete HTPattern;
		return -1;
	}

	theHTPattern = HTPattern;

	/*
	// use TCL_Eval to evaluate the list of load and single point constraint commands
	if (Tcl_Eval(interp, argv[argc - 1]) != TCL_OK) {
		opserr << "WARNING - error reading load pattern information in { } ";
		return -1;
	}
	*/
	return 0;
}


//Add Fire Model
int
OPS_addFireModel()
{

	// checking the TclHTModule exists or not
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - HTPattern\n";
		return -1;
	}
	// checking the HTDomain exists or not
	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain -  HTPattern\n";
		return -1;
	}

	//create a pointer to base class
	FireModel* theFireModel = 0;
	int FireModelTag = 0;
	int numData = 1;
	//get the fireModel tag;
	if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "WARNING:: invalid fireModel input\n";
		return -1;
	}

	const char* option = OPS_GetString();

	if (OPS_GetIntInput(&numData, &FireModelTag)<0) {
		opserr << "WARNING:: invalid fireModel Tag\n";
		return -1;
	}
	  //standard fire curve;
	if (strcmp(option, "standard") == 0 || strcmp(option, "Standard") == 0) {
		theFireModel = new NorminalFireEC1(FireModelTag, 1);

	}
	//hydroCarbon fire curve;
	else if (strcmp(option, "hydroCarbon") == 0 || strcmp(option, "HydroCarbon") == 0) {
		theFireModel = new NorminalFireEC1(FireModelTag, 3); // hydrocarbon fire tag is 3;
	}
	else if (strcmp(option, "ASTM") == 0 || strcmp(option, "ASTME119") == 0) {
		theFireModel = new NorminalFireEC1(FireModelTag, 4); // hydrocarbon fire tag is 3;
	}
	else if (strcmp(option, "Exponent") == 0 || strcmp(option, "Exp") == 0) {
		theFireModel = new NorminalFireEC1(FireModelTag, 5); // hydrocarbon fire tag is 3;
	}
	else if (strcmp(option, "UserDefined") == 0 || strcmp(option, "userDefined") == 0) {
		if (OPS_GetNumRemainingInputArgs() < 2) {
			opserr << "WARNING:: insufficient input for userdefined fireModel\n";
			return -1;
		}
		option = OPS_GetString();
		const char* filename = 0;
		int dataTypeTag = 1;
		if (strcmp(option, "file") == 0 || strcmp(option, "-file" )|| strcmp(option, "File") == 0) {
			filename = OPS_GetString();
			theFireModel = new UserDefinedFire(FireModelTag, filename, 1);
		}
		if (OPS_GetNumRemainingInputArgs() > 0) {
			const char* dataType = OPS_GetString();
			if (strcmp(dataType, "-type") == 0 || strcmp(dataType, "dataType") == 0 || strcmp(dataType, "type") == 0) {
				
				if (OPS_GetIntInput(&numData, &dataTypeTag) <0) {
					opserr << "WARNING:: invalid Data Type Tag for user defined fire model " << dataType << "\n";
					return -1;
				}
			}

		}



		theFireModel = new UserDefinedFire(FireModelTag, filename, dataTypeTag);
	}
	//Paramtetric fire
	else if (strcmp(option, "parametric") == 0 || strcmp(option, "Parametric") == 0) {
		double thi = 0; double avent = 0; double hvent = 0; 
		double atotal = 0; double afire = 0; double qfire = 0; double Tlim = 0;

		if (OPS_GetNumRemainingInputArgs()>6) {
			if (OPS_GetDoubleInput(&numData, &thi) <0) {
				opserr << "WARNING invalid thermal inertia of the compartment boundaries" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &avent) <0) {
				opserr << "WARNING invalid total area of vertial openings on walls" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &hvent) <0) {
				opserr << "WARNING invalid weighted average of window heights on walls" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &atotal) <0) {
				opserr << "WARNING invalid total area of the compartment(including walls)" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &afire) <0) {
				opserr << "WARNING invalid area of the floor with fire" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &qfire) <0) {
				opserr << "WARNING invalid total design fire" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
			if (OPS_GetDoubleInput(&numData, &Tlim) <0) {
				opserr << "WARNING invalid time levels corresponds to different fire growth rate" << endln;
				opserr << " for HeatTransfer fire model: " << FireModelTag << endln;
				return -1;
			}
		}
		else
			opserr << "WARNING:: Defining Parametric fire: " << FireModelTag << " recieved insufficient arguments" << "\n";

		theFireModel = new ParametricFireEC1(FireModelTag, thi, avent, hvent, atotal, afire, qfire, Tlim);
	}
	//localised fire curve;
	else if (strcmp(option, "localised") == 0 || strcmp(option, "Localised") == 0) {

		double crd1 = 0.0; double crd2 = 0.0; double crd3 = 0.0;
		double D = 0; double Q = 0; double H = 0; int lineTag = 0;

		if (OPS_GetNumRemainingInputArgs()<3)
		{
			opserr << "WARNING insufficient aguments" << endln;
			opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
			return -1;
		}

		//Add a tag for location of origin;
		option = OPS_GetString();
		if (strcmp(option, "-origin") == 0 || strcmp(option, "origin") == 0) {

			if (OPS_GetDoubleInput(&numData, &crd1) <0) {
				opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
			
			if (OPS_GetDoubleInput(&numData, &crd2) <0) {
				opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}

			if (OPS_GetDoubleInput(&numData, &crd3)<0) {
				//it it possible not to have a z loc for localised fire definiton
				opserr << "z loc is set as default 0" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				OPS_ResetCurrentInputArg(-1);
				crd3 = 0.0;
			}
		
		}
		else {
			OPS_ResetCurrentInputArg(-1);
		}

		//end of fire origin, waiting for firePars;
		option = OPS_GetString(); //in case no origin in input
		if (strcmp(option, "-firePars") == 0 || strcmp(option, "firePars") == 0) {

			if (OPS_GetDoubleInput(&numData, &D) <0) {
				opserr << "WARNING invalid diameter of the fire source" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
			
			if (OPS_GetDoubleInput(&numData, &Q) <0) {
				opserr << "WARNING invalid rate of heat release" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		
			if (OPS_GetDoubleInput(&numData, &H) <0) {
				opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
			//detect argument for linetag;
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&numData, &lineTag) <0) {
					opserr << "WARNING invalid central line tag " << endln;
					opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
					return -1;
				}
			}
			else {
				//opserr << "Central line tag for localised fire " << FireModelTag << "is set as default:3" << endln;
				lineTag = 3;
			}
		}
		else {
			opserr << "WARNING:: Defining Localised fire " << FireModelTag
				<< " expects tag:-firePars or firePars" << "\n";
		}

		theFireModel = new LocalizedFireEC1(FireModelTag, crd1, crd2, crd3, D, Q, H, lineTag);
	}
	//Travelling fire curve;
	else if (strcmp(option, "travelling") == 0 || strcmp(option, "Travelling") == 0 || strcmp(option, "NaturalFire") == 0) {

		double crd1 = 0.0; double crd2 = 0.0; double crd3 = 0.0;
		double D = 1.0; double Q = 1e6; double H = 3.0; double Ts = 293.15;  int lineTag = 2;
		PathTimeSeriesThermal* theSeries = 0;


		if (OPS_GetNumRemainingInputArgs() > 0)
		{
			//end of fire origin, waiting for firePars;
			option = OPS_GetString(); //in case no origin in input
			if (strcmp(option, "firePars") == 0 || strcmp(option, "-firePars") == 0) {
				option = OPS_GetString();
				if (strcmp(option, "source") == 0 || strcmp(option, "-source") == 0 || strcmp(option, "-file") == 0 || strcmp(option, "file") == 0) {
					if (OPS_GetNumRemainingInputArgs() < 1)
					{
						opserr << "WARNING insufficient aguments" << endln;
						opserr << " for HeatTransfer Travelling fire model: " << FireModelTag << endln;
						return -1;
					}
					const char* fileName = OPS_GetString();
					theSeries = new PathTimeSeriesThermal(1, fileName, 8, false);
					//8 parameters in total�� firelocs(3), Q, D, 
				}

				//To get line tag if available
				if (OPS_GetNumRemainingInputArgs() > 0) {
					if (OPS_GetIntInput(&numData, &lineTag) < 0) {
						opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
						opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
						lineTag =2;
					}

				}

			}
			else if (strcmp(option, "H") == 0 || strcmp(option, "-H") == 0) {
				if (OPS_GetDoubleInput(&numData, &H) < 0) {
					opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
					opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
					return -1;
				}

			}
		}

		if (theSeries != 0) {
			theFireModel = new NaturalFire(FireModelTag, D, Q, H, lineTag, Ts, theSeries);
		}
		else {
			theFireModel = new NaturalFire(FireModelTag, D, Q, H, lineTag, Ts);
		}
	}
	//else ----------
	else if (strcmp(option, "idealised") == 0 || strcmp(option, "Idealised") == 0) {

		double crd1 = 0.0; double crd2 = 0.0; double crd3 = 0.0;
		double Q = 0; double D1 = 0; double D2 = 0; double K1 = 0; double K2 = 0; int lineTag = 0;
		int typeTag = 1;

		if (OPS_GetNumRemainingInputArgs() < 3)
		{
			opserr << "WARNING insufficient aguments" << endln;
			opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
			return -1;
		}

		//Add a tag for location of origin;
		option = OPS_GetString();
		if (strcmp(option, "-origin") == 0 || strcmp(option, "origin") == 0) {
			if (OPS_GetDoubleInput(&numData, &crd1) <0) {
				opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		
			if (OPS_GetDoubleInput(&numData, &crd2) <0) {
				opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}

			
			if (OPS_GetDoubleInput(&numData, &crd3) < 0) {
				//it it possible not to have a z loc for localised fire definiton
				opserr << "WARNING invalid z axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				OPS_ResetCurrentInputArg(-1);
				crd3 = 0.0;
			}

			
		}
		//end of fire origin, waiting for firePars;
		option = OPS_GetString();
		if (strcmp(option, "-Q") == 0 || strcmp(option, "Q") == 0 || strcmp(option, "HRR") == 0) {

			if (OPS_GetDoubleInput(&numData, &Q) <0) {
				opserr << "WARNING invalid rate of heat release" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		}

		option = OPS_GetString();
		if (strcmp(option, "-Linear") == 0 || strcmp(option, "Linear") == 0 || strcmp(option, "-linear") == 0) {
			typeTag = 1;

			if (OPS_GetDoubleInput(&numData, &D1) <0) {
				opserr << "WARNING invalid range of pleatou" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}

			if (OPS_GetDoubleInput(&numData, &D2) <0) {
				opserr << "WARNING invalid length of decay" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		}
		else if (strcmp(option, "-quadratic") == 0 || strcmp(option, "Quadratic") == 0 || strcmp(option, "-quadratic") == 0) {
			typeTag = 2;
		
			if (OPS_GetDoubleInput(&numData, &D1) <0) {
				opserr << "WARNING invalid distance D1" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
			
			if (OPS_GetDoubleInput(&numData, &D2) <0) {
				opserr << "WARNING invalid distance D2" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		
			if (OPS_GetDoubleInput(&numData, &K1) <0) {
				opserr << "WARNING invalid distance ratio K1" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		}
		else if (strcmp(option, "-exp") == 0 || strcmp(option, "exponential") == 0 || strcmp(option, "-Exp") == 0) {
			typeTag = 3;

			if (OPS_GetDoubleInput(&numData, &D1) <0) {
				opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}

			if (OPS_GetDoubleInput(&numData, &D2) <0) {
				opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}

			if (OPS_GetDoubleInput(&numData, &K1) <0) {
				opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
				opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
				return -1;
			}
		}

		//detect argument for linetag;
		if (OPS_GetNumRemainingInputArgs()> 0) {
			option = OPS_GetString();
			if (strcmp(option, "-centreLine") == 0 || strcmp(option, "CentreLine") == 0 || strcmp(option, "-centreline") == 0) {
				if (OPS_GetIntInput(&numData, &lineTag) <0) {
					opserr << "WARNING invalid central line tag " << endln;
					opserr << " for HeatTransfer localised fire model: " << FireModelTag << endln;
					return -1;
				}
			}
		}
		else {
			opserr << "WARNING:: Defining idealised fire " << FireModelTag
				<< " expects tag:-firePars or firePars" << "\n";
		}
		if (typeTag == 1) {
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3, Q, D1, D2, lineTag);
		}
		else if (typeTag == 2) {
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3, Q, D1, D2, K1, K2, lineTag);
		}
		else if (typeTag == 3) {
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3, Q, D1, D2, K1, lineTag);
		}
		//Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3, double Q, double D1, double D2, double K1, double K2, int centerLineTag = 3);

		//Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3,double Q, double D1, double D2, double factor, int centerLineTag = 3);

	}
	else {
		opserr << "WARNING unknown fire model type" << option;
		opserr << "- for fire model " << FireModelTag;
		opserr << " valid types: standard, hydroCarbon,parametric, localised.. \n";
		return -1;
	}

	if (theFireModel != 0) {
		theHTModule->addFireModel(theFireModel);
#ifdef _DEBUG
		OPS_Stream* output = &opserr;
		theFireModel->Print(*output);
#endif
	}
	else
	{
		opserr << "WARNING: TclHTModule fail to add FireModel: " << option << endln;
	}

	return 0;

}


//Add HeatFluxBC
int
OPS_addHeatFluxBC()
{

	//HeatFluxBC -eleset 1 -facetag 3 -type -ConvecAndRad -HTConstants 1;

	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
		return -1;
	}

	//if (theTclHTBoudary == 0) {
	 // opserr << "WARNING no active HeatTransfer Simple_Boundary - addHeatFluxBC\n";
	  //return -1;
	//}
	//int PatternTag = 0;
	int HTEleSetID = 0;
	int HTEntityTag = 0;
	Simple_Entity* theHTEntity = 0;
	HTEleSet* theHTEleSet = 0;
	int FaceID = 0;
	ID EntityFaceID = 0;
	ID ElesRange = 0;
	int HeatFluxTypeTag = 0;
	int HTconstantsID = 0;
	int FireType = 0;

	int dataNum = 1;
	if (OPS_GetNumRemainingInputArgs() < 3) {
		opserr << "WARNING:: insufficient input for HeatfluxBC\n";
		return -1;
	}
	/*
	if (OPS_GetIntInput(&dataNum, &PatternTag) < 0)
	{
		opserr << "WARNING:: HeatFluxBC can not identify the assigned pattern tag\n";
		return -1;
	}
	*/
	const char* option = OPS_GetString();
	//if HTEntity tag is detected
	if (strcmp(option, "-HTEntity") == 0 || strcmp(option, "-Entity") == 0 || strcmp(option, "HTEntity") == 0) {


		if (OPS_GetIntInput(&dataNum, &HTEntityTag) <0) {
			opserr << "WARNING:: invalid HT Entity tag for defining HTEleSet: " << option << "\n";
			return -1;
		}
	

		theHTEntity = theHTModule->getHTEntity(HTEntityTag);
		if (theHTEntity == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested entity when defining simple mesh: " << option << "\n";
			return -1;
		}

		option = OPS_GetString();
		//if face tag is detected
		if (strcmp(option, "-face") == 0 || strcmp(option, "-Face") == 0 || strcmp(option, "-faceTag") == 0) {


			//-----for geting uncertain number of double data.
			int numRargs = OPS_GetNumRemainingInputArgs();
			int count = 0; int ArgEnd = 0;
			int data;
			while (count < numRargs && ArgEnd == 0) {
				if (OPS_GetIntInput(&dataNum, &data) < 0) {
					ArgEnd = count;
					OPS_ResetCurrentInputArg(-1);
				}
				else
					count++;
			}
			//~ detecting the remianing number of input
			OPS_ResetCurrentInputArg(-count);
			EntityFaceID.resize(count);
			if (count>0) {
				for (int i = 0; i <count; i++) {
					OPS_GetIntInput(&dataNum, &data);
					EntityFaceID(i) = data;
				}
			}
#ifdef _DEBUG
			opserr << "Detect " << count << " Faces: " << EntityFaceID << endln;
#endif
			// for geting uncertain number of doubel values
		}
	}
	else if (strcmp(option, "-HTEleSet") == 0 || strcmp(option, "-EleSet") == 0 || strcmp(option, "eleset") == 0) {

		if (OPS_GetIntInput(&dataNum, &HTEleSetID) <0) {
			opserr << "WARNING:: invalid HTEleSet tag for defining HeatFluxBC: " << "\n";
			return -1;
		}
		
		theHTEleSet = theHTModule->getHTEleSet(HTEleSetID);

		if (theHTEleSet == 0) {
			opserr << "WARNING:: invalid HTEleSet tag for defining HeatFluxBC: " << "\n";
			return -1;
		}
		//if face tag is detected
		option = OPS_GetString();
		if (strcmp(option, "-face") == 0 || strcmp(option, "-Face") == 0 || strcmp(option, "-faceTag") == 0) {

			if (OPS_GetIntInput(&dataNum, &FaceID) <0) {
				opserr << "WARNING:: invalid face tag for defining FaceID: " << option << "\n";
				return -1;
			}

		}

		if (FaceID == 0) {
			opserr << "WARNING::no faceID specified for defining heat flux BC" << endln;
			opserr << " valid tags: -face, -Face, -faceTag " << endln;
			return -1;
		}

	}
	else {
		opserr << "WARNING TclHTCommand HeatFluxBC wants  -HTEntity, HTEntity, -Entity or -HTElesET, -EleSet, eleset " << endln;
		return -1;
	}
	//Simple_Boundary::GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType)

	//end of obtaining face or entity information

	//type is identified
	if (OPS_GetNumRemainingInputArgs() > 0) {
		option = OPS_GetString();
		if (strcmp(option, "-type") == 0 || strcmp(option, "-Type") == 0 || strcmp(option, "type") == 0) {
			if (OPS_GetNumRemainingInputArgs() <1) {
				opserr << "WARNING unknown heat flux type" << option;
				opserr << " valid types: -Convec, -Radiation,-ConvecAndRad,and -Prescribed \n";
				return -1;
			}
			option = OPS_GetString();
			if (strcmp(option, "-convec") == 0 || strcmp(option, "-Convec") == 0 || strcmp(option, "Convec") == 0) {
				HeatFluxTypeTag = 1;
			}
			else if (strcmp(option, "-Radiation") == 0 || strcmp(option, "-Rad") == 0 || strcmp(option, "Rad") == 0) {
				HeatFluxTypeTag = 2;
			}
			else if (strcmp(option, "-Prescribed") == 0 || strcmp(option, "-Pres") == 0 || strcmp(option, "Pres") == 0) {
				HeatFluxTypeTag = 3;
			}
			else if (strcmp(option, "-ConvecAndRad") == 0 || strcmp(option, "-convecrad") == 0 || strcmp(option, "ConvecRad") == 0) {
				HeatFluxTypeTag = 4;
			}
			else {
				opserr << "WARNING unknown heat flux type" << option;
				opserr << " valid types: -Convec, -Radiation,-ConvecAndRad,and -Prescribed \n";
				return -1;
			}

		}

	}
	

	if (HeatFluxTypeTag == 0) {
		opserr << "WARNING::no HeatFlux type specified for defining heat flux BC" << endln;
		opserr << " valid tags: -type, -Type, type " << endln;
	}

	//if Tag HTConstants is detected
	if (OPS_GetNumRemainingInputArgs() < 0) {
		opserr << "WARNING unknown HTConstants for HeatfluxBC" << option;
		return -1;
	}
	option = OPS_GetString();
	if (strcmp(option, "-HTConstants") == 0 || strcmp(option, "-htConstants") == 0 || strcmp(option, "-Constants") == 0) {

		if (OPS_GetIntInput(&dataNum, &HTconstantsID) <0) {
			opserr << "WARNING:: invalid constants tag for defining HTEleSet: " << option << "\n";
			return -1;
		}
	}

	if (HTconstantsID == 0) {
		opserr << "WARNING::no HTConstants found for defining heat flux BC" << endln;
		opserr << " valid tags: -HTConstants, -htConstants, -Constants " << endln;
	}

	//Create boudary pattern pointed to theHTPattern;
	BoundaryPattern* thePattern = theHTPattern;
	if (theHTPattern == 0) {
		opserr << "WARNING:: HeatTransferModule failed to find Boundary pattern" << endln;
	}
	//get the pattern tag
	int PatternTag = thePattern->getTag();


	//obtain heat flux Constants;
	HTConstants* theHTConstants = theHTModule->getHTConstants(HTconstantsID);
	if (theHTConstants == 0) {
		opserr << "WARNING::HeatTransferModule failed to find required HTConstants with tag: " << HTconstantsID << endln;
		return -1;
	}
	Vector HeatFluxConstants = theHTConstants->getConstants();

	if (HeatFluxConstants == 0) {
		opserr << "WARNING TclHTModule:addHeatFluxBC failed to get HTConstants " << HTconstantsID << endln;
	}

	//------------------Applying heat flux BC--------------------------//
	//-----------------------------------------------------------------//
	if (HTEntityTag != 0) {

		int EntityMeshTag = theHTEntity->getMeshTag();
		Simple_Mesh* theHTMesh = theHTModule->getHTMesh(EntityMeshTag);

		int NumFaces = EntityFaceID.Size();
		//loop over all the faces;

		for (int i = 0; i < NumFaces; i++) {
			//if face tag is detected
			ID ElesRange(0);
			int EntFaceID = EntityFaceID(i);
			int FaceID;
			theHTMesh->SelectingElesbyFace(ElesRange, EntFaceID, FaceID);

			//opserr << "return face " << FaceID << endln;
			//detect the num of existing HeatfluxBCs;
			int ExistingHeatFluxBCs = thePattern->getNumHeatFluxBCs();

			//check the number of heat flux BCs being applied
			int NumHeatFluxBCs = ElesRange.Size();
			if (ElesRange == 0 || NumHeatFluxBCs == 0)
			{
				opserr << "WARNING:: HeatTransferModule: No HeatFluxBc will be defined for Boundary pattern" << endln;
			}


			if (HeatFluxTypeTag == 1)
			{
				for (int i = 0; i < NumHeatFluxBCs; i++) {
					Convection* Convec_BC = new Convection(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, HeatFluxConstants(0), HeatFluxConstants(1));
					theHTDomain->addHeatFluxBC(Convec_BC, PatternTag);
				}
			}
			else if (HeatFluxTypeTag == 2)
			{

				for (int i = 0; i < NumHeatFluxBCs; i++) {

					Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, HeatFluxConstants(2), 5.67e-8, HeatFluxConstants(3), HeatFluxConstants(4));
					theHTDomain->addHeatFluxBC(Radiat_BC, PatternTag);

				}
			}
			else if (HeatFluxTypeTag == 3)
			{
				int NumNPF = (theHTDomain->getElement(ElesRange(0)))->getNumNodesperFace();
				for (int i = 0; i < NumHeatFluxBCs; i++) {
					PrescribedSurfFlux* PreSurfFlux_BC = new PrescribedSurfFlux(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, NumNPF, FireType);
					theHTDomain->addHeatFluxBC(PreSurfFlux_BC, PatternTag);
				}
			}
			else if (HeatFluxTypeTag == 4)
			{
				for (int i = 0; i < NumHeatFluxBCs; i++) {
					Convection* Convec_BC = new Convection(ExistingHeatFluxBCs + 2 * i, ElesRange(i), FaceID, HeatFluxConstants(0), HeatFluxConstants(1));
					theHTDomain->addHeatFluxBC(Convec_BC, PatternTag);

					Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs + 2 * i + 1, ElesRange(i), FaceID, HeatFluxConstants(2), 5.67e-8, HeatFluxConstants(3), HeatFluxConstants(4));
					theHTDomain->addHeatFluxBC(Radiat_BC, PatternTag);
				}
			}
#ifdef _DEBUG
			opserr << "HeatTransferModule:: " << NumHeatFluxBCs << " eles are applied with HeatFlux(Type: "
				<< HeatFluxTypeTag << " ) at the face " << EntFaceID << " ,Entity :" << HTEntityTag << endln;
#endif
			//end of heat flux type
		}
		//end of the loop over all faces

	}
	//end of if elements defined by entity face
	else if (HTEleSetID!= 0) {

		ElesRange = theHTEleSet->getEleID();
		BoundaryPattern* thePattern = theHTPattern;
		int PatternTag = thePattern->getTag();
		int NumHeatFluxBCs = ElesRange.Size();
		if (ElesRange == 0 || NumHeatFluxBCs == 0)
		{
			opserr << "WARNING:: HeatTransferModule: No HeatFluxBc will be defined for Boundary pattern" << endln;
		}

		int ExistingHeatFluxBCs = thePattern->getNumHeatFluxBCs();

		Vector HeatFluxConstants = (theHTModule->getHTConstants(HTconstantsID))->getConstants();

		if (HeatFluxConstants == 0) {
			opserr << "WARNING TclHTModule:addHeatFluxBC failed to get HTConstants " << HTconstantsID << endln;
		}

		if (HeatFluxTypeTag == 1)
		{
			for (int i = 0; i < NumHeatFluxBCs; i++) {
				Convection* Convec_BC = new Convection(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, HeatFluxConstants(0), HeatFluxConstants(1));
				theHTDomain->addHeatFluxBC(Convec_BC, PatternTag);
			}
		}
		else if (HeatFluxTypeTag == 2)
		{

			for (int i = 0; i < NumHeatFluxBCs; i++) {

				Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, HeatFluxConstants(2), HeatFluxConstants(3), HeatFluxConstants(4), HeatFluxConstants(5));
				theHTDomain->addHeatFluxBC(Radiat_BC, PatternTag);

			}
		}
		else if (HeatFluxTypeTag == 3)
		{
			int NumNPF = (theHTDomain->getElement(ElesRange(0)))->getNumNodesperFace();
			for (int i = 0; i < NumHeatFluxBCs; i++) {
				PrescribedSurfFlux* PreSurfFlux_BC = new PrescribedSurfFlux(ExistingHeatFluxBCs + i, ElesRange(i), FaceID, NumNPF, FireType);
				theHTDomain->addHeatFluxBC(PreSurfFlux_BC, PatternTag);
			}
		}
		else if (HeatFluxTypeTag == 4)
		{
			for (int i = 0; i < NumHeatFluxBCs; i++) {
				Convection* Convec_BC = new Convection(ExistingHeatFluxBCs + 2 * i, ElesRange(i), FaceID, HeatFluxConstants(0), HeatFluxConstants(1));
				theHTDomain->addHeatFluxBC(Convec_BC, PatternTag);

				Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs + 2 * i + 1, ElesRange(i), FaceID, HeatFluxConstants(2), HeatFluxConstants(3), HeatFluxConstants(4), HeatFluxConstants(5));
				theHTDomain->addHeatFluxBC(Radiat_BC, PatternTag);
			}
		}
#ifdef _DEBUG
		opserr << "HeatTransferModule:: " << NumHeatFluxBCs << " eles are applied with HeatFlux(Type: "
			<< HeatFluxTypeTag << " ) for the HTeleSet :" << HTEleSetID << endln;
#endif
	}
	return 0;
	//Simple_Boundary::GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType)

}




//Add HeatFluxBC
int
OPS_addMPTemperatureBC()
{

	//HeatFluxBC -eleset 1 -facetag 3 -type -ConvecAndRad -HTConstants 1;

	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
		return -1;
	}
	int masterIDtag = 0;
	int slaveIDtag = 0;
	ID masterID;
	ID slaveID;
	int dataNum = 1;

	int count = 1;
	const char* option = OPS_GetString();
	//if HTEntity tag is detected
	if (strcmp(option, "-HTNodeSet") == 0 || strcmp(option, "-NodeSet") == 0 || strcmp(option, "nodeset") == 0)
	{

		if (OPS_GetIntInput(&dataNum, &masterIDtag) < 0) {
			opserr << "WARNING:: invalid masterID tag for coupling temperature: " << option << "\n";
			return -1;
		}


		if (OPS_GetIntInput(&dataNum, &slaveIDtag) < 0) {
			opserr << "WARNING:: invalid masterID tag for coupling temperature: " << option << "\n";
			return -1;
		}


		HTNodeSet* theMasterNodeSet = theHTModule->getHTNodeSet(masterIDtag);
		if (theMasterNodeSet == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested master NodeSet for coupling temperature: " << option << "\n";
			return -1;
		}

		HTNodeSet* theSlaveNodeSet = theHTModule->getHTNodeSet(slaveIDtag);
		if (theSlaveNodeSet == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested slave NodeSet for coupling temperature: " << option << "\n";
			return -1;
		}

		masterID = theMasterNodeSet->getNodeID();
		slaveID = theSlaveNodeSet->getNodeID();

	}
	// END OF OBTAINING NODES FORE COUPLING

	  //------APPLYING MP_TEMPERATUREBC

	if (masterID.Size() == 0 || slaveID.Size() == 0) {
		opserr << "Waring:No node will be coupled for tmperature" << endln;
	}
	else if (masterID.Size() != slaveID.Size())
	{
		opserr << "WARNING::The size of master node set " << masterIDtag
			<< " doesn't match the slave node set " << slaveIDtag << endln;
		return -1;
	}
	else
	{
		for (int i = 0; i < masterID.Size(); i++) {
			//Heat Transfer has only one DOF
			MP_TemperatureBC* MP_TempBC = new MP_TemperatureBC(masterID(i), slaveID(i));
			if (MP_TempBC == 0)
				opserr << "WARNING: ran out of memory for coupling temperature" << endln;

			if (theHTDomain->addMP_TemperatureBC(MP_TempBC) == false)
			{
				opserr << "Warning: could not add MP_TemperatureBC to domain" << endln;
				delete MP_TempBC;
			}
		}
	}

	return 0;
	//end of adding mp temperature BC;
}
//



//Add HT_transient analysis
int 
OPS_HTAnalysis()
{
	// make sure at least one other argument to contain type of system
	if (OPS_GetNumRemainingInputArgs()< 1) {
		opserr << "WARNING need to specify an analysis type (Static, Transient)\n";
		return -1;
	}
	if (theHTAnalysis != 0) {
		delete theHTAnalysis;
		theHTAnalysis = 0;
	}
	const char* option = OPS_GetString();
	if (strcmp(option, "HeatTransfer") == 0) {
		// make sure all the components have been built,
		// otherwise print a warning and use some defaults
		if (theAnalysisModel == 0)
			theAnalysisModel = new HT_AnalysisModel();

		if (theTest == 0) {
			opserr << "WARNING analysis Transient - no convergence test yet specified, \n";
			opserr << " CTestNormTempIncr default will be used\n";
			theTest = new CTestNormTempIncr(1e-4, 1000, 0);

		}

		if (theAlgorithm == 0) {
			opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
			opserr << " NewtonMethod default will be used\n";

			theAlgorithm = new NewtonMethod(*theTest);
		}
		if (theHandler == 0) {
			opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
			opserr << " yet specified, PenaltyBC_Handler default will be used\n";
			theHandler = new PenaltyBC_Handler(1.0e10, 1.0e10);
		}
		if (theNumberer == 0) {
			opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
			opserr << " RCM default will be used\n";
			RCM* theRCM = new RCM(false);
			theNumberer = new HT_DOF_Numberer(*theRCM);
		}
		if (theTransientIntegrator == 0) {
			opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
			opserr << " Backward difference default will be used\n";
			theTransientIntegrator = new BackwardDifference();
		}
		if (theSOE == 0) {
			opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
			opserr << " ProfileSPDLinSOE default will be used\n";
			BandGenLinLapackSolver* theSolver;
			theSolver = new BandGenLinLapackSolver();
			theSOE = new BandGenLinSOE(*theSolver);
		}

		theHTAnalysis = new HT_TransientAnalysis(*theHTDomain,
			*theHandler,
			*theNumberer,
			*theAnalysisModel,
			*theAlgorithm,
			*theSOE,
			*theTransientIntegrator,
			theTest);

	}
	return 0;
}


//analyze HT_transient model
int OPS_HTAnalyze()
{
	int result = 0;
	int numData = 1;
	double monitortime = 0;
	double Ltime = lasttime;

	if (theHTAnalysis != 0) {

		if (OPS_GetNumRemainingInputArgs() < 2) {
			opserr << "WARNING heat transfer analysis: analysis numIncr? deltaT?\n";
			return -1;
		}

		int numIncr;
		if (OPS_GetIntInput(&numData, &numIncr) < 0) {
			opserr << "WARNING heat transfer analysis: recieving incorret number of increment\n";
			return -1;
		}
			

		double dT;
		if (OPS_GetDoubleInput(&numData, &dT)<0) {
			opserr << "WARNING heat transfer analysis: recieving incorret time step\n";
			return -1;
		}

		if (OPS_GetNumRemainingInputArgs() > 0) {

			const char* option = OPS_GetString();
			if (strcmp(option, "-monitor") == 0) {
				if (OPS_GetDoubleInput(&numData, &monitortime) < 0) {
					opserr << "WARNING heat transfer analysis: recieving incorret time step\n";
					return -1;
				}
			}

		}


		result = theHTAnalysis->analyze(numIncr, dT,Ltime, monitortime);
		lasttime = Ltime;
		if (result < 0)
			opserr << "OpenSees > analyze failed, returned: " << result << " error flag\n";


	}
	//Output for the final result of heat transfer analysis;

#ifdef _DEBUG
	opserr <<"The heat transfer analysis:"<< result << endln;
#endif

	return 0;
}



int 
OPS_HTRecorder()
{
	int numData = 1;
	Vector* theRecVec = 0;
	if (theHTModule == 0) {
		opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
		return -1;
	}

	if (theHTDomain == 0) {
		opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
		return -1;
	}
	OPS_Stream* theOutputStream = 0;
	const char* fileName = 0;
	ID* theNodes = 0;
	HTRecorder* theHTRecorder = 0;
	double xloc = 0;

	//if file tag is detected
	const char* option = OPS_GetString();
	if (strcmp(option, "-file") == 0 || strcmp(option, "file") == 0 || strcmp(option, "-File") == 0)
	{
		fileName = OPS_GetString();
		//const char* pwd = getInterpPWD(interp);
		//simulationInfo.addOutputFile(fileName, pwd);
		
	}
	else {
		opserr << "WARNING::HTRecorder has not specified a file name" << endln;
		//OPS_GetNumRemainingInputArgs(-1);
		return - 1;
	}
	
	theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 0);  //The last argument is to determine CSV for using space 

	//if xloc tag is detected
	option = OPS_GetString();
	if (strcmp(option, "-xloc") == 0 || strcmp(option, "xLoc") == 0 || strcmp(option, "-xLoc") == 0)
	{
	
		if (OPS_GetDoubleInput(&numData, &xloc) ==0) {
			opserr << "WARNING invalid xloc" << endln;
			opserr << " for HeatTransfer recorder: " << endln;
			return -1;
		}
		
	
		//ending of xloc;
		option = OPS_GetString();
		if (strcmp(option, "-yloc") == 0 || strcmp(option, "yLoc") == 0 || strcmp(option, "-yLoc") == 0)
		{
			int numRArgs = 0;
			Vector yloc;

			if (OPS_GetNumRemainingInputArgs() > 0) {
				numRArgs = OPS_GetNumRemainingInputArgs();
				yloc.resize(numRArgs);
			}
			else {
				opserr << "WARNING:: HTRecorder can not find yloc\n";
				return -1;
			}
			//-----for geting uncertain number of double data.
			double data;

			for (int i = 0; i < numRArgs; i++) {
				OPS_GetDoubleInput(&numData, &data);
				yloc(i) = data;
			}


			int NumRows = yloc.Size();
			theRecVec = new Vector(NumRows + 1);
			theRecVec->Zero();

			for (int i = 0; i <= NumRows; i++) {
				if (i == 0) {
					(*theRecVec)(i) = xloc;
				}
				else {
					(*theRecVec)(i) = yloc(i - 1);
				}
			}
		}

#ifdef _DEBUG
		opserr << "TclHeatTransferModule::HTRecorder, theRecMatrix " << *theRecVec << endln;
#endif
		HTReorderTag++;
		theHTRecorder = new HTRecorderToStru(HTReorderTag, *theRecVec, *theHTDomain, *theOutputStream);
	}
	//end of if xloc is found;
	else if (strcmp(option, "-nodeSet") == 0 || strcmp(option, "NodeSet") == 0 || strcmp(option, "-NodeSet") == 0)
	{

		int NodeSetID = 0;
		if (OPS_GetIntInput(&numData, &NodeSetID) <0) {
			opserr << "WARNING:: invalid nodeSet tag for defining HTNodeRecorder : " << "\n";
			return -1;
		}

		HTNodeSet* theRecNodeSet = theHTModule->getHTNodeSet(NodeSetID);
		if (theRecNodeSet == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested NodeSet for HTNodeRecorder: " << "\n";
			return -1;
		}
		ID RecNodeID = 0;
		RecNodeID = theRecNodeSet->getNodeID();

#ifdef _DEBUG
		opserr << "TclHeatTransferModule::HTRecorder, theRecNodeID " << RecNodeID << endln;
#endif
		theHTRecorder = new HTNodeRecorder(HTReorderTag++, &RecNodeID, *theHTDomain, *theOutputStream);

	}
	//end of nodeset
	else if (strcmp(option, "-HTEntity") == 0 || strcmp(option, "htEntity") == 0 || strcmp(option, "-htentity") == 0)
	{
		int EntityID = 0; int DimTag = 0;
		if (OPS_GetIntInput(&numData, &EntityID) <0) {
			opserr << "WARNING:: invalid Entity tag for defining HTNodeRecorder : " << "\n";
			return -1;
		}
		//The following code is to specify the tag for num of dimensions
		if (OPS_GetNumRemainingInputArgs()>0) {
			option = OPS_GetString();
			if (strcmp(option, "-dim") == 0 || strcmp(option, "dim") == 0 || strcmp(option, "Dim") == 0) {
				if (OPS_GetIntInput(&numData, &DimTag) <0) {
					opserr << "WARNING:: invalid Entity tag for defining HTNodeRecorder : " << "\n";
					return -1;
				}
			}

			
		}

		Simple_Entity* theRecEntity = theHTModule->getHTEntity(EntityID);
		int theMshTag = theRecEntity->getMeshTag();
		Simple_Mesh* theHTMesh = theHTModule->getHTMesh(theMshTag);

		if (theHTMesh == 0) {
			opserr << "WARNING:: TclHTModule failed to get the requested Entity for HTNodeRecorder: " << "\n";
			return -1;
		}

		ID RecNodesID = 0;
		theHTMesh->GetNodesForRecorder(RecNodesID, DimTag);

#ifdef _DEBUG
		opserr << "TclHeatTransferModule::HTRecorder, theRecNodeID " << RecNodesID << endln;
#endif
		theHTRecorder = new HTNodeRecorder(HTReorderTag++, &RecNodesID, *theHTDomain, *theOutputStream);

	}


	//HTRecorderToStru::HTRecorderToStru(int tag, const Matrix& theCrds, HeatTransferDomain &theDom, OPS_Stream &theOutputHandler, double tolerance)
// for geting uncertain number of doubel values 

	if (theHTRecorder != 0) {
		theHTDomain->addRecorder(*theHTRecorder);
	}
	else {
		opserr << "WARNING::HeatTransfer Recoder is not defined" << endln;
	}
	// HTNodeRecorder(int tag, const ID* theNodes, HeatTransferDomain& theDomain,OPS_Stream &theOutputHandle); 
	
	return 0;
}


int OPS_getHTtime()
{
	int numdata = 1;
	if (theHTDomain == 0) return -1;
	double factor = theHTDomain->getCurrentTime();

	if (OPS_SetDoubleOutput(&numdata, &factor) < 0) {
		opserr << "WARNING failed to set load factor\n";
		return -1;
	}

	return 0;
}

int OPS_HTOutput()
{
	if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING HTOutput: insufficient argument for retruning a value\n";
		return -1;
	}
	int dataNum = 1;
	const char* option = 0;
	option = OPS_GetString();
	
	if (strcmp(option, "firemodel") == 0 || strcmp(option, "fireModel") == 0 || strcmp(option, "-fire") == 0)
	{
		
		int FireModelTag = 0;
		int FireParTag = 0;

		if (OPS_GetIntInput(&dataNum, &FireModelTag) < 0) {
			opserr << "WARNING:: invalid FireModel tag for HTOutput: " << "\n";
			return -1;
		}
		if (OPS_GetNumRemainingInputArgs() < 1)
		{
			opserr << "WARNING:: insufficient arguments for fire model HTOutput: " << "\n";
			return -1;
		}
		option = OPS_GetString();
		if (strcmp(option, "firepars") == 0 || strcmp(option, "firePars") == 0 || strcmp(option, "-firePars") == 0)
		{
			
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetIntInput(&dataNum, &FireParTag) < 0) {
					opserr << "WARNING:: invalid FireModel tag for HTOutput: " << "\n";
					return -1;
				}
			}

		}
		FireModel* thefire = 0;
		thefire = theHTModule->getFireModel(FireModelTag);
		if (thefire == 0) {
			opserr << "WARNING:: HToutput can not find the fire model " << endln;
			return -1;
		}
		double firePar = thefire->getFirePars(FireParTag);
		//opserr << "fireparTag " << FireParTag << "firepar " << firePar << endln;
		if (OPS_SetDoubleOutput(&dataNum, &firePar) < 0) {
			opserr << "WARNING failed to return fire pars for fire model "<< FireModelTag<<"\n";
			return -1;
		}
	}
	//end of firemodel
	else if (strcmp(option, "-node") == 0 || strcmp(option, "Node") == 0 || strcmp(option, "node") == 0)
	{
		if (OPS_GetNumRemainingInputArgs() < 2)
		{
			opserr << "WARNING:: insufficient arguments for nodal HTOutput: " << "\n";
			return -1;
		}
		int nodetag = 0;
		if (OPS_GetIntInput(&dataNum, &nodetag) < 0) {
			opserr << "WARNING:: invalid node tag for HTOutput: " << "\n";
			return -1;
		}
		HeatTransferNode* theNode = theHTDomain->getNode(nodetag);

		if (OPS_GetNumRemainingInputArgs() > 0)
		{
			option = OPS_GetString();
		}
		else
			option = "temp";

		if (strcmp(option, "-temp") == 0 || strcmp(option, "temp") == 0 || strcmp(option, "Temp") == 0) {
			double NodeTemp = theNode->getTemperature()(0)-273.15;
			if (OPS_SetDoubleOutput(&dataNum, &NodeTemp) < 0) {
				opserr << "WARNING failed to return temperature for node " << nodetag << "\n";
				return -1;
			}
		}
	}
	else if (strcmp(option, "-nodeset") == 0 || strcmp(option, "NodeSet") == 0 || strcmp(option, "nodeset") == 0)
	{

		if (OPS_GetNumRemainingInputArgs() < 2)
		{
			opserr << "WARNING:: insufficient arguments for nodal HTOutput: " << "\n";
			return -1;
		}
		int nodesetTag = 0;
		if (OPS_GetIntInput(&dataNum, &nodesetTag) < 0) {
			opserr << "WARNING:: invalid nodeset tag for HTOutput: " << "\n";
			return -1;
		}
		HTNodeSet* theNodeset = theHTModule->getHTNodeSet(nodesetTag);

		if (OPS_GetNumRemainingInputArgs() > 0)
		{
			option = OPS_GetString();
		}
		else
			option = "temp";

		if (strcmp(option, "-temp") == 0 || strcmp(option, "temp") == 0 || strcmp(option, "Temp") == 0) {
			ID nodeID = theNodeset->getNodeID();
			int numNodes = nodeID.Size();
			double* data = new double[numNodes];

			for (int i = 0; i < numNodes; i++) {
				int nodeTag = nodeID(i);
				HeatTransferNode* theNode = theHTDomain->getNode(nodeTag);
				data[i]= theNode->getTemperature()(0) - 273.15;
			}

			if (OPS_SetDoubleOutput(&numNodes, data) < 0) {
				opserr << "WARNING failed to return temperature for node set" << nodesetTag << "\n";
				delete[] data;
				return -1;
			}
		}
	}
	


	return 0;
}

int OPS_SetFirePars() {
	if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING HTOutput: insufficient argument for retruning a value\n";
		return -1;
	}
	int dataNum = 1;
	const char* option = 0;
	option = OPS_GetString();

	if (strcmp(option, "firemodel") == 0 || strcmp(option, "-fireModel") == 0 || strcmp(option, "-fire") == 0)
	{

		int FireModelTag = 0;
		int FireParTag = 0;
		double xloc = 0;
		double yloc = 0;
		double zloc = 0;
		double q = 0;
		double d = 0;
		double Ts = 0;
		double addq = 0;
		if (OPS_GetIntInput(&dataNum, &FireModelTag) < 0) {
			opserr << "WARNING:: invalid FireModel tag for HTOutput: " << "\n";
			return -1;
		}
		if (OPS_GetNumRemainingInputArgs() < 1)
		{
			opserr << "WARNING:: insufficient arguments for fire model HTOutput: " << "\n";
			return -1;
		}
		option = OPS_GetString();
		if (strcmp(option, "fireloc") == 0 || strcmp(option, "fireLoc") == 0 || strcmp(option, "-fireLoc") == 0)
		{

			if (OPS_GetNumRemainingInputArgs() < 3) {
				opserr<< "WARNING:: insufficient arguments for set fireloc: " << "\n";
			}
			else{
				if (OPS_GetDoubleInput(&dataNum, &xloc) < 0) {
					opserr << "WARNING:: invalid xloc for set firePars " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &yloc) < 0) {
					opserr << "WARNING:: invalid yloc for set firePars " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &zloc) < 0) {
					opserr << "WARNING:: invalid zloc for set firePars " << "\n";
					return -1;
				}
			}
			Vector firelocs(3);
			firelocs(0) = xloc; firelocs(1) = yloc; firelocs(2) = zloc;
			FireModel* thefire = theHTModule->getFireModel(FireModelTag);
			double thecurrentTime = theHTDomain->getCurrentTime();
			thefire->setFirePars(thecurrentTime, firelocs);

		}
		else if (strcmp(option, "firePars") == 0 || strcmp(option, "firepars") == 0 || strcmp(option, "-firepars") == 0 || strcmp(option, "-FirePars") == 0)
		{
			int numPars = 0;
			if (OPS_GetNumRemainingInputArgs() < 5) {
				opserr << "WARNING:: insufficient arguments for set fireloc: x,y,z,q,d" << "\n";
			}
			else {
				if (OPS_GetDoubleInput(&dataNum, &xloc) < 0) {
					opserr << "WARNING:: invalid xloc for set firePars " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &yloc) < 0) {
					opserr << "WARNING:: invalid yloc for set firePars " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &zloc) < 0) {
					opserr << "WARNING:: invalid zloc for set firePars " << "\n";
					return -1;
				}
				if (OPS_GetNumRemainingInputArgs() > 0) {
					if (OPS_GetDoubleInput(&dataNum, &q) < 0) {
						opserr << "WARNING:: invalid q for set firePars" << "\n";
						return -1;
					}
					numPars = 4;
				}
				else {
					numPars = 3;
				}
				
				if (OPS_GetNumRemainingInputArgs() > 0) {
					if (OPS_GetDoubleInput(&dataNum, &d) < 0) {
						opserr << "WARNING:: invalid d for set firePars " << "\n";
						return -1;
					}
					numPars = 5;
				}
				
				if (OPS_GetNumRemainingInputArgs() > 0) {
					if (OPS_GetDoubleInput(&dataNum, &Ts) < 0) {
						opserr << "WARNING:: invalid T_smoke for set firePars " << "\n";
						return -1;
					}
					numPars = 6;
				}
				
				if (OPS_GetNumRemainingInputArgs() > 0) {
					if (OPS_GetDoubleInput(&dataNum, &addq) < 0) {
						opserr << "WARNING:: invalid maximum incident q for set firePars " << "\n";
						return -1;
					}
					numPars = 7;
				}

			}
			Vector firepars(numPars);
			firepars(0) = xloc; firepars(1) = yloc; firepars(2) = zloc; 
			if (numPars == 4) {
				firepars(3) = q; 
			}
			else if (numPars == 5) {
				firepars(3) = q; firepars(4) = d;
			}
			else if (numPars == 6) {
				firepars(3) = q; firepars(4) = d; firepars(5) = Ts+273.15;
			}
			else if (numPars == 7) {
				firepars(3) = q; firepars(4) = d; firepars(5) = Ts + 273.15; firepars(6) = addq;
			}
			
			FireModel* thefire = 0;
			thefire = theHTModule->getFireModel(FireModelTag);
			if (thefire == 0) {
				opserr << "HeatTransferModule::SetFirePars failed to obtain the fire model" << endln;
				return -1;
			}
			double thecurrentTime = theHTDomain->getCurrentTime();
		
			thefire->setFirePars(thecurrentTime, firepars);
		}
		
		
	}


}


int OPS_GetFireOut() {
	if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING HTOutput: insufficient argument for retruning a value\n";
		return -1;
	}
	int dataNum = 1;
	const char* option = 0;
	option = OPS_GetString();

	if (strcmp(option, "firemodel") == 0 || strcmp(option, "fireModel") == 0 || strcmp(option, "-fire") == 0)
	{

		int FireModelTag = 0;
		int FireParTag = 0;
		double xloc = 0;
		double yloc = 0;
		double zloc = 0;
		double q = 0;
		double d = 0;
		double Ts = 0;
		double ptime = 0.0;
		if (OPS_GetIntInput(&dataNum, &FireModelTag) < 0) {
			opserr << "WARNING:: invalid FireModel tag for HTOutput: " << "\n";
			return -1;
		}
		if (OPS_GetNumRemainingInputArgs() < 1)
		{
			opserr << "WARNING:: insufficient arguments for fire model HTOutput: " << "\n";
			return -1;
		}
		option = OPS_GetString();
		if (strcmp(option, "loc") == 0 || strcmp(option, "Loc") == 0 || strcmp(option, "-Loc") == 0)
		{

			if (OPS_GetNumRemainingInputArgs() < 3) {
				opserr << "WARNING:: insufficient arguments for set loc: " << "\n";
			}
			else {
				if (OPS_GetDoubleInput(&dataNum, &xloc) < 0) {
					opserr << "WARNING:: invalid xloc for set xloc " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &yloc) < 0) {
					opserr << "WARNING:: invalid yloc for set yloc " << "\n";
					return -1;
				}
				if (OPS_GetDoubleInput(&dataNum, &zloc) < 0) {
					opserr << "WARNING:: invalid zloc for set zloc " << "\n";
					return -1;
				}
			}
			Vector locs(3);
			locs(0) = xloc; locs(1) = yloc; locs(2) = zloc;
			FireModel* thefire = theHTModule->getFireModel(FireModelTag);
			double thecurrentTime = theHTDomain->getCurrentTime();
#ifdef _DEBUG
			opserr << "HTModule Checking the fire model " << FireModelTag << "location:" << locs << endln;
#endif

			double incq = thefire->getFireOut(60.0, locs);

			if (OPS_SetDoubleOutput(&dataNum, &incq) < 0) {
				opserr << "WARNING failed to return fire pars for fire model " << FireModelTag << "\n";
				return -1;
			}

		}

	}


}

//Add HTWipe
int OPS_HTReset()
{
	theHTDomain->clearAll();

	theHTDomain = 0;
	theHTModule = 0;
	theHTPattern = 0;
	return 0;
}

//Add HTWipe
int OPS_WipeHT()
{
	if (theHTModule != 0) {
		delete theHTModule;
		theHTModule = 0;
	}

	
	return 0;
}



//////////////////////////////////////////////
/////// Python HTModule functions  ////////////
/////////////////////////////////////////////
static PyObject* Py_ops_HeatTransfer(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HeatTransfer() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHTMaterial(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTMaterial() < 0) {
		opserr << (void*)0;
		//return NULL;
	}

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHTEntity(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTEntity() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHTMesh(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTMesh() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTMeshAll(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTMeshAll() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_SetInitialT(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_SetInitialT() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTRefineMesh(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTRefineMesh() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHTConstants(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTConstants() < 0) return NULL;

	return theWrapper->getResults();
}


static PyObject* Py_ops_HTNodeSet(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTNodeSet() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTEleSet(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTEleSet() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHTPattern(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHTPattern() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addFireModel(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addFireModel() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_addHeatFluxBC(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addHeatFluxBC() < 0) return NULL;

	return theWrapper->getResults();
}



static PyObject* Py_ops_addMPTemperatureBC(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_addMPTemperatureBC() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTAnalysis(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTAnalysis() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTAnalyze(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTAnalyze() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTRecorder(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTRecorder() < 0) return NULL;

	return theWrapper->getResults();
}


static PyObject* Py_ops_HTtime(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_getHTtime() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HToutput(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTOutput() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_SetFirePars(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_SetFirePars() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_GetFireOut(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_GetFireOut() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_HTreset(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_HTReset() < 0) return NULL;

	return theWrapper->getResults();
}

static PyObject* Py_ops_wipeHT(PyObject* self, PyObject* args)
{
	theWrapper->resetCommandLine(PyTuple_Size(args), 1, args);

	if (OPS_WipeHT() < 0) return NULL;

	return theWrapper->getResults();
}


//-------------------------------------------------------------------------------//
//------------------------Procedure of Heat transfer commmands-----------------------------//


int OPS_addHTCommands(PythonWrapper* thewrapper)
{
	theWrapper = thewrapper;
	theWrapper->addCommand("HeatTransfer", &Py_ops_HeatTransfer);
	theWrapper->addCommand("HTMaterial", &Py_ops_addHTMaterial);
	theWrapper->addCommand("HTEntity", &Py_ops_addHTEntity);
	theWrapper->addCommand("HTMesh", &Py_ops_addHTMesh);
	
	theWrapper->addCommand("HTMeshAll", &Py_ops_HTMeshAll);
	theWrapper->addCommand("SetInitialT", &Py_ops_SetInitialT);
	theWrapper->addCommand("HTRefineMesh", &Py_ops_HTRefineMesh);
	theWrapper->addCommand("HTConstants", &Py_ops_addHTConstants);
	
	theWrapper->addCommand("HTEleSet", &Py_ops_HTEleSet);
	theWrapper->addCommand("HTNodeSet", &Py_ops_HTNodeSet);
	theWrapper->addCommand("HTPattern", &Py_ops_addHTPattern);
	theWrapper->addCommand("FireModel", &Py_ops_addFireModel);

	theWrapper->addCommand("HToutput", &Py_ops_HToutput);
	theWrapper->addCommand("HTTime", &Py_ops_HTtime);
	theWrapper->addCommand("SetFirePars", &Py_ops_SetFirePars);
	theWrapper->addCommand("GetFireOut", &Py_ops_GetFireOut);

	theWrapper->addCommand("HeatFluxBC", &Py_ops_addHeatFluxBC);
	theWrapper->addCommand("HTCouple", &Py_ops_addMPTemperatureBC);

	theWrapper->addCommand("HTAnalysis", &Py_ops_HTAnalysis);
	theWrapper->addCommand("HTRecorder", &Py_ops_HTRecorder);
	theWrapper->addCommand("HTAnalyze", &Py_ops_HTAnalyze);
	
	theWrapper->addCommand("HTReset", &Py_ops_HTreset);
	theWrapper->addCommand("WipeHT", &Py_ops_wipeHT);

	return 0;
}


