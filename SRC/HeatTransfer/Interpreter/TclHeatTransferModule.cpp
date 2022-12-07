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
// Written by  Liming Jiang (Liming.Jiang@ed.ac.uk)

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
#include <ModifiedNewtonMethod.h>
#include <CTestNormTempIncr.h>
#include <CTestNormResidual.h>
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
#include < NaturalFire.h>
#include <AlpertCeilingJetModel.h>
#include <UserDefinedFire.h>
#include <BoundaryPattern.h>
#include <FireImposedPattern.h>

//include recorders
#include <HTNodeRecorder.h>
#include <HTRecorderToStru.h>
#include <HTElementRecorder.h>

#include <elementAPI.h>
#include <TclHeatTransferModule.h>

 #include <StandardStream.h>
 #include <DataFileStream.h>

 #include <SimulationInformation.h>
 #include  <PathTimeSeriesThermal.h>

 extern SimulationInformation simulationInfo;
 extern const char * getInterpPWD(Tcl_Interp *interp);  // commands.cpp

	
static HeatTransferDomain* theHTDomain = 0;
static TclHeatTransferModule* theTclHTModule =0;
static BoundaryPattern* theTclHTPattern=0;

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

//Declaration of Heat transfer Tcl commands
int TclHeatTransferCommand_addHTMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addHTEntity(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addHTMesh(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv); 
int TclHeatTransferCommand_HTMeshAll(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_SetInitialT(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addHTConstants(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_HTNodeSet(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_HTEleSet(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addHTPattern(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addFireModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_setFirePars(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv);
int TclHeatTransferCommand_addHeatFluxBC(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_addMPTemperatureBC(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_HTRefineMesh(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


int TclHeatTransferCommand_HTReset(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int TclHeatTransferCommand_HTRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_HTTest(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv);
int TclHeatTransferCommand_HTAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclHeatTransferCommand_HTAnalyze(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int TclHeatTransferCommand_PrintNodes(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv);
int TclHeatTransferCommand_getHTTime(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv);

TclHeatTransferModule::TclHeatTransferModule(int ndm, Tcl_Interp* interp)
{
  // Set the interpreter (the destructor needs it to delete commands)
  theInterp = interp;
  
 // Tcl_CreateCommand(interp, "HTNode",	TclHeatTransferModule_addHTNode,(ClientData)NULL, NULL);
  opserr<<"               ----HeatTransfer Module is developed by UoE Group----       "<<endln;


  theTclHTModule = this;

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

  Tcl_CreateCommand(interp, "HTMaterial", (Tcl_CmdProc* )TclHeatTransferCommand_addHTMaterial,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTEntity", (Tcl_CmdProc* )TclHeatTransferCommand_addHTEntity,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTMesh", (Tcl_CmdProc* )TclHeatTransferCommand_addHTMesh,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTMeshAll", (Tcl_CmdProc* )TclHeatTransferCommand_HTMeshAll,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SetInitialT", (Tcl_CmdProc* )TclHeatTransferCommand_SetInitialT,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTRefineMesh", (Tcl_CmdProc* )TclHeatTransferCommand_HTRefineMesh,(ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "HTConstants", (Tcl_CmdProc* )TclHeatTransferCommand_addHTConstants,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTEleSet", (Tcl_CmdProc* )TclHeatTransferCommand_HTEleSet,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTNodeSet", (Tcl_CmdProc* )TclHeatTransferCommand_HTNodeSet,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTPattern", (Tcl_CmdProc* )TclHeatTransferCommand_addHTPattern,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "FireModel", (Tcl_CmdProc* )TclHeatTransferCommand_addFireModel,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SetFirePars", (Tcl_CmdProc*)TclHeatTransferCommand_setFirePars,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HeatFluxBC", (Tcl_CmdProc* )TclHeatTransferCommand_addHeatFluxBC,(ClientData)NULL, NULL);
  //Tcl_CreateCommand(interp, "HTFixT", (Tcl_CmdProc*)TclHeatTransferCommand_addSPTemperatureBC, (ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTCoupleT", (Tcl_CmdProc* )TclHeatTransferCommand_addMPTemperatureBC,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTAnalysis", (Tcl_CmdProc* )TclHeatTransferCommand_HTAnalysis,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTTest", (Tcl_CmdProc*)TclHeatTransferCommand_HTTest, (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "HTRecorder", (Tcl_CmdProc* )TclHeatTransferCommand_HTRecorder,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTAnalyze", (Tcl_CmdProc* )TclHeatTransferCommand_HTAnalyze,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "HTReset", (Tcl_CmdProc* )TclHeatTransferCommand_HTReset,(ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "HTPrintNodes", (Tcl_CmdProc*)TclHeatTransferCommand_PrintNodes, (ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "getHTTime", (Tcl_CmdProc*)TclHeatTransferCommand_getHTTime, (ClientData)NULL, NULL);
  //Tcl_CreateCommand(interp, "HTMaterial", (Tcl_ObjCmdProc*) TclHeatTransferCommand_addHTMaterial,(ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);



}

TclHeatTransferModule::~TclHeatTransferModule()
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
 theTclHTModule =0;
 theTclHTPattern=0;
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
  

	Tcl_DeleteCommand(theInterp, "HTMaterial");
	Tcl_DeleteCommand(theInterp, "HTEntity");
	Tcl_DeleteCommand(theInterp, "HTMesh");
  Tcl_DeleteCommand(theInterp, "HTMeshAll");
  Tcl_DeleteCommand(theInterp, "SetInitialT");
   Tcl_DeleteCommand(theInterp, "HTRefineMesh");

  Tcl_DeleteCommand(theInterp, "HTConstants");
  Tcl_DeleteCommand(theInterp, "HTEleSet");
   Tcl_DeleteCommand(theInterp, "HTNodeSet");
  Tcl_DeleteCommand(theInterp, "HTPattern");
  Tcl_DeleteCommand(theInterp, "FireModel");
  Tcl_DeleteCommand(theInterp, "SetFirePars");

  Tcl_DeleteCommand(theInterp, "HeatFuxBC");
  Tcl_DeleteCommand(theInterp, "HTCoupleT");

  Tcl_DeleteCommand(theInterp, "HTAnalysis");
  Tcl_DeleteCommand(theInterp, "HTTest");
  Tcl_DeleteCommand(theInterp, "HTAnalyze");

}


int 
TclHeatTransferModule::addHTMaterial(HeatTransferMaterial &theHTMaterial)
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
TclHeatTransferModule::getHTMaterial(int tag)
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
TclHeatTransferModule::addHTEntity(Simple_Entity &theHTEntity)
{
  bool result = theHTEntities->addComponent(&theHTEntity);
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addHTEntity() - failed to add Entity: " << theHTEntity;
    return -1;
  }
}

Simple_Entity *
TclHeatTransferModule::getHTEntity(int tag)
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
TclHeatTransferModule::addHTMesh(Simple_Mesh *theHTMesh)
{
  
  // first check if a mesh with a similar tag exists in model
  int tag = theHTMesh->getTag();
  TaggedObject *other = theHTMeshes->getComponentPtr(tag);
  if (other != 0) {
    opserr << "TclHeatTransferModule::addHTMesh - cannot add as HTMesh with tag" <<
    tag << "already exists in model\n";
    
    return false;
  }
  
  // now we add the load pattern to the container for load pattrens
  bool result = theHTMeshes->addComponent(theHTMesh);
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addHTMesh() - failed to add Mesh: " << theHTMesh;
    return -1;
  }
}

Simple_Mesh *
TclHeatTransferModule::getHTMesh(int tag)
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
TclHeatTransferModule::addHTConstants(HTConstants* theHTConstants)
{
  
  // now we add the load pattern to the container for load pattrens
  bool result = theHTCons->addComponent(theHTConstants);
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addHTConstants() - failed to add HTConstants: " << theHTConstants;
    return -1;
  }
}

// To get the pointer back
HTConstants*
TclHeatTransferModule::getHTConstants(int tag)
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
TclHeatTransferModule::addHTEleSet(HTEleSet* theHTEleSet)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theHTEleSets->addComponent(theHTEleSet);
  
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addHTEleSet() - failed to add HTEleSet: " << theHTEleSet;
    return -1;
  }
}

HTEleSet*
TclHeatTransferModule::getHTEleSet(int tag)
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
TclHeatTransferModule::addHTNodeSet(HTNodeSet* theHTNodeSet)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theHTNodeSets->addComponent(theHTNodeSet);
  
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addHTNodeSet() - failed to add HTNodeSet: " << theHTNodeSet;
    return -1;
  }
}

HTNodeSet*
TclHeatTransferModule::getHTNodeSet(int tag)
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
TclHeatTransferModule::addFireModel(FireModel* theFireModel)
{
  // now we add the load pattern to the container for load pattrens
  bool result = theFireModels->addComponent(theFireModel);
  
  if (result == true)
    return 0;
  else {
    opserr << "TclHeatTransferModule::addFireModel() - failed to add FireModel: " << theFireModel;
    return -1;
  }
}

FireModel*
TclHeatTransferModule::getFireModel(int tag)
{
  TaggedObject *mc = theFireModels->getComponentPtr(tag);
  if (mc == 0)
    return 0;
  
  // otherweise we do a cast and return
  FireModel *result = (FireModel *)mc;
  return result;
}

///------------------------------------------------------------------------//


Simple_MeshIter &
TclHeatTransferModule::getHTMeshs()
{
  theHTMeshIter->reset();
  return *theHTMeshIter;
}


HeatTransferDomain*
TclHeatTransferModule::getHeatTransferDomain(){
  
  return theHTDomain;
}

void 
TclHeatTransferModule::clearAll()
{
    if(theHTDomain!=0)
        theHTDomain->clearAll();

    if (theHTMaterials != 0)
        delete theHTMaterials;

    if (theHTEntities != 0)
        delete theHTEntities;

    if (theHTCons != 0) {
        theHTCons->clearAll();
        delete theHTCons;
    }

    if (theHTEleSets != 0) {
        theHTEleSets->clearAll();
        delete theHTEleSets;
    }

    if (theHTNodeSets != 0) {
        theHTNodeSets->clearAll();
        delete theHTNodeSets;
    }

    if (theFireModels != 0) {
        theFireModels->clearAll();
        delete theFireModels;
    }

    if (theHTMeshes != 0) {
        theHTMeshes->clearAll();
        delete theHTMeshes;
    }

    if (theHTMeshIter != 0)
        delete theHTMeshIter;


    theHTDomain = 0;
    theTclHTModule = 0;
    theTclHTPattern = 0;
    theHTMaterials = 0;
    theHTEntities = 0;
    theHTCons = 0;
    theHTEleSets = 0;
    theHTNodeSets = 0;
    theFireModels = 0;
    theHTMeshes = 0;
    theHTMeshIter = 0;
    if (theHTAnalysis != 0) {
        theHTAnalysis->clearAll();
    }

    //Analysis classes
    theAnalysisModel = 0;
    theAlgorithm = 0;
    theHandler = 0;
    theNumberer = 0;
    theSOE = 0;
    theHTAnalysis = 0;
    theTransientIntegrator = 0;
    theTest = 0;
    HTReorderTag = 0;
}
//-------------------------------------------------------------------------------//
//------------------------Procedure of tcl commmands-----------------------------//

int
TclHeatTransferCommand_addHTMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}

	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
	HeatTransferMaterial* theHTMaterial=0;
	int HTMaterialTag = 0;

	if(argc<3){
	opserr << "WARNING:: Creating HeatTransfer material requires at least 3 arguments." << "\n";	
	return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &HTMaterialTag) != TCL_OK) {
	opserr << "WARNING:: invalid material tag for defining HeatTransfer material: " << argv[1] << "\n";	
	return TCL_ERROR;
    }

    //Adding CarbonSteelEC3
	if (strcmp(argv[1],"SteelEC3") == 0||strcmp(argv[1], "CarbonSteelEC3")== 0) {

		theHTMaterial = new CarbonSteelEC3(HTMaterialTag);

	}
	//Adding ConcreteEC2
	else if(strcmp(argv[1],"ConcreteEC2") == 0){
	  
		double moisture=0;
		bool isLower =false;

		if(argc==3){
			moisture=0;
		}
		else if(argc==4){
			if (Tcl_GetDouble (interp, argv[3], &moisture) != TCL_OK) {
			opserr << "WARNING invalid mositure" << endln;
			opserr << " for HeatTransfer material: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		 }
		else if(argc==5){
			if (Tcl_GetDouble (interp, argv[3], &moisture) != TCL_OK) {
			opserr << "WARNING invalid mositure" << endln;
			opserr << " for HeatTransfer material: " << argv[1] << endln;	    
			return TCL_ERROR;
			}

			if(strcmp(argv[4],"Lower") == 0||strcmp(argv[4],"lower") == 0||strcmp(argv[4],"-lower") == 0){
				isLower = true;
			}
		 }
		else
			opserr<< "WARNING:: Defining HeatTransfer material: "<<argv[1]<<" recieved more than 4 arguments." << "\n";

		theHTMaterial = new ConcreteEC2(HTMaterialTag, moisture, isLower);
  }
  else if(strcmp(argv[1],"SFRMCoating") == 0||strcmp(argv[1],"SFRM") == 0){
    
    int typeTag =0;
    
    if(argc==3){
      typeTag = 1;
    }
    else if(argc==4){
      if (Tcl_GetInt (interp, argv[3], &typeTag) != TCL_OK) {
        opserr << "WARNING invalid typeTag" << endln;
        opserr << " for HeatTransfer material: " << argv[1] << endln;
        return TCL_ERROR;
      }
    }
    else
      opserr<< "WARNING:: Defining HeatTransfer material: "<<argv[1]<<" recieved more than 4 arguments." << "\n";
    
    theHTMaterial = new SFRMCoating(HTMaterialTag, typeTag);

	}
    //add timber HT material
  else if (strcmp(argv[1], "Timber") == 0 || strcmp(argv[1], "timber") == 0) {

       int typeTag = 0;
       Vector Pars = 0;
       Matrix thePars = Matrix();
       const char* Parsfilename = 0;

        if (argc == 3) {
            typeTag = 1;
        }
        else if (argc >3) {
            if (Tcl_GetInt(interp, argv[3], &typeTag) != TCL_OK) {
                opserr << "WARNING invalid typeTag" << endln;
                opserr << " for HeatTransfer material: " << argv[1] << endln;
                return TCL_ERROR;
            }
            int count =4;
            if (typeTag != 0) {
                //Non-EC timber material model, reading parameters from external file
                if (strcmp(argv[count], "-file") == 0 || strcmp(argv[count], "-File") == 0 || strcmp(argv[count], "file") == 0) {
                    count++;
                    Parsfilename = argv[count];


                    //-------------------------------------------------------
                    // determine the number of data points
                    int numDataPoints = 0;
                    int numRows = 0;
                    double dataPoint;
                    ifstream theFile;

                    // first open and go through file containg path
                    theFile.open(Parsfilename, ios::in);
                    if (theFile.bad() || !theFile.is_open()) {
                        opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
                        opserr << " - could not open file " << Parsfilename << endln;
                    }
                    else {
                        while (theFile >> dataPoint)
                            numDataPoints++;
                    }

                    theFile.close();

                    numRows = numDataPoints / 4;
                    // check number of data entries in both are the same

                    if (numDataPoints != 0) {
                        // now create the two vector
                        thePars.resize(numRows, 4);

                        // ensure did not run out of memory creating copies
                        if (thePars.noRows() == 0) {
                            opserr << "WARNING UserDefinedFire::UserDefinedFire() - out of memory\n ";
                        }

                        // first open the file for temperature/flux and read in the data

                        theFile.open(Parsfilename, ios::in);
                        // read in the path data and then do the time
                        int count = 0;
                        while (theFile >> dataPoint) {
                            thePars(count, 0) = dataPoint;
                            theFile >> dataPoint;
                            thePars(count, 1) = dataPoint;
                            theFile >> dataPoint;
                            thePars(count, 2) = dataPoint;
                            theFile >> dataPoint;
                            thePars(count, 3) = dataPoint;
                            count++;
                        }

                        // finally close the file
                        theFile.close();
                    }
#ifdef _DEBUG
                    opserr << "Timber Material properties: " << thePars << endln;
#endif // _DEBUG

                    //----End of reading parameters--------------------------------------------------

                    count++;
                }

            }
            
            
            
            if ((argc - count) > 0) {
                Pars.resize(argc - count);
            }
            else {
                opserr << "WARNING:: no parameter is defined for Timber HT Material: " << argv[1] << "\n";
                return TCL_ERROR;
            }


            //-----for geting uncertain number of double data.
            int ArgStart = count;
            int ArgEnd = argc;
            double data;

            if (ArgStart != ArgEnd) {
                for (int i = ArgStart; i < ArgEnd; i++) {
                    Tcl_GetDouble(interp, argv[i], &data);
                    Pars(i - ArgStart) = data;
                }
            }
#ifdef _DEBUG
            opserr << "HTTimber Material" << argv[1] << ": Parameters " << Pars << endln;
#endif

        }
        else
            opserr << "WARNING:: Defining HeatTransfer material: " << argv[1] << " recieved more than 4 arguments." << "\n";
        if(Parsfilename ==0)
            theHTMaterial = new TimberHTMaterial(HTMaterialTag, typeTag, theHTDomain, Pars);
        else
            theHTMaterial = new TimberHTMaterial(HTMaterialTag, typeTag, theHTDomain, thePars, Pars);

    }
  else if(strcmp(argv[1],"GenericMaterial") == 0){
	  
		double density=0;
		double cp =0;
		double conduct=0;

		if(argc==3){
			opserr << "Default parameters are used for generic material" << endln;

		}
		else if(argc==6){
			if (Tcl_GetDouble (interp, argv[3], &density) != TCL_OK) {
			opserr << "WARNING invalid density" << endln;
			opserr << " for HeatTransfer material: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble (interp, argv[4], &cp) != TCL_OK) {
			opserr << "WARNING invalid specific heat" << endln;
			opserr << " for HeatTransfer material: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble (interp, argv[5], &conduct) != TCL_OK) {
			opserr << "WARNING invalid conductivity" << endln;
			opserr << " for HeatTransfer material: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		 }
		else
			opserr<< "WARNING:: Defining HeatTransfer material: "<<argv[1]<<" recieved more than 4 arguments." << "\n";
        //SimpleMaterial(int tag, double rho, double cp, double kc); 
		theHTMaterial = new SimpleMaterial(HTMaterialTag, density,cp,conduct);

	}

	if(theHTMaterial!=0){
		theTclHTModule->addHTMaterial(*theHTMaterial);
	}
	else
		opserr<<"WARNING: TclHTModule fail to add HeatTransfer Material: "<<argv[1]<<endln;

	return TCL_OK;	

}

//Adding simpleEntity
int
TclHeatTransferCommand_addHTEntity(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}

	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
	Simple_Entity* theHTEntity=0;
	int HTEntityTag = 0;

	if(argc<5){
	opserr << "WARNING:: Creating an entity requires at least 7 arguments." << "\n";	
	return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &HTEntityTag) != TCL_OK) {
	opserr << "WARNING:: invalid entity tag for defining HeatTransfer entity: " << argv[1] << "\n";	
	return TCL_ERROR;
    }

    //Adding 2D entity:Block
  if (strcmp(argv[1],"Line") == 0||strcmp(argv[1],"Line1D") == 0||strcmp(argv[1],"Line1d") == 0)
  {
    double centerX=0.0; double lengthX =0.0;
    if (Tcl_GetDouble (interp, argv[3], &centerX) != TCL_OK) {
      opserr << "WARNING invalid center" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble (interp, argv[4], &lengthX) != TCL_OK) {
      opserr << "WARNING invalid Length" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
    //Simple_Line(int tag, double centerX, double lengthX);
    theHTEntity = new Simple_Line(HTEntityTag, centerX,lengthX);
    
  }
	else if (strcmp(argv[1],"Block") == 0||strcmp(argv[1],"Block2D") == 0||strcmp(argv[1],"Block2d") == 0)
  {
      int ArgStart = 3;
      int ArgEnd = argc;
      if (argc - 3 == 4) {
          double centerX = 0.0; double centerY = 0.0; double BreadthX = 0.0; double HeightY = 0.0;

          if (Tcl_GetDouble(interp, argv[3], &centerX) != TCL_OK) {
              opserr << "WARNING invalid centerX" << endln;
              opserr << " for HeatTransfer entity: " << argv[1] << endln;
              return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[4], &centerY) != TCL_OK) {
              opserr << "WARNING invalid centerY" << endln;
              opserr << " for HeatTransfer entity: " << argv[1] << endln;
              return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[5], &BreadthX) != TCL_OK) {
              opserr << "WARNING invalid centerY" << endln;
              opserr << " for HeatTransfer entity: " << argv[1] << endln;
              return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[6], &HeightY) != TCL_OK) {
              opserr << "WARNING invalid centerY" << endln;
              opserr << " for HeatTransfer entity: " << argv[1] << endln;
              return TCL_ERROR;
          }

          //Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
          theHTEntity = new Simple_Block(HTEntityTag, centerX, centerY, BreadthX, HeightY);
      }
      else if (argc - 3 == 8) {
          Vector data(8);
          int ArgStart = 3;
          int ArgEnd = argc;
          double inputdata;

          if (ArgStart != ArgEnd) {
              for (int i = ArgStart; i < ArgEnd; i++) {
                  Tcl_GetDouble(interp, argv[i], &inputdata);
                  data(i - ArgStart) = inputdata;
              }
          }
          theHTEntity = new Simple_Block(HTEntityTag, data(0), data(1), data(2), data(3), data(4), data(5), data(6), data(7));
      }
      
 }
	//Adding 2D entity:Isection
else if(strcmp(argv[1],"Isection") == 0||strcmp(argv[1],"Isection2D") == 0||strcmp(argv[1],"Isection2d") == 0)
  {
	  
		double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB;

		if (Tcl_GetDouble (interp, argv[3], &HTI_centerX) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[4], &HTI_centerY) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		 if (Tcl_GetDouble (interp, argv[5], &HTI_BF) != TCL_OK) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble (interp, argv[6], &HTI_HB) != TCL_OK) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[7], &HTI_Tw) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[8], &HTI_Tf) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
	    
		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Isection(HTEntityTag, HTI_centerX,HTI_centerY,HTI_BF, HTI_Tf,HTI_Tw,HTI_HB-2*HTI_Tf);

	}

	else if(strcmp(argv[1],"ProtectedIsection") == 0||strcmp(argv[1],"CoatIsection2D") == 0||strcmp(argv[1],"ProtectedIsection2d") == 0)
  {
	  
		double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB ,HTI_coat;

		if (Tcl_GetDouble (interp, argv[3], &HTI_centerX) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[4], &HTI_centerY) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		 if (Tcl_GetDouble (interp, argv[5], &HTI_BF) != TCL_OK) {
			opserr << "WARNING invalid HTI_HB" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble (interp, argv[6], &HTI_HB) != TCL_OK) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[7], &HTI_Tw) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[8], &HTI_Tf) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
	    if (Tcl_GetDouble (interp, argv[9], &HTI_coat) != TCL_OK) {
			opserr << "WARNING invalid HTI_coat" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		//Simple_IsecProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw, double HTI_coat)
		theHTEntity = new Simple_IsecProtected(HTEntityTag, HTI_centerX,HTI_centerY,HTI_BF, HTI_Tf,HTI_Tw,HTI_HB-2*HTI_Tf, HTI_coat );

	}

	//adding BRICK3D
	else if(strcmp(argv[1],"Brick") == 0||strcmp(argv[1],"Brick3D") == 0||strcmp(argv[1],"Brick3d") == 0)
  {
	  
		//Simple_Brick(int tag, double CenterX, double CenterY, double CenterZ, double Breadth_X, double Height_Y,double Length_Z)
		double CenterX, CenterY,CenterZ, Breadth_X,Height_Y,Length_Z;

		if (Tcl_GetDouble (interp, argv[3], &CenterX) != TCL_OK) {
			opserr << "WARNING invalid centerX" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[4], &CenterY) != TCL_OK) {
			opserr << "WARNING invalid centerY" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[5], &CenterZ) != TCL_OK) {
			opserr << "WARNING invalid centerZ" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[6], &Breadth_X) != TCL_OK) {
			opserr << "WARNING invalid Breadth_X" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[7], &Height_Y) != TCL_OK) {
			opserr << "WARNING invalid Height_Y" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[8], &Length_Z) != TCL_OK) {
			opserr << "WARNING invalid Length_Z" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Brick(HTEntityTag, CenterX, CenterY,  CenterZ, Breadth_X, Height_Y, Length_Z);

	}
	else if(strcmp(argv[1],"Isection3D") == 0||strcmp(argv[1],"Isection3d") == 0)
  {
	  
		//Simple_Isection3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
                              //double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw, double HTI_Len);
		double HTI_centerX,HTI_centerY, HTI_centerZ, HTI_BF, HTI_Tf, HTI_HB, HTI_Tw, HTI_Len;
        if (argc < 11) {
            opserr << "WARNING insufficient parameters" << endln;
            opserr << " for HeatTransfer entity: " << argv[1] << endln;
            return TCL_ERROR;
        }
          
		if (Tcl_GetDouble (interp, argv[3], &HTI_centerX) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerX" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[4], &HTI_centerY) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerY" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[5], &HTI_centerZ) != TCL_OK) {
			opserr << "WARNING invalid HTI_centerZ" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
   if (Tcl_GetDouble (interp, argv[6], &HTI_BF) != TCL_OK) {
			opserr << "WARNING invalid HTI_BF" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		if (Tcl_GetDouble (interp, argv[7], &HTI_HB) != TCL_OK) {
      opserr << "WARNING invalid HTI_HB" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
		
		if (Tcl_GetDouble(interp, argv[8], &HTI_Tw) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tw" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble (interp, argv[9], &HTI_Tf) != TCL_OK) {
			opserr << "WARNING invalid HTI_Tf" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		
		
		if (Tcl_GetDouble (interp, argv[10], &HTI_Len) != TCL_OK) {
			opserr << "WARNING invalid HTI_Len" << endln;
			opserr << " for HeatTransfer entity: " << argv[1] << endln;	    
			return TCL_ERROR;
			}
		//Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
		theHTEntity = new Simple_Isection3D(HTEntityTag, HTI_centerX, HTI_centerY, HTI_centerZ,
                             HTI_BF, HTI_Tf, HTI_HB-2*HTI_Tf, HTI_Tw,  HTI_Len);

	}
  //Adding 2D entity:Composite Section
  else if(strcmp(argv[1],"Composite2D") == 0||strcmp(argv[1],"composite2D") == 0||strcmp(argv[1],"Composite2d") == 0)
  {
    
    double HTI_centerX, HTI_centerY, HTI_BF, HTI_Tf, HTI_Tw, HTI_HB, Slab_H,Slab_W;
    
    if (Tcl_GetDouble (interp, argv[3], &HTI_centerX) != TCL_OK) {
      opserr << "WARNING invalid HTI_centerX" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble (interp, argv[4], &HTI_centerY) != TCL_OK) {
      opserr << "WARNING invalid HTI_centerY" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
	if (Tcl_GetDouble(interp, argv[5], &HTI_BF) != TCL_OK) {
		opserr << "WARNING invalid HTI_BF" << endln;
		opserr << " for HeatTransfer entity: " << argv[1] << endln;
		return TCL_ERROR;
	}
    if (Tcl_GetDouble (interp, argv[6], &HTI_HB) != TCL_OK) {
      opserr << "WARNING invalid HTI_HB" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
	if (Tcl_GetDouble(interp, argv[7], &HTI_Tw) != TCL_OK) {
		opserr << "WARNING invalid HTI_Tw" << endln;
		opserr << " for HeatTransfer entity: " << argv[1] << endln;
		return TCL_ERROR;
	}
    
    if (Tcl_GetDouble (interp, argv[8], &HTI_Tf) != TCL_OK) {
      opserr << "WARNING invalid HTI_Tf" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
  
    if (Tcl_GetDouble (interp, argv[9], &Slab_W) != TCL_OK) {
      opserr << "WARNING invalid Slab_w" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble (interp, argv[10], &Slab_H) != TCL_OK) {
      opserr << "WARNING invalid Slab_H" << endln;
      opserr << " for HeatTransfer entity: " << argv[1] << endln;
      return TCL_ERROR;
    }
    
    //Simple_Composite2D(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf,
             //double HTI_Tf, double HTI_Tw, double HTI_Hw,double slabW, double slabH)
    theHTEntity = new Simple_Composite2D(HTEntityTag, HTI_centerX,HTI_centerY,HTI_BF, HTI_Tf,HTI_Tw,HTI_HB-2*HTI_Tf,Slab_W,Slab_H);
    
  }
  
  
	if(theHTEntity!=0){
		theTclHTModule->addHTEntity(*theHTEntity);
	}
	else
		opserr<<"WARNING: TclHTModule fail to add HeatTransfer Material: "<<argv[1]<<endln;

	return TCL_OK;	

}

//Adding Simple Mesh
int
TclHeatTransferCommand_addHTMesh(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}

	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
	int HTMeshTag = 0;
	int HTEntityTag=0;
  
	int HTMaterialTag =0;
	int HTMaterialTag1 =0;
	int PhaseChangeTag= 0;
	int PhaseChangeTag1= 0;
    bool numCtrl = false;
  Vector* SectionLocs = new Vector(0);
  
	Simple_Entity* theHTEntity=0;
	Simple_Mesh* theHTMesh=0;
	HeatTransferMaterial* theHTMaterial=0;
	HeatTransferMaterial* theHTMaterial1=0;
	Vector MeshCtrls =0;
	

	/*if(argc<3){
	opserr << "WARNING:: Creating HeatTransfer material requires at least 3 arguments." << "\n";	
	return TCL_ERROR;
	}*/
  
  int count = 1;

	if (Tcl_GetInt(interp, argv[count], &HTMeshTag) != TCL_OK) {
	opserr << "WARNING:: invalid mesh tag for defining simple mesh: " << argv[1] << "\n";	
	return TCL_ERROR;
    }
  
  count++;

    if (Tcl_GetInt(interp, argv[count], &HTEntityTag) != TCL_OK) {
	opserr << "WARNING:: invalid entity tag for defining simple mesh: " << argv[1] << "\n";	
	return TCL_ERROR;
    }
  
  count++;

	 if (Tcl_GetInt(interp, argv[count], &HTMaterialTag) != TCL_OK) {
	opserr << "WARNING:: invalid HT material tag for defining simple mesh: " << argv[1] << "\n";	
	return TCL_ERROR;
    }

  count++;
  
  //if second material tag is detected
  if(strcmp(argv[count],"-SecondMat") == 0||strcmp(argv[count],"-secondmaterial") == 0){
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTMaterialTag1) != TCL_OK) {
      opserr << "WARNING:: invalid HT material tag for defining simple mesh: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
  }
  
  //if phasechange tag is detected
  if(strcmp(argv[count],"-PhaseChange") == 0||strcmp(argv[count],"-phasechange") == 0
	  ||strcmp(argv[count],"-phaseChange") == 0){
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &PhaseChangeTag) != TCL_OK) {
      opserr << "WARNING:: invalid phasechange tag for defining simple mesh: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    //if second material is defined
    if (HTMaterialTag1!=0) {
      if (Tcl_GetInt(interp, argv[count], &PhaseChangeTag1) != TCL_OK) {
        opserr << "WARNING:: invalid Phasechange tag for defining simple mesh: " << argv[1] << "\n";
        return TCL_ERROR;
      }else{
        count++;
      }
    }
    
  }
  
  //if SectionLoc identified
  if(strcmp(argv[count],"-SectionLoc") == 0||strcmp(argv[count],"-sectionLoc") == 0 ||strcmp(argv[count],"SectionLoc") == 0)
  {
    count++;
    double Loc1, Loc2;
    
    if (Tcl_GetDouble(interp, argv[count], &Loc1) != TCL_OK) {
      opserr << "WARNING:: invalid section location for defining simple mesh: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    SectionLocs = new Vector(1);
    (*SectionLocs)(0)=Loc1;
    //if second material is defined
    if (Tcl_GetDouble(interp, argv[count], &Loc2) == TCL_OK) {
      count++;
      SectionLocs->resize(2);
      (*SectionLocs)(0) = Loc1;
       (*SectionLocs)(1)=Loc2;
    }
    
    
  }
  
  //if meshctrls tag is detected
  if(strcmp(argv[count],"-MeshCtrls") == 0||strcmp(argv[count],"-MeshCtrl") == 0 || strcmp(argv[count], "-NumCtrl") == 0 || strcmp(argv[count], "-numctrl") == 0){
    count++;

    if ((argc - count)>0) {
      MeshCtrls.resize(argc-count);
    } else {
      opserr << "WARNING:: no meshCtrl defined for simple mesh: " << argv[1] << "\n";
      return TCL_ERROR;
    }
//-----for geting uncertain number of double data.
    if (strcmp(argv[count-1], "-MeshCtrls") == 0 || strcmp(argv[count-1], "-MeshCtrl") == 0) {
        int ArgStart = count;
        int ArgEnd = argc;
        double data;

        if (ArgStart != ArgEnd) {
            for (int i = ArgStart; i < ArgEnd; i++) {
                Tcl_GetDouble(interp, argv[i], &data);
                MeshCtrls(i - ArgStart) = data;
            }
        }
    } 
    else if (strcmp(argv[count-1], "-NumCtrl") == 0 || strcmp(argv[count-1], "-numctrl") == 0) {
        numCtrl = true;
        int ArgStart = count;
        int ArgEnd = argc;
        int data;

        if (ArgStart != ArgEnd) {
            for (int i = ArgStart; i < ArgEnd; i++) {
                Tcl_GetInt(interp, argv[i], &data);
                MeshCtrls(i - ArgStart) = data;
            }
        }
    }
    else {
        opserr << "WARNING:TclHTModule- MeshCtrls not found for Mesh " << argv[1] << endln;
        return TCL_ERROR;
    }
    
#ifdef _DEBUG
	  opserr<<"OriginLocs "<< *SectionLocs<<" MeshCtrls "<< MeshCtrls<<endln;
#endif
// for geting uncertain number of doubel values  
  }
  else {
	  opserr<<"WARNING:TclHTModule- MeshCtrls not found for Mesh "<<argv[1]<<endln;
	  return TCL_ERROR;
  }
  
	theHTEntity = theTclHTModule->getHTEntity(HTEntityTag);
	if(theHTEntity==0){
	opserr<< "WARNING:: HTMesh failed to get the requested entity: " <<argv[1] <<"\n";
	return TCL_ERROR;
	}

	theHTMaterial = theTclHTModule->getHTMaterial(HTMaterialTag);
	if(theHTMaterial==0){
	opserr<< "WARNING:: TclHTModule failed to get the requested HT material when defining Simple Mesh: " <<argv[1] <<"\n";
	return TCL_ERROR;
	}

	if (HTMaterialTag1!=0) {
		theHTMaterial1 = theTclHTModule->getHTMaterial(HTMaterialTag1);
		if(theHTEntity==0){
		opserr<< "WARNING:: TclHTModule failed to get the requested HT material when defining Simple Mesh: " <<argv[1] <<"\n";
		return TCL_ERROR;
		}
	}
  
  ID PHaseIDs=0;
  if (PhaseChangeTag1!=0) {
    PHaseIDs.resize(2);
    PHaseIDs(0)=PhaseChangeTag;
    PHaseIDs(1)=PhaseChangeTag1;
  }
  else
  {
    if (PhaseChangeTag!=0) {
      PHaseIDs.resize(1);
      PHaseIDs(0)=PhaseChangeTag;
    }
    
  }

    theHTMesh = new Simple_Mesh(HTMeshTag,theHTEntity,theHTDomain,theHTMaterial,MeshCtrls,theHTMaterial1,numCtrl);


	if(theHTMesh!=0){
    if (PHaseIDs!=0) {
      theHTMesh->SetEleParameters(PHaseIDs);
    }
	if(SectionLocs->Size()!=0)
	theHTMesh->SetOriginLocs(*SectionLocs);
 
	theTclHTModule->addHTMesh(theHTMesh);
	}
	else
		opserr<<"WARNING: TclHTModule fail to add Simple Mesh: "<<argv[1]<<endln;


	return TCL_OK;	

}


//TclHeatTransferCommand_HTRefineMesh
int
TclHeatTransferCommand_HTRefineMesh(ClientData clientData, Tcl_Interp *interp, int argc,
                                 TCL_Char **argv)
{
    if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}

	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
    int count=1;
    int HTEntityTag = 0;
    Simple_Entity* theHTEntity = 0;
	int SeedTag = 0 ;
  int SeedTag1 = 0;
  Vector MeshVec = 0;
	//if HTEntity tag is detected
  if(strcmp(argv[count],"-HTEntity") == 0||strcmp(argv[count],"-Entity") == 0||strcmp(argv[count],"HTEntity") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTEntityTag) != TCL_OK) {
      opserr << "WARNING:: invalid HT Entity tag for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
     theHTEntity = theTclHTModule->getHTEntity(HTEntityTag);
    if(theHTEntity==0){
      opserr<< "WARNING:: TclHTModule failed to get the requested entity when defining simple mesh: " <<argv[1] <<"\n";
      return TCL_ERROR;
    }
    
  }
  
  
  //Entity obtained
  
  if(strcmp(argv[count],"-SeedTag") == 0||strcmp(argv[count],"-seedTag") == 0||strcmp(argv[count],"SeedTag") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &SeedTag) != TCL_OK) {
      opserr << "WARNING:: invalid HT Entity tag for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    if (Tcl_GetInt(interp, argv[count], &SeedTag1) == TCL_OK) {
      count++;
      
    }
    else{
      SeedTag1 = 0;
    }
    
    //obtaining second seed tag if necessary
    
  }
  //Seed tag
  
  if(strcmp(argv[count],"-space") == 0||strcmp(argv[count],"-Space") == 0||strcmp(argv[count],"Space") == 0){
    
    count++;
    if ((argc - count)>0) {
      MeshVec.resize(argc-count);
    } else {
      opserr << "WARNING:: no Seed info defined for refining mesh: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    //-----for geting uncertain number of double data.
    int ArgStart = count;
    int ArgEnd = argc;
    double data;
    /*while (count < argc && ArgEnd == 0) {
     if (Tcl_GetDouble(interp, argv[count], &data) != TCL_OK)
     ArgEnd = count;
     else
     count++;
     }*/
    if (ArgStart != ArgEnd) {
      for (int i=ArgStart; i<ArgEnd; i++) {
        Tcl_GetDouble(interp, argv[i], &data);
        MeshVec(i-ArgStart) = data;
      }
    }
    // for geting uncertain number of doubel values
    

  }
  
  //space info obtained
  if(SeedTag!=0){
    if(theHTEntity->RefineSeeds ( SeedTag, MeshVec)!=0){
      opserr<<"WARNING::the Entity "<<HTEntityTag<< "failed to refine Seeds" << HTEntityTag << endln;
    }
  }
  if(SeedTag1!=0){
    if(theHTEntity->RefineSeeds ( SeedTag1, MeshVec)!=0){
      opserr<<"WARNING::the Entity "<<HTEntityTag<< "failed to refine Seeds "<< SeedTag1 << endln;
    }
  }

  #ifdef _DEBUG
  opserr << "HTEntity " << HTEntityTag << " has refined mesh as:" << MeshVec << endln;
  #endif

  return TCL_OK;
 
////////////////////////////////////////////
}

//MeshAll
int
TclHeatTransferCommand_HTMeshAll(ClientData clientData, Tcl_Interp *interp, int argc,
                                 TCL_Char **argv)
{
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}

	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
    }

  Simple_Mesh *theHTMesh;
  Simple_MeshIter &theHTMeshes  = theTclHTModule->getHTMeshs();
  while((theHTMesh = theHTMeshes())!=0){
    theHTMesh->GeneratingNodes();
    theHTMesh->GeneratingEles();
  }

return TCL_OK;

}


//SetInitialT
int
TclHeatTransferCommand_SetInitialT(ClientData clientData, Tcl_Interp *interp, int argc,
                                 TCL_Char **argv)
{
  
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
  }
  
  double initialT;
  
  if (Tcl_GetDouble (interp, argv[1], &initialT) != TCL_OK) {
    opserr << "WARNING invalid initial Temperature" << endln;
    opserr << " for TclHeatTransferModule: setInitialT " << argv[1] << endln;
    return TCL_ERROR;
  }
  
  theHTDomain->setInitial(initialT);
  
  return TCL_OK;
  
}


//HTConstants
int
TclHeatTransferCommand_addHTConstants(ClientData clientData, Tcl_Interp *interp, int argc,
                                   TCL_Char **argv)
{

   //Recieving heat transfer coefficients: convective coefficient, ambient temperature, absorption radiative coefficient,
    //                                     steffan-B constant ,fire radiative cofficient, irradiation
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTMaterial\n";
    return TCL_ERROR;
  }
  
  HTConstants* theHTConstants=0;
	int HTConstantTag = 0;
  
  if (Tcl_GetInt(interp, argv[1], &HTConstantTag) != TCL_OK) {
    opserr << "WARNING:: invalid constants tag for defining HT constants: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  
  int count=2;
  Vector constants(5);
  if((argc-count!=4)&&(argc-count!=5)){
	  opserr<<"WARNING:: Non-sufficent or oversized data for defining heat transfer contants: "<<argv[1]<<"\n";
	  return TCL_ERROR;
  }
  //-----for geting uncertain number of double data.
    int ArgStart = count;
      int ArgEnd = argc;
      double data; 
      /*while (count < argc && ArgEnd == 0) {
		if (Tcl_GetDouble(interp, argv[count], &data) != TCL_OK)
		ArgEnd = count;
		else
		 count++;
	}*/
      if (ArgStart != ArgEnd) {
		for (int i=ArgStart; i<ArgEnd; i++) {
		Tcl_GetDouble(interp, argv[i], &data);
		constants(i-ArgStart) = data;
		}
      }
// for geting uncertain number of doubel values 
  if(argc-count==4){
	  constants(4)= 5.67e-8;

	  //calculating radiation from the ambient air;
  }

#ifdef _DEBUG
  opserr << "TclHeatTransfer:: AddHTconstants: " << constants << endln;
#endif

  theHTConstants = new HTConstants(HTConstantTag,constants);
  
  if(theHTConstants!=0){
		theTclHTModule->addHTConstants(theHTConstants);
	}
	else
		opserr<<"WARNING: TclHTModule fail to add HTConstants: "<<argv[1]<<endln;
  
	return TCL_OK;

}



//HTEleSet
int
TclHeatTransferCommand_HTEleSet(ClientData clientData, Tcl_Interp *interp, int argc,
                                   TCL_Char **argv)
{
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTEleSet\n";
    return TCL_ERROR;
  }
  
  if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTNodeSet\n";
    return TCL_ERROR;
   }
  
  HTEleSet* theHTEleSet=0;
	int HTEleSetTag = 0;
  int HTEntityTag = 0;
  int FaceID =0;
  
  if (Tcl_GetInt(interp, argv[1], &HTEleSetTag) != TCL_OK) {
    opserr << "WARNING:: HTEleSet failed to identify set tag: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  
  int count=2;
  
  //if HTEntity tag is detected
  if(strcmp(argv[count],"-HTEntity") == 0||strcmp(argv[count],"-Entity") == 0||strcmp(argv[count],"HTEntity") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTEntityTag) != TCL_OK) {
      opserr << "WARNING::HTEleSet failed to identify HTEntity tag: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
  }else{
	opserr<<"WARNING TclHTCommand HTEleSet wants tag -HTEntity, HTEntity, -Entity "<<endln;
  }
  
  
  Simple_Entity* theHTEntity = theTclHTModule->getHTEntity(HTEntityTag);
  if(theHTEntity==0){
    opserr<< "WARNING:: HTEleSet failed to get the HTEntity: " <<argv[1] <<"\n";
    return TCL_ERROR;
  }
  
  int EntityMeshTag = theHTEntity->getMeshTag();
  Simple_Mesh* theHTMesh = theTclHTModule->getHTMesh(EntityMeshTag);
  
  
  ID EleRange(0);
  
  //if face tag is detected
  if(strcmp(argv[count],"-face") == 0||strcmp(argv[count],"-Face") == 0||strcmp(argv[count],"face") == 0)
  {
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &FaceID) != TCL_OK) {
      opserr << "WARNING:: invalid face tag for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    int eleFaceID;
    theHTMesh->SelectingElesbyFace(EleRange, FaceID,eleFaceID);
  }
  else if(strcmp(argv[count],"-NodeSet") == 0||strcmp(argv[count],"-nodest") == 0||strcmp(argv[count],"NodeSet") == 0)
  {
    
    count++;
    
    int NodeSetID;
    if (Tcl_GetInt(interp, argv[count], &NodeSetID) != TCL_OK) {
      opserr << "WARNING:: invalid NodesetID for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    
    HTNodeSet* theHTNodeSet = theTclHTModule->getHTNodeSet(NodeSetID);
    
    if(theHTNodeSet ==0){
      opserr<< "WARNING:: invalid NodesetID for defining HTEleSet: "<< argv[1] << "\n";
      return TCL_ERROR;
    }
    
    ID NodesRange = theHTNodeSet->getNodeID();
    
    if(strcmp(argv[count],"-face") == 0||strcmp(argv[count],"-Face") == 0||strcmp(argv[count],"face") == 0){
      count++;
      
      if (Tcl_GetInt(interp, argv[count], &FaceID) != TCL_OK) {
      opserr << "WARNING:: invalid face tag for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
	  }else{
      count++;
      }
	}

    theHTMesh->SelectingEles(EleRange, NodesRange, FaceID);
    
  }
  
    
    
    
    
    
  theHTEleSet = new HTEleSet(HTEleSetTag);
  
  if (EleRange!=0) {
    theHTEleSet->addEleID(EleRange);
  }

  opserr << "HTEleSet " << HTEleSetTag << " selects elements:" << EleRange << endln;
  
  if(theHTEleSet!=0){
		theTclHTModule->addHTEleSet(theHTEleSet);
	}
	else
		opserr<<"WARNING: TclHTModule failed to add HTEleSet: "<<argv[1]<<endln;
  
	return TCL_OK;
  
}

//--------------------
//HTNodeSet
int
TclHeatTransferCommand_HTNodeSet(ClientData clientData, Tcl_Interp *interp, int argc,
                                   TCL_Char **argv)
{
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTNodeSet\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTNodeSet\n";
    return TCL_ERROR;
  }
  
  HTNodeSet* theHTNodeSet = 0;
  Simple_Entity* theHTEntity = 0;
  Simple_Mesh* theHTMesh = 0;
  int HTNodeSetTag = 0;
  int HTEntityTag = 0;
  int FaceID =0;
  ID NodeRange(0);
  
  if (Tcl_GetInt(interp, argv[1], &HTNodeSetTag) != TCL_OK) {
    opserr << "WARNING:: invalid entity tag for defining HTNodeSet: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  
  int count=2;

  //if directly giving the node tags
if (strcmp(argv[count], "-HTNodes") == 0 ||strcmp(argv[count], "-nodes") == 0 ||strcmp(argv[count], "-HTNodeSet") == 0 || strcmp(argv[count], "-NodeSet") == 0 || strcmp(argv[count], "nodeset") == 0)
{
    count++;

    //-----for geting uncertain number of integer data.
    int ArgStart = count;
    int ArgEnd = 0;
    int data;
    ID NodeSetRange(0);

    while (count < argc && ArgEnd == 0) {
        if (count == argc - 1)
            ArgEnd = count + 1;
        else if (Tcl_GetInt(interp, argv[count], &data) != TCL_OK)
        {
            opserr << "Error creating NodeSet " << argv[1] << endln;
            opserr << "\"" << argv[count] << "\" detected after using -NodeSet flag." << endln;
            opserr << "Can only have integers after the -NodeSet flag when creating a new combined nodeset." << endln;
            return TCL_ERROR;
        }
        else
            count++;
    }
    if (ArgEnd - ArgStart <= 0) {
        opserr << "Error creating NodeSet " << argv[1] << endln;
        opserr << "Not enough arguments to create a combined of NodeSets. Need at least 1 valid component NodeSet." << endln;
        return TCL_ERROR;
    }
    //~ detecting the remaining number of input
    NodeSetRange.resize(ArgEnd - ArgStart);
    if (ArgStart != ArgEnd) {
        for (int i = ArgStart; i < ArgEnd; i++) {
            Tcl_GetInt(interp, argv[i], &data);
            NodeSetRange(i - ArgStart) = data;
        }
    }
    vector <int> NodeRangeV;
    for (int i = 0; i < NodeSetRange.Size(); i++) {
        HTNodeSet* theNodeSet = theTclHTModule->getHTNodeSet(NodeSetRange(i));
        if (theNodeSet == 0)
        {
            opserr << "Cannot create combined NodeSet " << argv[1] << endln;
            opserr << "Component NodeSet at position " << i + 1 << " with ID " << NodeSetRange(i) << " not found." << endln;
            return TCL_ERROR;
        }
        ID tempNodeSetID = theNodeSet->getNodeID();
        if (tempNodeSetID == 0) {
            opserr << "Cannot create combined NodeSet " << argv[1] << endln;
            opserr << "Component NodeSet at position " << i + 1 << " with ID " << NodeSetRange(i) << " is empty." << endln;
            return TCL_ERROR;
        }
        for (int j = 0; j < tempNodeSetID.Size(); j++) {
            NodeRangeV.push_back(tempNodeSetID(j));
        }
    }

    NodeRange.resize(NodeRangeV.size());
    for (int i = 0; i < NodeRangeV.size(); i++)
        NodeRange(i) = NodeRangeV[i];

}
else
{
    //if HTEntity tag is detected, then we definitely expect faceID
    if (strcmp(argv[count], "-HTEntity") == 0 || strcmp(argv[count], "-Entity") == 0 || strcmp(argv[count], "HTEntity") == 0) {

        count++;

        if (Tcl_GetInt(interp, argv[count], &HTEntityTag) != TCL_OK) {
            opserr << "WARNING:: invalid HT Entity tag for defining HTNodeSet: " << argv[1] << "\n";
            return TCL_ERROR;
        }
        else {
            count++;
        }


        theHTEntity = theTclHTModule->getHTEntity(HTEntityTag);
        if (theHTEntity == 0) {
            opserr << "WARNING:: HTNodeSet failed to get the HTEntity: " << argv[1] << "\n";
            return TCL_ERROR;
        }
        //HTEntityTag obtained

        int EntityMeshTag = theHTEntity->getMeshTag();
        theHTMesh = theTclHTModule->getHTMesh(EntityMeshTag);

        //to obtain faceTag
        if (strcmp(argv[count], "-face") == 0 || strcmp(argv[count], "-Face") == 0 || strcmp(argv[count], "face") == 0) {
            count++;
            if (Tcl_GetInt(interp, argv[count], &FaceID) != TCL_OK) {
                opserr << "WARNING:: invalid face tag for defining HTEleSet: " << argv[1] << "\n";
                return TCL_ERROR;
            }
            else {
                count++;
            }

            theHTMesh->SelectingNodesbyFace(NodeRange, FaceID);

        }

    }

  //end of HTEntity tag 

  //Now search node in the range of xloc
  //check whether it is the end of arguments
    if (argc - count > 0) {
        if (strcmp(argv[count], "-Locx") == 0 || strcmp(argv[count], "Locx") == 0 || strcmp(argv[count], "locx") == 0 || strcmp(argv[count], "-locx") == 0) {
            count++;


            double xlocLB, xlocUB;

            if (Tcl_GetDouble(interp, argv[count], &xlocLB) != TCL_OK) {
                opserr << "WARNING invalid xloc" << endln;
                opserr << " for adding node to HTNodeSet: " << endln;
                return TCL_ERROR;
            }
            count++;
            if (argc - count > 0) {
                if (Tcl_GetDouble(interp, argv[count], &xlocUB) == TCL_OK) {
                    count++;
                }
                else
                    xlocUB = xlocLB;
            }
            else
                xlocUB = xlocLB;


            if (xlocUB < xlocLB) {
                opserr << "WARNING::TclHeatTransfer::xlocUB " << xlocUB << "should be greater than xlocLb " << xlocLB;
                return TCL_ERROR;
            }
            if (theHTMesh == 0) {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within the domain, and inside the interval x = (" << xlocLB << ", " << xlocUB << ")." << endln;
#endif  
                theHTDomain->SelectingNodes(NodeRange, 0, xlocLB, xlocUB);
            }
            else {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within mesh, and inside the interval x = (" << xlocLB << ", " << xlocUB << ")." << endln;
#endif
                theHTMesh->SelectingNodes(NodeRange, 0, xlocLB, xlocUB);
            }
            // for geting uncertain number of double values
            if (NodeRange == 0) {
                opserr << "WARNING: There are no nodes to add to the HTNodeSet at interval x = (" << xlocLB << ", " << xlocUB << ")." << endln;
                opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << argv[1] << endln;
                return -1;
            }

        }
    }
    //search node in the range of yloc
    if (argc - count > 0) {
        if (strcmp(argv[count], "-Locy") == 0 || strcmp(argv[count], "Locy") == 0 || strcmp(argv[count], "locy") == 0) {
            count++;



            double ylocLB, ylocUB;

            if (Tcl_GetDouble(interp, argv[count], &ylocLB) != TCL_OK) {
                opserr << "WARNING invalid yloc" << endln;
                opserr << " for adding node to HTNodeSet: " << endln;
                return TCL_ERROR;
            }
            count++;
            if (argc - count > 0) {
                if (Tcl_GetDouble(interp, argv[count], &ylocUB) == TCL_OK) {
                    count++;
                }
                else
                    ylocUB = ylocLB;
            }
            else
                ylocUB = ylocLB;

            if (ylocUB < ylocLB) {
                opserr << "WARNING::TclHeatTransfer::ylocUB " << ylocUB << "should be greater than ylocLb " << ylocLB;
                return TCL_ERROR;
            }
            if (theHTMesh == 0) {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within the domain, and inside the interval y = (" << ylocLB << ", " << ylocUB << ")." << endln;
#endif  
                theHTDomain->SelectingNodes(NodeRange, 1, ylocLB, ylocUB);
            }
            else {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within mesh, and inside the interval y = (" << ylocLB << ", " << ylocUB << ")." << endln;
#endif
                theHTMesh->SelectingNodes(NodeRange, 1, ylocLB, ylocUB);
            }

            // for geting uncertain number of doubel values
            if (NodeRange == 0) {
                opserr << "WARNING: There are no nodes to add to the HTNodeSet at interval y = (" << ylocLB << ", " << ylocUB << ")." << endln;
                opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << argv[1] << endln;
                return -1;
            }
        }
    }
    //search node in the range of zloc
    if (argc - count > 0) {
        if (strcmp(argv[count], "-Locz") == 0 || strcmp(argv[count], "Locz") == 0 || strcmp(argv[count], "locz") == 0) {
            count++;

            double zlocLB, zlocUB;

            if (Tcl_GetDouble(interp, argv[count], &zlocLB) != TCL_OK) {
                opserr << "WARNING invalid zloc" << endln;
                opserr << " for adding node to HTNodeSet: " << endln;
                return TCL_ERROR;
            }
            count++;
            if (argc - count > 0) {
                if (Tcl_GetDouble(interp, argv[count], &zlocUB) == TCL_OK) {
                    count++;
                }
                else
                    zlocUB = zlocLB;
            }
            else
                zlocUB = zlocLB;

            if (zlocUB < zlocLB) {
                opserr << "WARNING::TclHeatTransfer::zlocUB " << zlocUB << "should be greater than zlocLb " << zlocLB;
                return TCL_ERROR;
            }
            if (theHTMesh == 0) {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within the domain, and inside the interval z = (" << zlocLB << ", " << zlocUB << ")." << endln;
#endif  
                theHTDomain->SelectingNodes(NodeRange, 2, zlocLB, zlocUB);
            }
            else {
#ifdef _DEBUG
                opserr << "Nodeset " << argv[1] << " is selecting nodes from within mesh, and inside the interval z = (" << zlocLB << ", " << zlocUB << ")." << endln;
#endif
                theHTMesh->SelectingNodes(NodeRange, 2, zlocLB, zlocUB);
            }

            // for geting uncertain number of double values 
            if (NodeRange == 0) {
                opserr << "WARNING: There are no nodes to add to the HTNodeSet at interval z = (" << zlocLB << ", " << zlocUB << ")." << endln;
                opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << argv[1] << endln;
                return -1;
            }
        }
    }
}
  
  
  if (NodeRange!=0) {
    theHTNodeSet = new HTNodeSet(HTNodeSetTag);
    theHTNodeSet->addNodeID(NodeRange);
  }
  else{
    opserr << "WARNING: There are no nodes to add to the HTNodeSet." << endln;
    opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << argv[1] << endln;
    return -1;
  }
  
  if (theHTNodeSet != 0) {
      theTclHTModule->addHTNodeSet(theHTNodeSet);
      opserr << "NodeSet " << HTNodeSetTag << " selects nodes:" << NodeRange << endln;
      return TCL_OK;
  }
  else {
      opserr << "WARNING: TclHTModule failed to add HTNodeSet: " << argv[1] << endln;
      return -1;
  }
}

//HTPattern
int
TclHeatTransferCommand_addHTPattern(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char **argv)
{
  if (argc<3) {
    opserr<<"WARNING invalid command-want:HTPattern type\n";
    opserr<<" valid types: AmbientBC, fireAction(fire) \n";
    return TCL_ERROR;
  }
  
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTPattern\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - HTPattern\n";
    return TCL_ERROR;
  }
  
  BoundaryPattern* theHTPattern=0;
	int HTPatternTag = 0;
  
  if (Tcl_GetInt(interp, argv[2], &HTPatternTag) != TCL_OK) {
    opserr << "WARNING:: invalid entity tag for defining HTPattern: " << argv[2] << "\n";
    return TCL_ERROR;
  }
  
  int  commandEndMarker = 0;
  int count=2;
  //if face tag is detected
  if(strcmp(argv[1],"AmbientBC") == 0||strcmp(argv[1],"ambientBC") == 0||strcmp(argv[1],"ambient") == 0){
    
    theHTPattern = new BoundaryPattern(HTPatternTag);
    
    count++;
    
  }
  //if fireAction tag is detected here
  else if(strcmp(argv[1],"fireAction") == 0||strcmp(argv[1],"FireAction") == 0||strcmp(argv[1],"fire") == 0){
    
    int FireModelID = 0;
    theHTPattern = new FireImposedPattern(HTPatternTag);
    
    count++;
    
   if(strcmp(argv[count],"model") == 0||strcmp(argv[count],"Model") == 0||strcmp(argv[count],"-model") == 0){
     count++;
     if (Tcl_GetInt(interp, argv[count], &FireModelID) != TCL_OK) {
       opserr << "WARNING:: invalid entity tag for defining HTPattern: " << argv[2] << "\n";
       return TCL_ERROR;
     }
     count++;
   }
    
    FireModel *theFireModel = theTclHTModule->getFireModel(FireModelID);
    if (theFireModel ==0) {
      opserr << "WARNING:: no fire model is defined for heat flux BC: "  << "\n";
      return TCL_ERROR;
    }
   ((FireImposedPattern*)theHTPattern)->setFireModel(theFireModel);
    
  }
  //END OF FIRE ACTION PATTERN
  else {
    opserr<<"WARNING unknown pattern type"<< argv[1];
    opserr<<"- want : HTPattern patternType "<<HTPatternTag;
    opserr<<" valid types: AmbientBC, fireAction(fire) \n";
    return TCL_ERROR;
  }
  
  //Add boundaryPattern into the heat transfer domain
  if(theHTDomain->addBoundaryPattern(theHTPattern)==false) {
    opserr<< "WARNING could not add boundary pattern to the heat transfer domain" << *theHTPattern;
    delete theHTPattern;
    return TCL_ERROR;
  }
  
  theTclHTPattern = theHTPattern;
  
  
  // use TCL_Eval to evaluate the list of load and single point constraint commands
  if (Tcl_Eval(interp, argv[argc-1]) != TCL_OK) {
    opserr << "WARNING - error reading load pattern information in { } ";
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


//Add Fire Model
int
TclHeatTransferCommand_addFireModel(ClientData clientData, Tcl_Interp *interp, int argc,
                                    TCL_Char **argv)
{
  
// checking the TclHTModule exists or not
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - HTPattern\n";
    return TCL_ERROR;
	}
// checking the HTDomain exists or not
 if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain -  HTPattern\n";
    return TCL_ERROR;
  }
  
//create a pointer to base class
  FireModel* theFireModel=0;
	int FireModelTag = 0;

//get the fireModel tag;
  if (Tcl_GetInt(interp, argv[2], &FireModelTag) != TCL_OK) {
    opserr << "WARNING:: invalid fireModel Tag: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  
	int count=2;
  //Remove the type tag to follow OpenSees Command convention
 
   
    //standard fire curve;
    if(strcmp(argv[1],"standard") == 0||strcmp(argv[1],"Standard") == 0){
      theFireModel = new NorminalFireEC1(FireModelTag , 1);
  
    }
    //hydroCarbon fire curve;
    else if(strcmp(argv[1],"hydroCarbon") == 0||strcmp(argv[1],"HydroCarbon") == 0){
      theFireModel = new NorminalFireEC1(FireModelTag , 3); // hydrocarbon fire tag is 3;
    }
    else if (strcmp(argv[1], "external") == 0 || strcmp(argv[1], "External") == 0) {
        theFireModel = new NorminalFireEC1(FireModelTag, 2); // hydrocarbon fire tag is 2;
    }
	else if(strcmp(argv[1],"ASTM") == 0||strcmp(argv[1],"ASTME119") == 0){
      theFireModel = new NorminalFireEC1(FireModelTag , 4); // hydrocarbon fire tag is 3;
    }
	else if (strcmp(argv[1], "Exponent") == 0 || strcmp(argv[1], "Exp") == 0) {
		theFireModel = new NorminalFireEC1(FireModelTag, 5); // hydrocarbon fire tag is 3;
	}
    else if (strcmp(argv[1], "UserDefined") == 0 || strcmp(argv[1], "userDefined") == 0) {
        count++;
        char* filename = 0;
        int dataType = 1;
        if (argc < 3) {
            opserr << "WARNING insufficient input for UserDeifined Fire Model " << FireModelTag<< endln;
            return TCL_ERROR;
        }

        if (strcmp(argv[count], "-file") == 0 || strcmp(argv[1], "File") == 0 ||strcmp(argv[1], "file") == 0) {
            count++;
            filename = argv[count];
            count++;
        }
        if (argc - count > 0) {
            if (strcmp(argv[count], "-type") == 0 || strcmp(argv[count], "dataType") == 0 || strcmp(argv[count], "type") == 0) {
                count++;
                if (Tcl_GetInt(interp, argv[count], &dataType) != TCL_OK) {
                    opserr << "WARNING:: invalid fireModel Tag: " << argv[1] << "\n";
                    return TCL_ERROR;
                }
            }
        }

        theFireModel = new UserDefinedFire(FireModelTag, filename, dataType);

    }
    //Paramtetric fire
	else if(strcmp(argv[1],"parametric") == 0||strcmp(argv[1],"Parametric") == 0){
		double thi=0; double avent=0; double hvent=0; double atotal=0; double afire=0; double qfire=0; double Tlim=0;
		count++;
		if(argc==10){
			if (Tcl_GetDouble(interp, argv[count], &thi) != TCL_OK) {
			opserr << "WARNING invalid thermal inertia of the compartment boundaries" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+1], &avent) != TCL_OK) {
			opserr << "WARNING invalid total area of vertial openings on walls" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+2], &hvent) != TCL_OK) {
			opserr << "WARNING invalid weighted average of window heights on walls" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+3], &atotal) != TCL_OK) {
			opserr << "WARNING invalid total area of the compartment(including walls)" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+4], &afire) != TCL_OK) {
			opserr << "WARNING invalid area of the floor with fire" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+5], &qfire) != TCL_OK) {
			opserr << "WARNING invalid total design fire" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			if (Tcl_GetDouble(interp, argv[count+6], &Tlim) != TCL_OK) {
			opserr << "WARNING invalid time levels corresponds to different fire growth rate" << endln;
			opserr << " for HeatTransfer fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		}
		else
			opserr<< "WARNING:: Defining Parametric fire: "<<argv[2]<<" recieved insufficient arguments" << "\n";

		theFireModel = new ParametricFireEC1(FireModelTag, thi, avent, hvent, atotal, afire, qfire, Tlim);
	}
	//localised SFPE fire curve;
    else if(strcmp(argv[1],"localisedSFPE") == 0||strcmp(argv[1],"LocalisedSFPE") == 0){
		
		count++; //count should be updated
		double crd1=0.0; double crd2=0.0; double crd3=0.0; 
		double D=0; double Q=0; double HB=0; double HC = 0; int lineTag=0;
		
		if(argc-count<=0)
		{
			opserr << "WARNING invalid aguments" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
		}

		//Add a tag for location of origin;
		if(strcmp(argv[count],"-origin") == 0||strcmp(argv[count],"origin") == 0){
			count++;
			if (Tcl_GetDouble(interp, argv[count], &crd1) != TCL_OK) {
			opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &crd2) != TCL_OK) {
			opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			
			if (Tcl_GetDouble(interp, argv[count], &crd3) == TCL_OK) {
				//if the z loc is successfully recieved, count should be added with 1;
			count++;
			}
			else 
			{
				//it it possible not to have a z loc for localised fire definiton
				opserr << "WARNING invalid z axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	
				crd3=0.0;
			}
		}
		//end of fire origin, waiting for firePars;
		if(strcmp(argv[count],"-firePars") == 0||strcmp(argv[count],"firePars") == 0){
			count++;
			
			if (Tcl_GetDouble(interp, argv[count], &D) != TCL_OK) {
			opserr << "WARNING invalid diameter of the fire source" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &Q) != TCL_OK) {
			opserr << "WARNING invalid rate of heat release" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &HC) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
            if (Tcl_GetDouble(interp, argv[count], &HB) != TCL_OK) {
                opserr << "WARNING invalid distance between the fire source and the beam" << endln;
                opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                return TCL_ERROR;
            }
            count++;
			//detect argument for linetag;
			if(argc-count>0){
				if (Tcl_GetInt(interp, argv[count], &lineTag) != TCL_OK) {
				opserr << "WARNING invalid central line tag " << endln;
				opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
				return TCL_ERROR;
				}
			}
			else {
				opserr << "Central line tag for localised fire "<< argv[2]<<"is set as default:3" << endln;
				lineTag=3;
			}
		}
		else {
			opserr<< "WARNING:: Defining Localised fire "<<argv[2]
			      <<" expects tag:-firePars or firePars" << "\n";
		}

		theFireModel = new LocalizedFireSFPE(FireModelTag, crd1, crd2, crd3, D, Q, HC, HB, lineTag);
    } 
    //travelling
    else if (strcmp(argv[1], "Travelling") == 0 || strcmp(argv[1], "travelling") == 0 || strcmp(argv[1], "NaturalFire") == 0) {

    count++; //count should be updated
    double crd1 = 0.0; double crd2 = 0.0; double crd3 = 0.0;
    double D = 0; double Q = 0; double H = 0; int lineTag = 0;
    double smokeT = 0;

    if (argc - count <= 0)
    {
        opserr << "WARNING invalid aguments" << endln;
        opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
        return TCL_ERROR;
    }

    //end of fire origin, waiting for firePars;
    if (strcmp(argv[count], "-firePars") == 0 || strcmp(argv[count], "firePars") == 0) {
        count++;

        if (strcmp(argv[count], "-file") == 0 || strcmp(argv[count], "file") == 0 || strcmp(argv[count], "-File") == 0) {
            //Definition of Fire parameters from external file
            count++;
            const char* fileName = argv[count];
            PathTimeSeriesThermal* theSeries = new PathTimeSeriesThermal(1, fileName,8, false);
            count++;

            //For lineTag if available
            if (argc - count > 0) {
                if (Tcl_GetInt(interp, argv[count], &lineTag) != TCL_OK) {
                    opserr << "WARNING invalid ceiling height" << endln;
                    opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                    return TCL_ERROR;
                }
            }
            else
                lineTag = 2;

            theFireModel = new NaturalFire(FireModelTag, lineTag, theSeries);


        }
        else {
            //Definition of Fire parameters from commandline
            if (Tcl_GetDouble(interp, argv[count], &D) != TCL_OK) {
                opserr << "WARNING invalid diameter of the fire source" << endln;
                opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                return TCL_ERROR;
            }
            count++;
            if (Tcl_GetDouble(interp, argv[count], &Q) != TCL_OK) {
                opserr << "WARNING invalid rate of heat release" << endln;
                opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                return TCL_ERROR;
            }
            count++;
            if (Tcl_GetDouble(interp, argv[count], &H) != TCL_OK) {
                opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
                opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                return TCL_ERROR;
            }
            count++;
            //detect argument for linetag;
            if (argc - count > 0) {
                if (Tcl_GetInt(interp, argv[count], &lineTag) != TCL_OK) {
                    opserr << "WARNING invalid central line tag " << endln;
                    opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                    return TCL_ERROR;
                }
                count++;
            }


            if (argc - count > 0) {
                if (Tcl_GetDouble(interp, argv[count], &smokeT) != TCL_OK) {
                    opserr << "WARNING invalid smoke temperature between the fire source and the ceiling" << endln;
                    opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                    return TCL_ERROR;
                }
            }
            theFireModel = new NaturalFire(FireModelTag, D, Q, H, lineTag, smokeT);

        }
        
    }
    
    }
    //localised EC1
    else if (strcmp(argv[1], "localised") == 0 || strcmp(argv[1], "Localised") == 0 || strcmp(argv[1], "LocalisedEC") == 0) {

    count++; //count should be updated
    double crd1 = 0.0; double crd2 = 0.0; double crd3 = 0.0;
    double D = 0; double Q = 0; double H = 0; int lineTag = 0;

    if (argc - count <= 0)
    {
        opserr << "WARNING invalid aguments" << endln;
        opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
        return TCL_ERROR;
    }

    //Add a tag for location of origin;
    if (strcmp(argv[count], "-origin") == 0 || strcmp(argv[count], "origin") == 0) {
        count++;
        if (Tcl_GetDouble(interp, argv[count], &crd1) != TCL_OK) {
            opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;
        if (Tcl_GetDouble(interp, argv[count], &crd2) != TCL_OK) {
            opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;

        if (Tcl_GetDouble(interp, argv[count], &crd3) == TCL_OK) {
            //if the z loc is successfully recieved, count should be added with 1;
            count++;
        }
        else
        {
            //it it possible not to have a z loc for localised fire definiton
            opserr << "WARNING invalid z axis coordinate of fire origin" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            crd3 = 0.0;
        }
    }
    //end of fire origin, waiting for firePars;
    if (strcmp(argv[count], "-firePars") == 0 || strcmp(argv[count], "firePars") == 0) {
        count++;

        if (Tcl_GetDouble(interp, argv[count], &D) != TCL_OK) {
            opserr << "WARNING invalid diameter of the fire source" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;
        if (Tcl_GetDouble(interp, argv[count], &Q) != TCL_OK) {
            opserr << "WARNING invalid rate of heat release" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;
        if (Tcl_GetDouble(interp, argv[count], &H) != TCL_OK) {
            opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;
        //detect argument for linetag;
        if (argc - count > 0) {
            if (Tcl_GetInt(interp, argv[count], &lineTag) != TCL_OK) {
                opserr << "WARNING invalid central line tag " << endln;
                opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
                return TCL_ERROR;
            }
        }
        else {
            opserr << "Central line tag for localised fire " << argv[2] << "is set as default:3" << endln;
            lineTag = 3;
        }
    }
    else {
        opserr << "WARNING:: Defining Localised fire " << argv[2]
            << " expects tag:-firePars or firePars" << "\n";
    }

    theFireModel = new LocalizedFireEC1(FireModelTag, crd1, crd2, crd3, D, Q, H, lineTag);
    }
	//else ----------
	else if(strcmp(argv[1],"idealised") == 0||strcmp(argv[1],"Idealised") == 0){
		
		count++; //count should be updated
		double crd1=0.0; double crd2=0.0; double crd3=0.0; 
		double Q=0; double D1=0; double D2 =0 ; double K1 = 0; double K2 =0; int lineTag=0;
		int typeTag =1;
		if(argc-count<=0)
		{
			opserr << "WARNING invalid aguments" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
		}

		//Add a tag for location of origin;
		if(strcmp(argv[count],"-origin") == 0||strcmp(argv[count],"origin") == 0){
			count++;
			if (Tcl_GetDouble(interp, argv[count], &crd1) != TCL_OK) {
			opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &crd2) != TCL_OK) {
			opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			
			if (Tcl_GetDouble(interp, argv[count], &crd3) == TCL_OK) {
				//if the z loc is successfully recieved, count should be added with 1;
			count++;
			}
			else 
			{
				//it it possible not to have a z loc for localised fire definiton
				opserr << "WARNING invalid z axis coordinate of fire origin" << endln;
				opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	
				crd3=0.0;
			}
		}
		//end of fire origin, waiting for firePars;

		if(strcmp(argv[count],"-q") == 0||strcmp(argv[count],"q")==0||strcmp(argv[count],"-Q") == 0){
			count++;
			if (Tcl_GetDouble(interp, argv[count], &Q) != TCL_OK) {
			opserr << "WARNING invalid rate of heat release" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
		}
		if(strcmp(argv[count],"-Linear") == 0||strcmp(argv[count],"Linear")==0||strcmp(argv[count],"-linear") == 0){
			typeTag =1;
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D1) != TCL_OK) {
			opserr << "WARNING invalid range of pleatou" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D2) != TCL_OK) {
			opserr << "WARNING invalid length of decay" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
		}
		else if(strcmp(argv[count],"-quadratic") == 0||strcmp(argv[count],"Quadratic")==0||strcmp(argv[count],"-quadratic") == 0){
			typeTag =2;
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D1) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D2) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &K1) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
		}
		else if(strcmp(argv[count],"-exp") == 0||strcmp(argv[count],"exponential")==0||strcmp(argv[count],"-Exp") == 0){
			typeTag =3;
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D1) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &D2) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &K1) != TCL_OK) {
			opserr << "WARNING invalid distance between the fire source and the ceiling" << endln;
			opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
			count++;
		}
        else if (strcmp(argv[count], "-constant") == 0 || strcmp(argv[count], "uniform") == 0 || strcmp(argv[count], "-uniform") == 0) {
            count++;
            D1 = 0.0;
            D2 = 0.0;
            lineTag = 1;
            count++;

        }
			//detect argument for linetag;
		if(argc-count>0){
			if(strcmp(argv[count],"-centreLine") == 0||strcmp(argv[count],"CentreLine")==0||strcmp(argv[count],"-centreline") == 0){
				count++;
				if (Tcl_GetInt(interp, argv[count], &lineTag) != TCL_OK) {
				opserr << "WARNING invalid central line tag " << endln;
				opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;	    
				return TCL_ERROR;
				}
			}
		}
	
		if(typeTag ==1){
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3, Q,D1,  D2, lineTag);
		}
		else if(typeTag ==2){
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3,  Q, D1, D2,  K1,K2, lineTag);
		}
		else if(typeTag ==3){
			theFireModel = new Idealised_Local_Fire(FireModelTag, crd1, crd2, crd3,  Q, D1, D2,  K1, lineTag);
		}
		//Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3, double Q, double D1, double D2, double K1, double K2, int centerLineTag = 3);
  
		//Idealised_Local_Fire(int tag, double crd1, double crd2, double crd3,double Q, double D1, double D2, double factor, int centerLineTag = 3);

	}
	else{
		opserr<<"WARNING unknown fire model type"<< argv[1];
		opserr<<"- for fire model "<<FireModelTag;
		opserr<<" valid types: standard, hydroCarbon,parametric, localised.. \n";
		return TCL_ERROR;
   }
	
  if(theFireModel!=0){
		theTclHTModule->addFireModel(theFireModel);
#ifdef _DEBUG
		OPS_Stream* output = &opserr;
		theFireModel->Print(*output);
#endif
  }
  else
  {
		opserr<<"WARNING: TclHTModule fail to add FireModel: "<<argv[1]<<endln;
  }
  
  return TCL_OK;

}

//add setFIrePars
int
TclHeatTransferCommand_setFirePars(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    // checking the TclHTModule exists or not
    if (theTclHTModule == 0) {
        opserr << "WARNING current HeatTransfer Module has been destroyed - HTPattern\n";
        return TCL_ERROR;
    }
    // checking the HTDomain exists or not
    if (theHTDomain == 0) {
        opserr << "WARNING no active HeatTransfer Domain -  HTPattern\n";
        return TCL_ERROR;
    }

    //create a pointer to base class
    FireModel* theFireModel = 0;
    int FireModelTag = 0;
    int count = 1;
    int FireParTag = 0;
    double crd1 = 0;
    double crd2 = 0;
    double crd3 = 0;
    double q = 0;
    double d = 0;
    double Ts = 0;
    double addq = 0;

    if (strcmp(argv[count], "firemodel") == 0 || strcmp(argv[count], "-fireModel") == 0 || strcmp(argv[count], "fire") == 0)
    {

        count++;

        if (Tcl_GetInt(interp, argv[count], &FireModelTag) != TCL_OK) {
            opserr << "WARNING:: invalid masterID tag for coupling temperature: " << argv[1] << "\n";
            return TCL_ERROR;
        }
        else {
            count++;
        }
    }

    if (argc - count < 3) {
        opserr << "WARNING:: insufficient fire parameters for SetFirePars for fire model: " << FireModelTag << "\n";
        return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[count], &crd1) != TCL_OK) {
        opserr << "WARNING invalid x axis coordinate of fire origin" << endln;
        opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
        return TCL_ERROR;
    }
    count++;
    if (Tcl_GetDouble(interp, argv[count], &crd2) != TCL_OK) {
        opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
        opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
        return TCL_ERROR;
    }
    count++;

    if (Tcl_GetDouble(interp, argv[count], &crd3) == TCL_OK) {
        //if the z loc is successfully recieved, count should be added with 1;
        count++;
    }

    int remaining = argc - count;
    double dataInput = 0;
    double data[4];
    for (int i = 0; i < remaining; i++)
    {
        if (Tcl_GetDouble(interp, argv[count], &dataInput) != TCL_OK) {
            opserr << "WARNING invalid y axis coordinate of fire origin" << endln;
            opserr << " for HeatTransfer localised fire model: " << argv[2] << endln;
            return TCL_ERROR;
        }
        count++;
        data[i] = dataInput;
    }
    int numPars = 3 + remaining;
    Vector firepars(numPars);
    firepars(0) = crd1; firepars(1) = crd2; firepars(2) = crd3;
    if (numPars == 4) {
        firepars(3) = data[0];
    }
    else if (numPars == 5) {
        firepars(3) = data[0]; firepars(4) = data[1];
    }
    else if (numPars == 6) {
        firepars(3) = data[0]; firepars(4) = data[1]; firepars(5) = data[2] ;
    }
    else if (numPars == 7) {
        firepars(3) = data[0]; firepars(4) = data[1]; firepars(5) = data[2]; firepars(6) = data[3];
    }

    FireModel* thefire = 0;
    thefire = theTclHTModule->getFireModel(FireModelTag);
    if (thefire == 0) {
        opserr << "HeatTransferModule::SetFirePars failed to obtain the fire model" << endln;
        return -1;
    }
    double thecurrentTime = theHTDomain->getCurrentTime();

    thefire->setFirePars(thecurrentTime, firepars);


    return TCL_OK;


}

//Add HeatFluxBC
int
TclHeatTransferCommand_addMPTemperatureBC(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  
  //HeatFluxBC -eleset 1 -facetag 3 -type -ConvecAndRad -HTConstants 1;
  
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
    return TCL_ERROR;
  }
  
  if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
    return TCL_ERROR;
  }
  int masterIDtag,slaveIDtag;
  ID masterID;
  ID slaveID;
  
  int count=1;
  
  //if HTEntity tag is detected
  if(strcmp(argv[count],"-HTNodeSet") == 0||strcmp(argv[count],"-NodeSet") == 0||strcmp(argv[count],"nodeset") == 0)
  {
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &masterIDtag) != TCL_OK) {
      opserr << "WARNING:: invalid masterID tag for coupling temperature: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    if (Tcl_GetInt(interp, argv[count], &slaveIDtag) != TCL_OK) {
      opserr << "WARNING:: invalid masterID tag for coupling temperature: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
    HTNodeSet* theMasterNodeSet = theTclHTModule->getHTNodeSet(masterIDtag);
    if(theMasterNodeSet==0){
      opserr<< "WARNING:: TclHTModule failed to get the requested master NodeSet for coupling temperature: " <<argv[1] <<"\n";
      return TCL_ERROR;
    }
    
    HTNodeSet* theSlaveNodeSet = theTclHTModule->getHTNodeSet(slaveIDtag);
    if(theSlaveNodeSet==0){
      opserr<< "WARNING:: TclHTModule failed to get the requested slave NodeSet for coupling temperature: " <<argv[1] <<"\n";
      return TCL_ERROR;
    }
    
    masterID=theMasterNodeSet->getNodeID();
    slaveID = theSlaveNodeSet->getNodeID();
    
  }
// END OF OBTAINING NODES FORE COUPLING
  
  //------APPLYING MP_TEMPERATUREBC
  
  if (masterID.Size()==0||slaveID.Size()==0){
    opserr<<"Waring:No node will be coupled for tmperature"<<endln;
  }
  else if (masterID.Size()!=slaveID.Size())
  {
    opserr<<"WARNING::The size of master node set "<<masterIDtag
    <<" doesn't match the slave node set "<<slaveIDtag<<endln;
    return TCL_ERROR;
  }
  else
  {
    for (int i =0; i<masterID.Size();i++) {
        //Heat Transfer has only one DOF
        MP_TemperatureBC* MP_TempBC = new MP_TemperatureBC(masterID(i),slaveID(i));
        if(MP_TempBC==0)
          opserr << "WARNING: ran out of memory for coupling temperature"<<endln;
       
		if(theHTDomain->addMP_TemperatureBC(MP_TempBC) ==false)
        {
          opserr<< "Warning: could not add MP_TemperatureBC to domain" <<endln;
          delete MP_TempBC;
        }
    }
  }
return TCL_OK;
//end of adding mp temperature BC;
}
//



//Add HeatFluxBC
int
TclHeatTransferCommand_addHeatFluxBC(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  
  //HeatFluxBC -eleset 1 -facetag 3 -type -ConvecAndRad -HTConstants 1;
  
  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
    return TCL_ERROR;
  }
  
  //if (theTclHTBoudary == 0) {
   // opserr << "WARNING no active HeatTransfer Simple_Boundary - addHeatFluxBC\n";
    //return TCL_ERROR;
  //}
  
  int HTEleSetID = 0;
  int HTEntityTag = 0;
  Simple_Entity* theHTEntity = 0;
  HTEleSet* theHTEleSet = 0;
  int FaceID =0;
  ID EntityFaceID =0;
  ID ElesRange = 0;
  int HeatFluxTypeTag =0;
  int HTconstantsID = 0;
  int FireType=0;
  
  
  int count=1;
  
  //if HTEntity tag is detected
  if(strcmp(argv[count],"-HTEntity") == 0||strcmp(argv[count],"-Entity") == 0||strcmp(argv[count],"HTEntity") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTEntityTag) != TCL_OK) {
      opserr << "WARNING:: invalid HT Entity tag for defining HTEleSet: " << argv[1] << "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
  
	theHTEntity = theTclHTModule->getHTEntity(HTEntityTag);
    if(theHTEntity==0){
      opserr<< "WARNING:: HeatFluxBC failed to get the HTEntity: " <<argv[1] <<"\n";
      return TCL_ERROR;
    }

    //if face tag is detected
    if(strcmp(argv[count],"-face") == 0||strcmp(argv[count],"-Face") == 0||strcmp(argv[count],"-faceTag") == 0){
      
      count++;
      
      //-----for geting uncertain number of double data.
      int ArgStart = count;
      int ArgEnd = 0;
      int data;
      while (count < argc && ArgEnd == 0) {
       if (Tcl_GetInt(interp, argv[count], &data) != TCL_OK)
       ArgEnd = count;
       else
       count++;
      }
      //~ detecting the remianing number of input
	    EntityFaceID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          EntityFaceID(i-ArgStart) = data;
        }
      }
      // for geting uncertain number of doubel values
    }
  }else if(strcmp(argv[count],"-HTEleSet") == 0||strcmp(argv[count],"-EleSet") == 0||strcmp(argv[count],"eleset") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTEleSetID) != TCL_OK) {
      opserr << "WARNING:: invalid HTEleSet tag for defining HeatFluxBC: "<< "\n";
      return TCL_ERROR;
    }else{
      count++;
    }
    
	theHTEleSet = theTclHTModule->getHTEleSet(HTEleSetID);

	if(theHTEleSet==0){
		opserr << "WARNING:: invalid HTEleSet tag for defining HeatFluxBC: "<< "\n";
		return TCL_ERROR;
	}
    //if face tag is detected
    if(strcmp(argv[count],"-face") == 0||strcmp(argv[count],"-Face") == 0||strcmp(argv[count],"-faceTag") == 0){
      
      count++;
      
      if (Tcl_GetInt(interp, argv[count], &FaceID) != TCL_OK) {
        opserr << "WARNING:: invalid face tag for defining FaceID: " << argv[1] << "\n";
        return TCL_ERROR;
      }else{
        count++;
      }
    }
    
    if(FaceID==0){
      opserr<<"WARNING::no faceID specified for defining heat flux BC"<<endln;
      opserr<<" valid tags: -face, -Face, -faceTag "<<endln;
	  return TCL_ERROR;
    }
    
  }
  else{
    opserr<<"WARNING TclHTCommand HeatFluxBC wants  -HTEntity, HTEntity, -Entity or -HTElesET, -EleSet, eleset "<<endln;
	return TCL_ERROR;
  }
  //Simple_Boundary::GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType)
  
  //end of obtaining face or entity information
  
  //type is identified
  if(strcmp(argv[count],"-type") == 0||strcmp(argv[count],"-Type") == 0||strcmp(argv[count],"type") == 0){
    count++;
    
    if(strcmp(argv[count],"-convec") == 0||strcmp(argv[count],"-Convec") == 0||strcmp(argv[count],"Convec") == 0){
      HeatFluxTypeTag = 1;
    }
    else if(strcmp(argv[count],"-Radiation") == 0||strcmp(argv[count],"-Rad") == 0||strcmp(argv[count],"Rad") == 0){
      HeatFluxTypeTag = 2;
    }
    else if(strcmp(argv[count],"-Prescribed") == 0||strcmp(argv[count],"-Pres") == 0||strcmp(argv[count],"Pres") == 0){
      HeatFluxTypeTag = 3;
    }
    else if(strcmp(argv[count],"-ConvecAndRad") == 0||strcmp(argv[count],"-convecAndrad") == 0||strcmp(argv[count],"ConvecAndRad") == 0){
      HeatFluxTypeTag = 4;
    }
    else{
      opserr<<"WARNING unknown heat flux type"<< argv[1];
      opserr<<" valid types: -Convec, -Radiation,-ConvecAndRad,and -Prescribed \n";
      return TCL_ERROR;
    }
    count++;
  }
  
  if(HeatFluxTypeTag==0){
    opserr<<"WARNING::no HeatFlux type specified for defining heat flux BC"<<endln;
    opserr<<" valid tags: -type, -Type, type "<<endln;
  }
  
  //if Tag HTConstants is detected
  if(strcmp(argv[count],"-HTConstants") == 0||strcmp(argv[count],"-htConstants") == 0||strcmp(argv[count],"-Constants") == 0){
    
    count++;
    
    if (Tcl_GetInt(interp, argv[count], &HTconstantsID) != TCL_OK) {
      opserr << "WARNING:: invalid constants tag for defining Heatflux BC: " << argv[1] << "\n";
      return TCL_ERROR;
    }
    count++;
  }

  if (argc - count > 0)
  {
      if (strcmp(argv[count], "-par") == 0 || strcmp(argv[count], "-Par") == 0 || strcmp(argv[count], "-firePar") == 0) {

          count++;

          if (Tcl_GetInt(interp, argv[count], &FireType) != TCL_OK) {
              opserr << "WARNING:: invalid par tag for defining extral fire parameter: " << argv[1] << "\n";
              return TCL_ERROR;
          }
          count++;
      }
  }
  
  
  if(HTconstantsID==0){
    opserr<<"WARNING::no HTConstants found for defining heat flux BC"<<endln;
    opserr<<" valid tags: -HTConstants, -htConstants, -Constants "<<endln;
  }
  
  //Create boudary pattern pointed to theTclHTPattern;
  BoundaryPattern* thePattern = theTclHTPattern;
  if(theTclHTPattern==0){
	  opserr<<"WARNING:: TclHeatTransferModule failed to find Boundary pattern"<<endln;
  }
  //get the pattern tag
  int PatternTag = thePattern->getTag();
  
  
  //obtain heat flux Constants;
  HTConstants* theHTConstants = theTclHTModule->getHTConstants(HTconstantsID);
  if(theHTConstants==0){
	 opserr<<"WARNING::TclHeatTransferModule failed to find required HTConstants with tag: "<<HTconstantsID<<endln;
	 return TCL_ERROR;
  }
  Vector HeatFluxConstants = theHTConstants->getConstants();


  opserr << HeatFluxConstants << endln;
  
  if(HeatFluxConstants==0){
    opserr<<"WARNING TclHTModule:addHeatFluxBC failed to get HTConstants "<<HTconstantsID<<endln;
  }
  double ambqir = HeatFluxConstants(1) * HeatFluxConstants(1) * HeatFluxConstants(1) * HeatFluxConstants(1) * HeatFluxConstants(4);
  
  //------------------Applying heat flux BC--------------------------//
  //-----------------------------------------------------------------//
  if(strcmp(argv[1],"-HTEntity") == 0||strcmp(argv[1],"-Entity") == 0||strcmp(argv[1],"HTEntity") == 0){
  
    int EntityMeshTag = theHTEntity->getMeshTag();
    Simple_Mesh* theHTMesh = theTclHTModule->getHTMesh(EntityMeshTag);

    int NumFaces = EntityFaceID.Size();
    //loop over all the faces;
    
    for(int i=0 ; i<NumFaces; i++){
		//if face tag is detected
      ID ElesRange(0);
      int EntFaceID = EntityFaceID(i);
	  int FaceID;
      theHTMesh->SelectingElesbyFace(ElesRange,EntFaceID,FaceID );
	  
	  //detect the num of existing HeatfluxBCs;
     int ExistingHeatFluxBCs = thePattern->getNumHeatFluxBCs();

	  //check the number of heat flux BCs being applied
	int NumHeatFluxBCs = ElesRange.Size();
	if (ElesRange==0||NumHeatFluxBCs==0)
	{
		opserr<<"WARNING:: TclHeatTransferModule: No HeatFluxBc will be defined for Boundary pattern"<<endln;
	}
      
    
    if(HeatFluxTypeTag==1)
    {
      for(int i= 0;i<NumHeatFluxBCs; i++){
        Convection* Convec_BC = new Convection(ExistingHeatFluxBCs+i,ElesRange(i),FaceID,HeatFluxConstants(0),HeatFluxConstants(1));
        theHTDomain->addHeatFluxBC(Convec_BC,PatternTag);
      }
    }
    else if(HeatFluxTypeTag==2)
    {
      
      for(int i= 0;i<NumHeatFluxBCs; i++){  
        Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs+i,ElesRange(i), FaceID, HeatFluxConstants(2), HeatFluxConstants(4), HeatFluxConstants(3), ambqir);
        theHTDomain->addHeatFluxBC(Radiat_BC,PatternTag);
        
      }
    }
    else if(HeatFluxTypeTag ==3)
    {
      int NumNPF = (theHTDomain->getElement(ElesRange(0)))->getNumNodesperFace();
      for(int i= 0;i<NumHeatFluxBCs; i++){
        PrescribedSurfFlux* PreSurfFlux_BC = new PrescribedSurfFlux(ExistingHeatFluxBCs+i,ElesRange(i),FaceID,NumNPF,FireType);
        theHTDomain->addHeatFluxBC(PreSurfFlux_BC,PatternTag);
      }
    }
    else if(HeatFluxTypeTag ==4)
    {
      for(int i= 0;i<NumHeatFluxBCs; i++){
        Convection* Convec_BC = new Convection(ExistingHeatFluxBCs+2*i,ElesRange(i),FaceID,HeatFluxConstants(0),HeatFluxConstants(1));
        theHTDomain->addHeatFluxBC(Convec_BC,PatternTag);
        
        Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs+2*i+1,ElesRange(i), FaceID, HeatFluxConstants(2), HeatFluxConstants(4), HeatFluxConstants(3), ambqir);
        theHTDomain->addHeatFluxBC(Radiat_BC,PatternTag);
      }
    }
#ifdef _DEBUG
	opserr<<"TclHeatTransferModule:: "<<NumHeatFluxBCs <<" eles are applied with HeatFlux(Type: "
		<< HeatFluxTypeTag <<" ) at the face "<< EntFaceID << " ,Entity :"<<HTEntityTag<< endln;
#endif
    //end of heat flux type
  }
    //end of the loop over all faces

 }
   //end of if elements defined by entity face
  else if(strcmp(argv[1],"-HTEleSet") == 0||strcmp(argv[1],"-EleSet") == 0||strcmp(argv[1],"eleset") == 0){
    
	ElesRange = theHTEleSet->getEleID();
    BoundaryPattern* thePattern = theTclHTPattern;
    int PatternTag = thePattern->getTag();
    int NumHeatFluxBCs = ElesRange.Size();
	if (ElesRange==0||NumHeatFluxBCs==0)
	{
		opserr<<"WARNING:: TclHeatTransferModule: No HeatFluxBc will be defined for Boundary pattern"<<endln;
	}
    
    int ExistingHeatFluxBCs = thePattern->getNumHeatFluxBCs();
    
    Vector HeatFluxConstants = (theTclHTModule->getHTConstants(HTconstantsID))->getConstants();
    
    
	  if(HeatFluxConstants==0){
    opserr<<"WARNING TclHTModule:addHeatFluxBC failed to get HTConstants "<<HTconstantsID<<endln;
	  }
    
    if(HeatFluxTypeTag==1)
    {
      for(int i= 0;i<NumHeatFluxBCs; i++){
        Convection* Convec_BC = new Convection(ExistingHeatFluxBCs+i,ElesRange(i),FaceID,HeatFluxConstants(0),HeatFluxConstants(1));
        theHTDomain->addHeatFluxBC(Convec_BC,PatternTag);
      }
    }
    else if(HeatFluxTypeTag==2)
    {
      
      for(int i= 0;i<NumHeatFluxBCs; i++){
        
        Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs+i,ElesRange(i), FaceID, HeatFluxConstants(2), HeatFluxConstants(4), HeatFluxConstants(3), ambqir);
        theHTDomain->addHeatFluxBC(Radiat_BC,PatternTag);
        
      }
    }
    else if(HeatFluxTypeTag ==3)
    {
      int NumNPF = (theHTDomain->getElement(ElesRange(0)))->getNumNodesperFace();
      for(int i= 0;i<NumHeatFluxBCs; i++){
        PrescribedSurfFlux* PreSurfFlux_BC = new PrescribedSurfFlux(ExistingHeatFluxBCs+i,ElesRange(i),FaceID,NumNPF,FireType);
        theHTDomain->addHeatFluxBC(PreSurfFlux_BC,PatternTag);
      }
    }
    else if(HeatFluxTypeTag ==4)
    {
      for(int i= 0;i<NumHeatFluxBCs; i++){
        Convection* Convec_BC = new Convection(ExistingHeatFluxBCs+2*i,ElesRange(i),FaceID,HeatFluxConstants(0),HeatFluxConstants(1));
        theHTDomain->addHeatFluxBC(Convec_BC,PatternTag);
        
        Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs+2*i+1,ElesRange(i),FaceID,HeatFluxConstants(2),HeatFluxConstants(4),HeatFluxConstants(3), ambqir);
        theHTDomain->addHeatFluxBC(Radiat_BC,PatternTag);
      }
    } 
	#ifdef _DEBUG
	opserr<<"TclHeatTransferModule:: "<<NumHeatFluxBCs <<" eles are applied with HeatFlux(Type: "
		<< HeatFluxTypeTag <<" ) for the HTeleSet :"<<HTEleSetID<< endln;
	#endif
  }
  return TCL_OK;
  //Simple_Boundary::GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType)
  
}

//Add HTWipe
int TclHeatTransferCommand_HTReset(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	theHTDomain->clearAll();

	theHTDomain = 0;
	theTclHTModule =0;
	theTclHTPattern=0;
	return TCL_OK;
}


//Add HT_transient analysis
int TclHeatTransferCommand_HTTest(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
        opserr << "WARNING need to specify an convergence test (Static, Transient)\n";
        return TCL_ERROR;
    }

    int count = 1;

            if (strcmp(argv[count], "Residual") == 0 || strcmp(argv[count], "residual") == 0) {
                count++;
                double testTol;
                if (Tcl_GetDouble(interp, argv[count], &testTol) != TCL_OK) {
                    opserr << "WARNING test object tolerance must be a double.\n";
                    return TCL_ERROR;
                }
                count++;
                int maxIterations;
                if (Tcl_GetInt(interp, argv[count], &maxIterations) != TCL_OK) {
                    opserr << "WARNING test object maximum allowable iteration must be an integer.\n";
                    return TCL_ERROR;
                }
                count++;
                int analysisFlag;
                if (Tcl_GetInt(interp, argv[count], &analysisFlag) != TCL_OK) {
                    opserr << "WARNING test object flag must be an integer (0,1,2,3).\n";
                    return TCL_ERROR;
                }
                count++;
                theTest = new CTestNormResidual(testTol, maxIterations, analysisFlag);
                opserr << "Using NormResidual test with tolerance = " << testTol << ", max iterations = " << maxIterations << " and analysis flag = " << analysisFlag << ".\n";
            }
            else if (strcmp(argv[count], "TempIncr") == 0 || strcmp(argv[count], "temperature") == 0) {
                count++;
                double testTol=0.0;
                if (Tcl_GetDouble(interp, argv[count], &testTol) != TCL_OK) {
                    opserr << "WARNING test object tolerance must be a double.\n";
                    return TCL_ERROR;
                }
                count++;
                int maxIterations =0;
                if (Tcl_GetInt(interp, argv[count], &maxIterations) != TCL_OK) {
                    opserr << "WARNING test object maximum allowable iteration must be an integer.\n";
                    return TCL_ERROR;
                }
                count++;
                int analysisFlag=0;
                if (Tcl_GetInt(interp, argv[count], &analysisFlag) != TCL_OK) {
                    opserr << "WARNING test object flag must be an integer (0,1,2,3).\n";
                    return TCL_ERROR;
                }
                

                theTest = new CTestNormTempIncr(testTol, maxIterations, analysisFlag);
                opserr << "Using NormTempIncr test with tolerance = " << testTol << ", max iterations = " << maxIterations << " and analysis flag = " << analysisFlag << ".\n";

            }


 
  
    return TCL_OK;
}

//Add HT_transient analysis
int TclHeatTransferCommand_HTAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	 // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING need to specify an analysis type (Static, Transient)\n";
	return TCL_ERROR;
    }   
	 if (theHTAnalysis != 0) {
	delete theHTAnalysis;
	theHTAnalysis = 0;
	 }
     int count = 1;
     if (strcmp(argv[count],"HeatTransfer") == 0) {
         count++;
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new HT_AnalysisModel();
    if (argc == count) {
        //using default algorithm
        //theTest = 0;
        theAlgorithm = 0;
    }
    else {

        if (strcmp(argv[count], "Newton") == 0 || strcmp(argv[count], "newton") == 0) {
            count++;
            theAlgorithm = new NewtonMethod(*theTest);
            opserr << "Using the NewtonMethod algorithm.\n";
        }
        else if (strcmp(argv[count], "ModifiedNewton") == 0 || strcmp(argv[count], "modifiedNewton") == 0 || strcmp(argv[count], "modifiednewton") == 0) {
            count++;
            theAlgorithm = new ModifiedNewtonMethod(*theTest);
            opserr << "Using the ModifiedNewtonMethod algorithm.\n";
        }

    }
   

    if (theTest == 0) {
        opserr << "WARNING analysis Transient - no convergence test yet specified, \n";
        opserr << " CTestNormTempIncr default will be used\n";
        // theTest = new CTestNormTempIncr(1e-3, 500,1);
        //theTest = new CTestNormResidual(1e-1, 2000, 1);
        theTest = new CTestNormTempIncr(1e-3, 2000, 0);
    }

	if (theAlgorithm == 0) {
	    opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	    opserr << " NewtonMethod default will be used\n";	    
#ifdef _DEBUG
        theAlgorithm = new NewtonMethod(*theTest);
#else
        theAlgorithm = new NewtonMethod(*theTest);
#endif
	   
	}
	if (theHandler == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
	    opserr << " yet specified, PenaltyBC_Handler default will be used\n";
	    theHandler = new PenaltyBC_Handler(1.0e10,1.0e10);    
	}
	if (theNumberer == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
	    opserr << " RCM default will be used\n";
	    RCM *theRCM = new RCM(false);	
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
	    BandGenLinLapackSolver *theSolver;
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
	return TCL_OK;
}

//analyze HT_transient model
int TclHeatTransferCommand_HTAnalyze(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int result = 0;

	if (theHTAnalysis != 0) {
    
	 if (argc < 3) {
      opserr << "WARNING heat transfer analysis: analysis numIncr? deltaT?\n";
      return TCL_ERROR;
     }
   
	 int numIncr;
    if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)	
      return TCL_ERROR;
    
	double dT;
	if (Tcl_GetDouble(interp, argv[2], &dT) != TCL_OK)	
      return TCL_ERROR;

	double time = 0;
    result = theHTAnalysis->analyze(numIncr, dT,time);

	if (result < 0) 
    opserr << "OpenSees > analyze failed, returned: " << result << " error flag\n";
	
	}
	//Output for the final result of heat transfer analysis;
	opserr<<result<<endln;
	return TCL_OK;
}

int TclHeatTransferCommand_HTRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

	Vector* theRecVec=0;
	  if (theTclHTModule == 0) {
    opserr << "WARNING current HeatTransfer Module has been destroyed - addHeatFluxBC\n";
    return TCL_ERROR;
	}
  
	if (theHTDomain == 0) {
    opserr << "WARNING no active HeatTransfer Domain - addHeatFluxBC\n";
    return TCL_ERROR;
   }
    OPS_Stream *theOutputStream = 0;
	TCL_Char *fileName = 0;
	ID* theNodes =0;
	HTRecorder* theHTRecorder =0;
	int count=1;
	double xloc = 0;
  
  //if file tag is detected
  if(strcmp(argv[count],"-file") == 0||strcmp(argv[count],"file") == 0||strcmp(argv[count],"-File") == 0)
	{
    
    fileName = argv[count+1];
	const char *pwd = getInterpPWD(interp);
	simulationInfo.addOutputFile(fileName, pwd);
	theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 0 );  //The last argument is to determine CSV for using space or 
	count += 2;
	}
  
  ///////-----------------------HTRecorderToStructural member-----------------------------------//
  if(strcmp(argv[count],"-xloc") == 0||strcmp(argv[count],"xLoc") == 0||strcmp(argv[count],"-xLoc") == 0)
  {
	 count++;
	 if (Tcl_GetDouble (interp, argv[count], &xloc) != TCL_OK) {
			opserr << "WARNING invalid xloc" << endln;
			opserr << " for HeatTransfer recorder: " << endln;	    
			return TCL_ERROR;
			}
	 count++;
    //ending of xloc;
  
	if(strcmp(argv[count],"-yloc") == 0||strcmp(argv[count],"yLoc") == 0||strcmp(argv[count],"-yLoc") == 0)
    {
	 count++;
	 Vector yloc(argc-count);
    //-----for geting uncertain number of double data.
     int ArgStart = count;
     int ArgEnd = argc;
     double data;
     /*while (count < argc && ArgEnd == 0) {
	 if (Tcl_GetDouble(interp, argv[count], &data) != TCL_OK)
		ArgEnd = count;
		else
		 count++;
	}*/
      if (ArgStart != ArgEnd) {
		for (int i=ArgStart; i<ArgEnd; i++) {
		Tcl_GetDouble(interp, argv[i], &data);
		yloc(i-ArgStart) = data;
		}
      }
	  //ending of yloc;
    
	  int NumRows = yloc.Size();
      theRecVec =  new Vector(NumRows+1);
	  theRecVec->Zero();

      for(int i=0;i<=NumRows;i++) {
        if (i==0) {
          (*theRecVec)(i)= xloc;
        } else {
          (*theRecVec)(i)= yloc(i-1);
        }
      }
   }
#ifdef _DEBUG
   opserr<< "TclHeatTransferModule::HTRecorder, theRecMatrix "<<*theRecVec<<endln;
#endif
	HTReorderTag++; 
   theHTRecorder = new HTRecorderToStru(HTReorderTag,*theRecVec,*theHTDomain,*theOutputStream);
   }

  ///////-----------------------HTNodeRecorder-----------------------------------//
  else if(strcmp(argv[count],"-nodeSet") == 0||strcmp(argv[count],"NodeSet") == 0||strcmp(argv[count],"-NodeSet") == 0)
  {
    count++;
	int NodeSetID =0;
	if (Tcl_GetInt(interp, argv[count], &NodeSetID) != TCL_OK) {
      opserr << "WARNING:: invalid nodeSet tag for defining HTNodeRecorder : " << "\n";
      return TCL_ERROR;
    }
	
	HTNodeSet* theRecNodeSet = theTclHTModule->getHTNodeSet(NodeSetID);
    if(theRecNodeSet==0){
      opserr<< "WARNING:: TclHTModule failed to get the requested NodeSet for HTNodeRecorder: " <<"\n";
      return TCL_ERROR;
    }
    ID RecNodeID = 0;
   RecNodeID = theRecNodeSet->getNodeID();
  
#ifdef _DEBUG
   opserr<< "TclHeatTransferModule::HTRecorder, theRecNodeID "<<RecNodeID<<endln;
#endif
   theHTRecorder = new HTNodeRecorder(HTReorderTag++,&RecNodeID,*theHTDomain,*theOutputStream);

  }
  //end of nodeset
  else if (strcmp(argv[count], "-node") == 0 || strcmp(argv[count], "Node") == 0 || strcmp(argv[count], "-Node") == 0)
  {
      count++;

      int numnodes = argc - count;
      ID RecNodeID = ID(numnodes);
      int NodeTag = 0;
      for (int i = 0; i < numnodes; i++) {

          if (Tcl_GetInt(interp, argv[count], &NodeTag) != TCL_OK) {
              opserr << "WARNING:: invalid nodeSet tag for defining HTNodeRecorder : " << "\n";
              return TCL_ERROR;
          }
          RecNodeID(i) = NodeTag;
          count++;
      }
     

#ifdef _DEBUG
      opserr << "TclHeatTransferModule::HTRecorder, theRecNodeID " << RecNodeID << endln;
#endif
      theHTRecorder = new HTNodeRecorder(HTReorderTag++, &RecNodeID, *theHTDomain, *theOutputStream);

  }
 else if(strcmp(argv[count],"-HTEntity") == 0||strcmp(argv[count],"htEntity") == 0||strcmp(argv[count],"-htentity") == 0)
  {
    count++;
	int EntityID =0; int DimTag = 0;
	if (Tcl_GetInt(interp, argv[count], &EntityID) != TCL_OK) {
      opserr << "WARNING:: invalid Entity tag for defining HTNodeRecorder : " << "\n";
      return TCL_ERROR;
    }

	count++;
	if(argc>count){
		if(strcmp(argv[count],"-dim") == 0||strcmp(argv[count],"dim") == 0||strcmp(argv[count],"Dim") == 0){
			count++;
		}

		if (Tcl_GetInt(interp, argv[count], &DimTag) != TCL_OK) {
			opserr << "WARNING:: invalid Entity tag for defining HTNodeRecorder : " << "\n";
			return TCL_ERROR;
		}
	}

	Simple_Entity* theRecEntity = theTclHTModule->getHTEntity(EntityID);
	int theMshTag = theRecEntity->getMeshTag();
	Simple_Mesh* theHTMesh = theTclHTModule->getHTMesh(theMshTag);

    if(theHTMesh==0){
      opserr<< "WARNING:: TclHTModule failed to get the requested Entity for HTNodeRecorder: " <<"\n";
      return TCL_ERROR;
    }

    ID RecNodesID = 0;
	theHTMesh->GetNodesForRecorder(RecNodesID, DimTag);
  
#ifdef _DEBUG
   opserr<< "TclHeatTransferModule::HTRecorder, theRecNodeID "<<RecNodesID<<endln;
#endif
   theHTRecorder = new HTNodeRecorder(HTReorderTag++,&RecNodesID,*theHTDomain,*theOutputStream);

  }
  ///////-----------------------HTEleRecorder-----------------------------------//
 else if (strcmp(argv[count], "-EleSet") == 0 || strcmp(argv[count], "eleset") == 0 || strcmp(argv[count], "EleSet") == 0)
  {
     count++;
     int EleSetID = 0;
    if (Tcl_GetInt(interp, argv[count], &EleSetID) != TCL_OK) {
      opserr << "WARNING:: invalid nodeSet tag for defining HTNodeRecorder : " << "\n";
      return TCL_ERROR;
    }

    HTEleSet* theRecNodeSet = theTclHTModule->getHTEleSet(EleSetID);
    if (theRecNodeSet == 0) {
      opserr << "WARNING:: TclHTModule failed to get the requested NodeSet for HTNodeRecorder: " << "\n";
      return TCL_ERROR;
    }
    ID RecEleID = 0;
     RecEleID = theRecNodeSet->getEleID();

    #ifdef _DEBUG
  opserr << "TclHeatTransferModule::HTRecorder, RecEleID " << RecEleID << endln;
    #endif
//    HTElementRecorder(const ID* eleID,const char** argv, int argc, bool echoTime, HeatTransferDomain& theDomain,  OPS_Stream& theOutputHandler)
    theHTRecorder = new HTElementRecorder(2, &RecEleID, argv,argc, true, *theHTDomain, *theOutputStream);



  }
	
	//HTRecorderToStru::HTRecorderToStru(int tag, const Matrix& theCrds, HeatTransferDomain &theDom, OPS_Stream &theOutputHandler, double tolerance)
// for geting uncertain number of doubel values 

  if(theHTRecorder!=0){
   theHTDomain->addRecorder(*theHTRecorder);
	  }
  else{
	  opserr<<"WARNING::HeatTransfer Recoder is not defined"<<endln;
	  }
// HTNodeRecorder(int tag, const ID* theNodes, HeatTransferDomain& theDomain,OPS_Stream &theOutputHandle); 

return TCL_OK;
}

int 
TclHeatTransferCommand_PrintNodes(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv) {

    if (theTclHTModule == 0) {
        opserr << "WARNING current HeatTransfer Module has been destroyed\n";
        return TCL_ERROR;
    }

    if (theHTDomain == 0) {
        opserr << "WARNING no active HeatTransfer Domain\n";
        return TCL_ERROR;
    }

    if (argc > 1) {
        opserr << "Too many arguments for HTPrintNodes, which takes no arguments at all!" << endln;
        return TCL_ERROR;
    }

    int NumNodes = 0;
    int NodeTag = 0;
    double Nodalx; 
    double Nodaly;
    double Nodalz;


    NumNodes = theHTDomain->getNumNodes();
    opserr << "There are " << NumNodes << " in the current domain." << endln;
    opserr << "Node Tag, x, y, z " << endln;
    for (int i = 0; i < NumNodes; i++) {
        NodeTag = i + 1;
        Nodalx = (theHTDomain->getNode(NodeTag)->getCrds())(0);
        Nodaly = (theHTDomain->getNode(NodeTag)->getCrds())(1);
        Nodalz = (theHTDomain->getNode(NodeTag)->getCrds())(2);
        if (abs(Nodalx) < 1e-6 || abs(Nodalx) > 1e6)
            Nodalx = 0;
        if (abs(Nodaly) < 1e-6 || abs(Nodaly) > 1e6)
            Nodaly = 0;
        if (abs(Nodalz) < 1e-6 || abs(Nodalz) > 1e6)
            Nodalz = 0;
        opserr << NodeTag << ", "<< Nodalx <<", "<< Nodaly<<", "<< Nodalz << endln;
    }

    return 0;
}
int
TclHeatTransferCommand_getHTTime(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv) {

    if (theTclHTModule == 0) {
        opserr << "WARNING current HeatTransfer Module has been destroyed\n";
        return TCL_ERROR;
    }

    if (theHTDomain == 0) {
        opserr << "WARNING no active HeatTransfer Domain\n";
        return TCL_ERROR;
    }

    if (argc > 1) {
        opserr << "Too many arguments for getHTTime, which takes no arguments at all!" << endln;
        return TCL_ERROR;
    }

    double time = theHTDomain->getCurrentTime();

    // get the display format
    char format[80];
    if (argc == 1) {
        //      strcpy(format,"%f");
        sprintf(format, "%f", time);
    }
    else if (argc == 2) {
        //      strcpy(format,argv[1]);
        sprintf(format, argv[1], time);
    }

    // now we copy the value to the tcl string that is returned
    //  sprintf(interp->result,format,time);
    Tcl_SetResult(interp, format, TCL_VOLATILE);
    return TCL_OK;
}