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
**                                                                    **                      **
** ****************************************************************** */
                                                                        
/**********************************************************************
** This project aims to provide an integrated computational 
** environment for modelling structural behaviours under fire action. **
** Developed by:                  `                                   **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $


//
// Written by Liming Jiang (Liming.Jiang@ed.ac.uk)
//

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

#include <StandardStream.h>
#include <DataFileStream.h>
#include <SimulationInformation.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>
	
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
// includes for the domain classes
#include <Simple_Brick.h>
#include <Simple_Block.h>
#include <Simple_Isection.h>
#include <Simple_Isection3D.h>
#include <Simple_Mesh.h>
#include <Simple_Boundary.h>
#include <Vector.h>
#include <Matrix.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <QuadFour.h>
#include <BrickEight.h>
#include <TemperatureBC.h>
#include <MP_TemperatureBC.h>
#include <Convection.h>

#include <CarbonSteelEC3.h>
#include <ConcreteEC2.h>
#include <BoundaryPattern.h>
#include <FireImposedPattern.h>

#include <SimpleMaterial.h>

// includes for the analysis classes
#include <HT_TransientAnalysis.h>
#include <BackwardDifference.h>
#include <HT_AnalysisModel.h>  
//#include <LinearAlgorithm.h>
#include <NewtonMethod.h>
//#include <ModifiedNewtonMethod.h>
#include <CTestNormTempIncr.h>
#include <PenaltyBC_Handler.h>
#include <RCM.h>
#include <HT_DOF_Numberer.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <ParametricFireEC1.h>
#include <LocalizedFireEC1.h>
#include <NorminalFireEC1.h>
#include <HTNodeRecorder.h>
#include <HTRecorderToStru.h>
#include <TimeSeries.h>

#include <TclSIFBuilder.h>
#include <SIFBuilderDomain.h>
#include <SIFJointIter.h>
#include <SIFJoint.h>
#include <SIFMaterial.h>
#include <SIFSection.h>

#include <SIFMember.h>
#include <SIFXBeam.h>
#include <SIFYBeam.h>
#include <SIFXBeamSec.h>
//#include <SIFYBeamSec.h>
#include <SIFColumn.h>
#include <SIFXWall.h>
#include <SIFYWall.h>
#include <SIFSlab.h>

#include <SIFJointIter.h>
#include <SIFXBeamIter.h>
#include <SIFYBeamIter.h>
#include <SIFColumnIter.h>
#include <SIFSlabIter.h>
#include <SIFSecXBeamIter.h>
#include <SIFCompartment.h>
#include <SIFfireAction.h>
#include <SIFHTforMember.h>

#include <NodeRecorder.h>

extern SimulationInformation simulationInfo;
extern const char * getInterpPWD(Tcl_Interp *interp);  // commands.cpp	
static SIFBuilderDomain* theSIFDomain = 0;       //static pointer to SIFBuilderDomain
static Domain *theDomain = 0;
//static TclSIFBuilder* theTclSIFBuilder  = 0;   //static pointer to
static Vector* XBayVec = 0;                      //static vetor pointer for xbay
static Vector* YBayVec = 0;                      //static vetor pointer for ybay
static Vector* StoreyVec = 0;                    //static vetor pointer for storey
static int SIFModelStatus = 0;                   //static index to detect whether the sif model beig created or not
static int BuilderType =0;						 //static index to define the builder type (structutal model)
static int SIFdisplayindex = 0;                  //static index to define the display need for SIFBuilder;
static Vector* BuilderTypeVec = new Vector(3);   //static vector for indentifying the builder type
static int SecXBeamTag = 0;
static int SecYBeamTag = 0;
static int theNewJtTag = 10000;                  //Additional joints are tagged as 10000+;

//Declaration of Tcl commands for SIFBuilder
int TclSIFBuilderCommand_addMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addSection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_assignSection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addXBay(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addYBay(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addSkew_Angle(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addStorey(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int TclSIFBuilderCommand_addSecBeam(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_delMember(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int TclSIFBuilderCommand_addFireAction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_applyAndanalyze(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_addFirePars(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_SetBC(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
 int TclSIFBuilderCommand_AddLoad(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_AddPartDamage(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


int TclSIFBuilderCommand_AddSIFRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int TclSIFBuilderCommand_BuildSIFModel();
int TclSIFBuilderCommand_BuildModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
//int TclSIFBuilderCommand_meshSIFModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

TclSIFBuilder::TclSIFBuilder(Domain &Domain, Tcl_Interp* interp, int displayindex)
{
  // Set the interpreter (the destructor needs it to delete commands)
  theInterp = interp;
  
  opserr<<"               ----SIFBuilder is developed by UoE Group----       "<<endln;


  if(theSIFDomain==0) {
	  theSIFDomain= new SIFBuilderDomain();
  }
  
  if(SIFModelStatus!=0)
	  SIFModelStatus=0;

  SIFdisplayindex = displayindex;

  theDomain= &Domain;
  theSIFDomain->SetStructureDomain(&Domain);

  Tcl_CreateCommand(interp, "AddMaterial", (Tcl_CmdProc* )TclSIFBuilderCommand_addMaterial,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddSection", (Tcl_CmdProc* )TclSIFBuilderCommand_addSection,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AssignSection", (Tcl_CmdProc* )TclSIFBuilderCommand_assignSection,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SIFXBay", (Tcl_CmdProc* )TclSIFBuilderCommand_addXBay,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddSecBeam", (Tcl_CmdProc* )TclSIFBuilderCommand_addSecBeam,(ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "DelMember", (Tcl_CmdProc*)TclSIFBuilderCommand_delMember, (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "SIFZBay", (Tcl_CmdProc* )TclSIFBuilderCommand_addYBay,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SIFStorey", (Tcl_CmdProc* )TclSIFBuilderCommand_addStorey,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddSkew_Angle", (Tcl_CmdProc* )TclSIFBuilderCommand_addSkew_Angle,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddFire", (Tcl_CmdProc* )TclSIFBuilderCommand_addFireAction,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddFirePars", (Tcl_CmdProc* )TclSIFBuilderCommand_addFirePars,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SIFAnalyze", (Tcl_CmdProc* )TclSIFBuilderCommand_applyAndanalyze,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "BuildModel", (Tcl_CmdProc* )TclSIFBuilderCommand_BuildModel,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SetBC", (Tcl_CmdProc* )TclSIFBuilderCommand_SetBC,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddLoad", (Tcl_CmdProc* )TclSIFBuilderCommand_AddLoad,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "AddPartDamage", (Tcl_CmdProc* )TclSIFBuilderCommand_AddPartDamage,(ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "SIFRecorder", (Tcl_CmdProc* )TclSIFBuilderCommand_AddSIFRecorder,(ClientData)NULL, NULL);
 // Tcl_CreateCommand(interp, "MeshSIFModel", (Tcl_CmdProc* )TclSIFBuilderCommand_meshSIFModel,(ClientData)NULL, NULL);

  
}

TclSIFBuilder::~TclSIFBuilder()
{
	theSIFDomain->clearAll();
	if(XBayVec!=0){
		delete XBayVec;
		XBayVec =0;
	}
	if(YBayVec!=0){
		delete YBayVec;
		YBayVec =0;
	}
	if(StoreyVec!=0){
		delete StoreyVec;
		StoreyVec =0;
	}
	
	 if(SIFModelStatus!=0)
	  SIFModelStatus=0;
	
	Tcl_DeleteCommand(theInterp, "AddMaterial");
	Tcl_DeleteCommand(theInterp, "SIFXBay");
	Tcl_DeleteCommand(theInterp, "SIFZBay");
	Tcl_DeleteCommand(theInterp, "SIFStorey");
	Tcl_DeleteCommand(theInterp, "AddSecBeam");
	Tcl_DeleteCommand(theInterp, "AddFire");
	Tcl_DeleteCommand(theInterp, "SIFAnalyze");
	Tcl_DeleteCommand(theInterp, "AddFirePars");
	Tcl_DeleteCommand(theInterp, "AddSkew_Angle");
	Tcl_DeleteCommand(theInterp, "AddSection");
	Tcl_DeleteCommand(theInterp, "AssignSection");
	Tcl_DeleteCommand(theInterp, "SetBC");
	Tcl_DeleteCommand(theInterp, "AddLoad");
	Tcl_DeleteCommand(theInterp, "AddPartDamage");
	Tcl_DeleteCommand(theInterp, "BuildModel");
	delete theSIFDomain;
	theSIFDomain = 0;
	SecXBeamTag =0;
	SecYBeamTag=0;
	theNewJtTag = 10000;
	
}


//--------------------SIFModel generation------------------------------------------
int 
TclSIFBuilderCommand_BuildSIFModel( ){
     //Build model consists of three steps:
	 //   1) Build SIFModel, which stores the building information
	 //      --SIFXjoint(s)
	 //      --SIFMembers(s):SIFXBeam, SIFYBeam,SIFColumn
	 //   2) Build Structural Model
	     

	int NumXBays;			//number of bays in x direction
	int NumYBays;			//number of bays in y direction
	int NumStoreys;			//number of NumStoreys in z direction

	if(XBayVec!=0){
		NumXBays = (XBayVec->Size())-1;
		(*BuilderTypeVec)(0) = 1;
	}
	else{
		(*BuilderTypeVec)(0) = 0;
		NumXBays = 1;		// Here the real num of XBays is 0, just for conviniece of defining virtual compartment
	}

	if(YBayVec!=0){
		NumYBays = (YBayVec->Size())-1;
		(*BuilderTypeVec)(2) = 1;
	}
	else{
		(*BuilderTypeVec)(2) = 0;
		NumYBays = 1;     // Here the real num of YBays is 0, just for conviniece of defining virtual compartment
	}

	if(StoreyVec!=0){
		NumStoreys = (StoreyVec->Size())-1;
		(*BuilderTypeVec)(1) = 1;
	}
	else{
		(*BuilderTypeVec)(1) = 0;
		NumStoreys = 1;     // Here the real num of storey is 0, just for conviniece of defining virtual compartment
	}

	//3D full frame
	if ((*BuilderTypeVec)(0) == 1 && (*BuilderTypeVec)(1) == 1 && (*BuilderTypeVec)(2) == 1) {
		BuilderType = 1;		//1:L M N (X-Z Plan, Y storey-3D)
	}

    //Plane frame
	else if( (*BuilderTypeVec)(0) == 1 && (*BuilderTypeVec)(1) == 0 && (*BuilderTypeVec)(2) == 1){

			BuilderType = 11;		//11:M 0 N (X-Z grillage with out of plane displacements - 3D)
	}
	else if( (*BuilderTypeVec)(0) == 1 && (*BuilderTypeVec)(1) == 1 && (*BuilderTypeVec)(2) == 0){

			BuilderType = 12;		//12:M N 0 (X-Y frame with in and out of plane displacements - 3D)
	}
	else if( (*BuilderTypeVec)(0) == 0 && (*BuilderTypeVec)(1) == 1 && (*BuilderTypeVec)(2) == 1){

			BuilderType = 13;		//13:0 M N (Y-Z frame with in and out of plane displacements - 3D)
	}
	
	//Single member or continuous member
	else if( (*BuilderTypeVec)(0) == 0 && (*BuilderTypeVec)(1) == 1 && (*BuilderTypeVec)(2) == 0){

			BuilderType = 14;		//14:0 N 0 (single storey or multi-storey column - 3D)
	}
	else if( (*BuilderTypeVec)(0) == 1 && (*BuilderTypeVec)(1) == 0 && (*BuilderTypeVec)(2) == 0){

			BuilderType = 15;		//15:N 0 0 (single or continuous beam in the X direction - 3D)
	}
	else if( (*BuilderTypeVec)(0) == 0 && (*BuilderTypeVec)(1) == 0 && (*BuilderTypeVec)(2) == 1){

			BuilderType = 16;		//16:0 0 N (single or continuous beam in the Z direction - 3D)
	}

	ID theBuilderInfo = ID(4);
	theBuilderInfo(0)=BuilderType; theBuilderInfo(1)= NumXBays ;
	theBuilderInfo(2)= NumYBays ;theBuilderInfo(3)= NumStoreys;
	theSIFDomain->setSIFBuilderInfo(theBuilderInfo);

	//if(NumXBays<0

	double xBaycrd = 0.0;		//xBaycrd is the coordinate of joint in x direction;
	double yBaycrd = 0.0;		//yBaycrd is the coordinate of joint in y direction;
	double zBaycrd = 0.0;		//zBaycrd is the coordinate of joint in z direction;

	//-------------------------------------------for SIFCompartment generation---------------------------------------------;
   if(SIFdisplayindex==1|| SIFdisplayindex == 2)
	   opserr << "SIFCompartment: " << endln;
   //begin to generate SIFCompartment;
	for	(int b = 1; b <= NumYBays; b++)
  {
    for	(int c = 1; c <=NumStoreys; c++)
    {
      for (int a = 1; a <= NumXBays; a++)
      {
        //To set the crds for the SIFJoint, which are obtained from the vectors
		int CompTag = 1000*a+100*b +c;   
		//it can be CompTag++;
        //joint tag: abc -> a th XBay, b th YBay, c th Storey
        ID CompInfo =  ID(3);
		CompInfo(0)= a ; CompInfo(1)= b; CompInfo(2) = c; 
		Vector CompOrigin = Vector(3);     
        
		if(BuilderType == 11||BuilderType ==15||BuilderType ==16){
			StoreyVec = new Vector(1);
			(*StoreyVec)(0) = 0;
			CompOrigin(2) = (*StoreyVec)(c-1);
		}
		else{
			CompOrigin(2)= (*StoreyVec)(c-1);
		}

		if(BuilderType == 12||BuilderType ==14||BuilderType ==15){
			YBayVec = new Vector(1);
			(*YBayVec)(0) = 0;
			CompOrigin(1) = (*YBayVec)(b-1);
		}
		else{
			CompOrigin(1)= (*YBayVec)(b-1);
		}

		if(BuilderType == 13||BuilderType ==14||BuilderType ==16){
			XBayVec = new Vector(1);
			(*XBayVec)(0) = 0;
			CompOrigin(0) = (*XBayVec)(a-1);
		}
		else{
			CompOrigin(0)= (*XBayVec)(a-1);
		}

		SIFCompartment *theComp = new SIFCompartment(CompTag, CompInfo, CompOrigin);
		theSIFDomain->addSIFCompartment(theComp);
		if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
			opserr << CompTag<<"   ";
		//To get reference just add a * to get the class, which is then parsed as reference
	  }
	}
  }

	
	//-----------------------------For SIFJoint generation------------------------------------------------------------;
	if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
		opserr <<endln<< "SIFJoint: " << endln;

	//if a sub strcture, a, b start from 0
	int aStart,bStart,cStart,aEnd,bEnd,cEnd;

	if(BuilderType ==1||BuilderType ==2){					//frame strcture
		aStart = 0; bStart = 0; cStart = 0;
		aEnd = NumXBays; bEnd=NumYBays;cEnd=NumStoreys;
	}

	else if(BuilderType ==3){								//sub frame structure
		aEnd = NumXBays+2; bEnd=NumXBays+2;cEnd=NumStoreys+2;
	}

	else if(BuilderType ==11){								//grillage3D structure

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = NumXBays; bEnd=NumYBays;cEnd=0;
	}

	else if(BuilderType ==12){								// 3D plane frame structure in X-Y coordinates

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = NumXBays; bEnd=0;cEnd=NumStoreys;
	}

	else if(BuilderType ==13){								// 3D plane frame structure in Y-Z coordinates

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = 0; bEnd=NumYBays;cEnd=NumStoreys;
	}

	else if(BuilderType ==14){								// 3D single/multi storey column in Y direction

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = 0; bEnd=0;cEnd=NumStoreys;
	}

	else if(BuilderType ==15){								// 3D single/continuous beam in X direction

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = NumXBays; bEnd=0;cEnd=0;
	}

	else if(BuilderType ==16){								// 3D single/continuous beam in Z direction

		aStart = 0; bStart = 0; cStart = 0;
		aEnd = 0; bEnd=NumYBays;cEnd=0;
	}

	for (int b = bStart; b <= bEnd; b++)
	{
		for	(int c = cStart; c <= cEnd; c++)
		{
			for (int a = aStart; a <= aEnd; a++)
			{
			//To set the crds for the SIFJoint, which are obtained from the vectors
			xBaycrd= (*XBayVec)(a);
			yBaycrd= (*YBayVec)(b);
			zBaycrd= (*StoreyVec)(c);
        
			int JointTag;

			if(BuilderType == 11||BuilderType ==15||BuilderType ==16){							

				JointTag = 1000*(a+1)+100*(b+1) +c+1; 
			}
			else {
			JointTag = 1000*(a+1)+100*(b+1) +c; 
			}

			//"Tag" starts in upper case
			//joint tag: abc -> a th XBay, b th YBay, c th Storey
       
			SIFJoint *joint = new SIFJoint(JointTag, xBaycrd, zBaycrd, yBaycrd);
			theSIFDomain->addSIFJoint(joint);
			if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
				opserr <<JointTag<< "  " ;
			//To get reference just add a * to get the class, which is then parsed as reference
			}
		}
	}



if(BuilderType ==1||BuilderType ==11||BuilderType ==12||BuilderType ==13){

	// for Generating SIF slab;
	if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
		opserr << endln << "SIFSlab: " << endln;

	for	(int b = 1; b <= NumYBays; b++)
	{
		for	(int c = 1; c <=NumStoreys; c++)
		{
			for (int a = 1; a <= NumXBays; a++)
			{
				int Jt1Tag = 1000*a+100*b+c;
				int Jt2Tag = 1000*(a+1)+100*b+c;
				int Jt3Tag = 1000*(a+1)+100*(b+1)+c;
				int Jt4Tag = 1000*a+100*(b+1)+c;
				int SlabTag = 1000*a+100*b+c;
		

				ID MemberInfo =  ID(3);
				MemberInfo(0)= a ; MemberInfo(1)= b; MemberInfo(2) = c; SIFSlab* theSlab;
				//Xbeam locates in the (a+1)XBay, the b bay grid, the c storey;
				if(BuilderType ==12){
					theSlab = new SIFSlab(SlabTag,3, Jt1Tag, Jt2Tag, MemberInfo);
				}
				else if(BuilderType ==13){
					theSlab = new SIFSlab(SlabTag,3, Jt1Tag, Jt4Tag, MemberInfo);
				}
				else{
					theSlab = new SIFSlab(SlabTag,3, Jt1Tag, Jt2Tag, Jt3Tag, Jt4Tag, MemberInfo);
				}
				theSIFDomain->addSIFSlab(theSlab);
				//column->assignSection(theSIFDomain->getSIFSection(1));
				if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
					opserr << SlabTag << "  ";

				SIFCompartment* theComp=0;
				int theCompTag = 1000*a+100*b+c;
				theComp = theSIFDomain->getSIFCompartment(theCompTag);
				theComp->AddSlab(SlabTag);
				//To get reference just add a * to get the class, which is then parsed as reference
			}
		}
	}
}
  //end of SIFSlab generation;
  
	//for SIFxbeam generation

	if(BuilderType ==13||BuilderType ==14||BuilderType ==16){								
		
	}
	else{								

	double gamma=0.0;							//skew angle

	if(BuilderType ==12||BuilderType ==15){								
		bEnd=0;
	}
	else{								
		bEnd=NumYBays;
	}

	for	(int b = 0; b <= bEnd; b ++)
	{
		for	(int c = 1; c <=NumStoreys; c ++)
		{
			for (int a = 0; a<NumXBays; a ++)
			{
				int Jt1Tag = 1000*(a+1) + 100*(b+1) +c;
				int Jt2Tag = 1000*(a+2) +100*(b+1) +c;
				int xBeamTag = 1000*(a+1)+ 100*(b+1)  + c ;
				
				ID MemberInfo =  ID(3);
				MemberInfo(0)= a+1 ; MemberInfo(1)= b+1; MemberInfo(2) = c;
	
				SIFXBeam *xBeam = new SIFXBeam(xBeamTag, 1, Jt1Tag, Jt2Tag, MemberInfo);
				//xBeam->assignSection(theSIFDomain->getSIFSection(1));
				
				//Xbeam locates in the (a+1)XBay, the b bay grid, the c storey;
				theSIFDomain->addSIFXBeam(xBeam);  

				int theCompID =0 ; SIFCompartment* theComp=0;
				if(b==0){
					theCompID = (a+1)*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddXBeam(xBeamTag);
				}
				else if(b==NumYBays){
					theCompID = (a+1)*1000+b*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddXBeam(xBeamTag);
				}
				else{
					theCompID = (a+1)*1000+b*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddXBeam(xBeamTag);
					//Now adding the adjacent compartment
					theCompID = (a+1)*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddXBeam(xBeamTag);
				}
				if(BuilderType ==1||BuilderType ==11||BuilderType ==12||BuilderType ==13){
				//Add connected slab to the beam
					int slabID = theComp->getConnectedSlabs()(0);
					if(slabID!=0)
						xBeam->setConnectedSlab(slabID);
				}
	
				 //Here you have to update the function in SIFBuilderDomain, and rememeber to add mapofTaggedObjects for different class storage.
			}
		}
	}
	}

	if(BuilderType == 15)
		return TCL_OK;

	//for SIFybeam generation;

	if(BuilderType ==12||BuilderType ==14){								

	}
	else{
		if(BuilderType ==13||BuilderType ==16){								
			aEnd=0;
		}
		else{								
			aEnd=NumXBays;
		}

		for	(int b=0; b<NumYBays; b++)			
		{
		for	(int c=1; c<=NumStoreys; c++)		
		{
			for (int a=0; a<=aEnd; a++)	
			{
				int Jt1Tag = 1000*(a+1)+100*(b+1)+c;
				int Jt2Tag = 1000*(a+1)+100*(b+2)+c;
				int yBeamTag = 1000*(a+1)+100*(b+1)+c;

				ID MemberInfo =  ID(3);
				MemberInfo(0)= a+1 ; MemberInfo(1)= b+1 ;MemberInfo(2) = c;
				//Xbeam locates in the (a+1)XBay, the b bay grid, the c storey;

				SIFYBeam *yBeam = new SIFYBeam(yBeamTag, 1, Jt1Tag, Jt2Tag, MemberInfo);
				theSIFDomain->addSIFYBeam(yBeam);
				//yBeam->assignSection(theSIFDomain->getSIFSection(1));
				//Refine later for assigning section to the members
				
				int theCompID =0 ; SIFCompartment* theComp=0;
				if(a==0){
					theCompID = (a+1)*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddYBeam(yBeamTag);
				}
				else if(a==NumXBays){
					theCompID = a*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddYBeam(yBeamTag);
				}
				else{
					theCompID = a*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddYBeam(yBeamTag);
					//Now adding the adjacent compartment
					theCompID = (a+1)*1000+(b+1)*100+c;
					theComp = theSIFDomain->getSIFCompartment(theCompID);
					theComp->AddYBeam(yBeamTag);
				}	

				if(BuilderType ==1||BuilderType ==11||BuilderType ==12||BuilderType ==13){
				//Add connected slab to the beam
					int slabID = theComp->getConnectedSlabs()(0);
					if(slabID!=0)
						yBeam->setConnectedSlab(slabID);
				}

			}
		}
	}
	}

	

	if(BuilderType == 11||BuilderType ==16)
		return TCL_OK;


    //for SIFcolumn generation;
	if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
		opserr << endln << "SIFColumn: " << endln;

	if(BuilderType ==12){								
		bEnd=0;
	}
	else{								
		bEnd=NumYBays;
	}

	if(BuilderType ==14){								
		aEnd=0; bEnd=0;
	}
	else if(BuilderType ==13){																
			aEnd=0;
	}
	else{								
		aEnd=NumXBays; 
	}

	for	(int b=0; b<=bEnd; b++)			
	{
		for	(int c=1; c<=NumStoreys; c++)		
		{
			for (int a=0; a<=aEnd; a++)	
			{
				int Jt1Tag = 1000*(a+1)+100*(b+1)+c-1;
				int Jt2Tag = 1000*(a+1)+100*(b+1)+c;
				int columnTag = 1000*(a+1)+100*(b+1)+c;

				ID MemberInfo =  ID(3);
				MemberInfo(0)= a+1 ; MemberInfo(1)= b+1; MemberInfo(2) = c;
				//Xbeam locates in the (a+1)XBay, the b bay grid, the c storey;

				SIFColumn *column = new SIFColumn(columnTag, 1, Jt1Tag, Jt2Tag,MemberInfo);
				theSIFDomain->addSIFColumn(column);
				//column->assignSection(theSIFDomain->getSIFSection(1));

				ID theCompID = ID() ; SIFCompartment* theComp=0;
				
				if((a!=0)&&(b!=0)&&(a!=NumXBays)&&(b!=NumYBays)){
					theCompID.resize(4);
					theCompID(0) = a*1000+b*100+c;
					theCompID(1) = (a+1)*1000+b*100+c;
					theCompID(2) = a*1000+(b+1)*100+c;
					theCompID(3) = (a+1)*1000+(b+1)*100+c;	
				}
				else{
					theCompID.resize(1);
					if((a==0)&&(b==0))
						theCompID(0)=  (a+1)*1000+(b+1)*100+c;
					else if((a==0)&&(b==NumYBays))
						theCompID(0)=  (a+1)*1000+b*100+c;
					else if((a==NumXBays)&&(b==0))
						theCompID(0)=  a*1000+(b+1)*100+c;
					else if((a==NumXBays)&&(b==NumYBays))
						theCompID(0)=  a*1000+b*100+c;
					else{
						theCompID.resize(2);
						if(a==0){
							theCompID(0)=  (a+1)*1000+b*100+c;
							theCompID(1)=  (a+1)*1000+(b+1)*100+c;
						}
						else if(a==NumXBays){
							theCompID(0)=  a*1000+b*100+c;
							theCompID(1)=  a*1000+(b+1)*100+c;
						}
						else if(b==0){
							theCompID(0)=  a*1000+(b+1)*100+c;
							theCompID(1)=  (a+1)*1000+(b+1)*100+c;
						}
						else if(b==NumYBays){
							theCompID(0)=  a*1000+b*100+c;
							theCompID(1)=  (a+1)*1000+b*100+c;
						}
					}
				}
				int NumComps=theCompID.Size();
				for(int k=0;k<NumComps;k++){
					int theCompTag = theCompID(k);
					theComp = theSIFDomain->getSIFCompartment(theCompTag);
					theComp->AddColumn(columnTag);
					
					if (SIFdisplayindex == 1 || SIFdisplayindex == 2)
						opserr << "  "<<columnTag;

				}

			}
		}
	}
  
	//end of generating SIFModel
	opserr <<endln;
	return TCL_OK;

}

//-------------------------------------------Build finite element model for structure-----------------------------------------
int 
TclSIFBuilderCommand_BuildModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv){
	
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }

	bool isElastic= false;
	bool isDynamic = false;
	bool isPinned = false;
	int NumCtrlx, NumCtrly, NumCtrlz;
	
	int count;
	count =1;
	if (argc>1){
		if (strcmp(argv[count],"Elastic") == 0 ||strcmp(argv[count],"elastic") == 0 ) {
		isElastic = true; 
		count++;
		}
		else if (strcmp(argv[count],"Dynamic") == 0 ||strcmp(argv[count],"dynamic") == 0 ) {
		isElastic = false; 
		isDynamic = true;
		count++;
		}
		else if (strcmp(argv[count],"DynamicElastic") == 0 ||strcmp(argv[count],"dynamicElastic") == 0 ) {
		isElastic = true; 
		isDynamic = true;
		count++;
		}
	}
	//Here define the mesh control for structural model
	if(isElastic){
		NumCtrlx =2 ; NumCtrly =2 ;NumCtrlz =2 ;
		}
	else {
		NumCtrlx =6 ; NumCtrly =6 ;NumCtrlz =6 ;
		}
	//customise mesh control
	if(argc-count>0){
		if (strcmp(argv[count],"MeshCtrl") == 0 ||strcmp(argv[count],"-MeshCtrl") == 0 ) {
		count++;
		
		if(argc-count<3){
			opserr<<"WARNING TclSIFBuilderCommand expects 3 arguments for mesh control"<<endln;
			return TCL_ERROR;
		}
		
		if (Tcl_GetInt(interp, argv[count], &NumCtrlx) != TCL_OK) {
			opserr << "WARNING:: invalid number control " << argv[count] << " for building model \n";	
			return TCL_ERROR;
		}
		count++;
		if (Tcl_GetInt(interp, argv[count], &NumCtrly) != TCL_OK) {
			opserr << "WARNING:: invalid number control " << argv[count] << " for building model \n";	
			return TCL_ERROR;
		}
		count++;
		if (Tcl_GetInt(interp, argv[count], &NumCtrlz) != TCL_OK) {
			opserr << "WARNING:: invalid number control " << argv[count] << " for building model \n";	
			return TCL_ERROR;
		}
		count++;
		}
		//end of recieving mesh control		
	}
	if (argc - count > 0) {
		if (strcmp(argv[count], "pinned") == 0 || strcmp(argv[count], "-pinned") == 0) {
			isPinned = true;
			count++;
		}
		//end of detecting if the beams are pinned to the columns
	}

	
	ID MeshPars = ID(3);
	MeshPars(0)= NumCtrlx;MeshPars(1)= NumCtrly;MeshPars(2)= NumCtrlz;
	 
	//Now generate structural model following the check of SIFModel
	if (SIFModelStatus == 0) {
		if (TclSIFBuilderCommand_BuildSIFModel() == TCL_OK) {
			SIFModelStatus++;
		}
		else
			opserr << "WARNING::SIFModel has not been built before adding secondary beams" << endln;
	}
	
	if(SIFModelStatus>0){
		//Now generate structural model
		theSIFDomain->GenStructuralModel(MeshPars,isElastic, isDynamic, isPinned);	
    }
	

return TCL_OK;

}


//-----------------Add SIFSecBeam----------------
int
TclSIFBuilderCommand_addSecBeam(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{

	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - HTMaterial\n";
    return TCL_ERROR;
    }

	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before adding secondary beams"<<endln;
   }

	int theCompTag;
	int SecBeamType;
	int numBeams;
	int count = 1; 
	int SectionTag;

	if (strcmp(argv[count],"XBeam") == 0 ||strcmp(argv[count],"xbeam") == 0 ) {
		SecBeamType = 1;
		count++;
	}
	else if (strcmp(argv[count],"ZBeam") == 0 ||strcmp(argv[count],"zbeam") == 0 ) {
		SecBeamType = 2;
		count++;
	}
	else 
		opserr<<"SIFBuilder::Secondary beam tag"<<argv[count]<<endln;

	
	if (Tcl_GetInt(interp, argv[count], &SectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid beam section for Adding SIFSecBeams: " << argv[count] << "\n";	
		return TCL_ERROR;
	}
	else
		count++;


	if (strcmp(argv[count],"Compartment") == 0 ||strcmp(argv[count],"-compartment") == 0 ||strcmp(argv[count],"-Compartment") == 0 ) {
		count++;
		if (Tcl_GetInt(interp, argv[count], &theCompTag) != TCL_OK) {
		opserr << "WARNING:: invalid compartment tag for Adding SIFSecBeams: " << argv[count] << "\n";	
		return TCL_ERROR;
		}
		count++;
	}

	if (strcmp(argv[count],"NumBeams") == 0 ||strcmp(argv[count],"-numBeams") == 0 ||strcmp(argv[count],"-numbeams") == 0 ) {
		count++;
		if (Tcl_GetInt(interp, argv[count], &numBeams) != TCL_OK) {
		opserr << "WARNING:: invalid number of beams for Adding SIFSecBeams: " << argv[count] << "\n";	
		return TCL_ERROR;
		};
		count++;
	}
	
	SIFSection* theSection = theSIFDomain->getSIFSection(SectionTag);

	SIFXBeamSec* theSecBeam = 0;
	SIFCompartment* theCompartment = theSIFDomain->getSIFCompartment(theCompTag);
	ID theCompInfo = theCompartment->getCompartmentInfo();
	

	ID XBeamTags = theCompartment->getConnectedXBeams();
#ifdef _DEBUG
	opserr<<"Compartment "<<theCompTag<<" XBeam "<<XBeamTags;
#endif
	SIFXBeam* theXBeam1 = theSIFDomain->getSIFXBeam(XBeamTags(0));
	SIFXBeam* theXBeam2 = theSIFDomain->getSIFXBeam(XBeamTags(1));
	SIFJoint* theJT1 = theSIFDomain->getSIFJoint(theXBeam1->getConnectedJoints()(0));
	SIFJoint* theJT2 = theSIFDomain->getSIFJoint(theXBeam1->getConnectedJoints()(1));

	SIFJoint* theJT3 = theSIFDomain->getSIFJoint(theXBeam2->getConnectedJoints()(0));
	SIFJoint* theJT4 = theSIFDomain->getSIFJoint(theXBeam2->getConnectedJoints()(1));

	Vector theCrds1 = theJT1->getCrds();
	Vector theCrds2 = theJT2->getCrds();
	Vector theCrds3 = theJT3->getCrds();
	Vector theCrds4 = theJT4->getCrds();

	int Dim = theCrds1.Size();
	

//opserr<<theCrds1<<theCrds2<<theCrds3<<theCrds4;

	for(int j=0; j<numBeams; j++){
		ID theMemInfo = ID(4);
		for(int i=0; i<3; i++)
			theMemInfo(i) = theCompInfo(i);

		theMemInfo(3) = j+1;
		double Crds1[3]; double Crds2[3];
		
		for(int k=0;k<Dim;k++){
			Crds1[k] = theCrds1(k) + (j+1.0)/(numBeams+1)*(theCrds3(k)-theCrds1(k));
			Crds2[k] = theCrds2(k) + (j+1.0)/(numBeams+1)*(theCrds4(k)-theCrds2(k));
		}
		
		SIFJoint* theNewJoint1 = new SIFJoint(theNewJtTag++, Crds1[0], Crds1[1], Crds1[2]);
		SIFJoint* theNewJoint2 = new SIFJoint(theNewJtTag++, Crds2[0], Crds2[1], Crds2[2]);
		
		theSIFDomain->addSIFJoint(theNewJoint1);
		theSIFDomain->addSIFJoint(theNewJoint2);

		if(SecBeamType==1){
			SecXBeamTag++;
			theSecBeam = new SIFXBeamSec(SecXBeamTag, 2, theNewJtTag-2,theNewJtTag-1 ,theMemInfo,j);
		}
		else if(SecBeamType ==2){
		 //SecBeamTag++;
		 //theSecBeam = new SIFYBeamSec(SecBeamTag,2, theCompTag, theMemInfo);
		}
		theSecBeam->assignSection(theSection);
		if(BuilderType ==1){
			//Add connected slab to the beam
			int slabID = theCompartment->getConnectedSlabs()(0);
			if(slabID!=0)
			theSecBeam->setConnectedSlab(slabID);
		}
		theCompartment->AddXBeamSec(SecXBeamTag);
		theSIFDomain->addSIFXBeamSec(theSecBeam);
		//opserr<<theMemInfo;
	}

	
	return TCL_OK;
}

//-----------------delete member----------------
int
TclSIFBuilderCommand_delMember(ClientData clientData, Tcl_Interp *interp, int argc,
	TCL_Char **argv)
{

	if (theSIFDomain == 0) {
		opserr << "WARNING no active SIFBuilder Domain - HTMaterial\n";
		return TCL_ERROR;
	}

	if (SIFModelStatus == 0) {
		if (TclSIFBuilderCommand_BuildSIFModel() == TCL_OK) {
			SIFModelStatus++;
		}
		else
			opserr << "WARNING::SIFModel has not been built before adding secondary beams" << endln;
	}

	int theCompTag;
	int Membertype;
	int numBeams;
	int count = 1;
	int SectionTag;

	if (strcmp(argv[count], "column") == 0 || strcmp(argv[count], "column") == 0) {
		Membertype = 1;
		count++;
	}
	else if (strcmp(argv[count], "xbeam") == 0 || strcmp(argv[count], "XBeam") == 0) {
		Membertype = 2;
		count++;
	}
	else if (strcmp(argv[count], "zbeam") == 0 || strcmp(argv[count], "ZBeam") == 0) {
		Membertype = 3;
		count++;
	}
	else
		opserr << "SIFBuilder::wrong member type tag" << argv[count] << endln;

	//select members based on ranges of bays and storeys
	if (strcmp(argv[count], "xbay") == 0 || strcmp(argv[count], "-xbay") == 0 || strcmp(argv[count], "-XBay") == 0) {
		count++;
		if (Tcl_GetInt(interp, argv[count], &theCompTag) != TCL_OK) {
			opserr << "WARNING:: invalid compartment tag for Adding SIFSecBeams: " << argv[count] << "\n";
			return TCL_ERROR;
		}
		count++;
	}

	return TCL_OK;


}

//-----------------Add SIFSection----------------
int
TclSIFBuilderCommand_addSection(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{

	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
	SIFSection* theSIFSection = 0;
	int SIFSectionTag = 0;
	int SIFSectionTypeTag = 0;
	int SIFSectionMatTag =0;
	int SIFSectionSecMatTag =0;
	Vector SectionPars = Vector();
	//if(argc<3){
	//opserr << "WARNING:: Creating HeatTransfer material requires at least 3 arguments." << "\n";	
	//return TCL_ERROR;
	//}
    int count = 1; 
    //AddMaterial steel 1 -type EC3 fy E0;
	if (strcmp(argv[count],"Rectangular") == 0 ||strcmp(argv[count],"Rect") == 0 ) {
		SIFSectionTypeTag = 1; 
		double xLoc =0; double yLoc=0;
		double Breadth = 0;
		double Height = 0;
		
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFSectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFSection Tag recieved
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFSectionMatTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
	    
		count++;

	   if(argc-count<2){
			opserr << "WARNING:Dimension should be speciefied for SIFSection "<<SIFSectionTag<<endln;
	   }
	   else{
	    if (Tcl_GetDouble (interp, argv[count], &Breadth) != TCL_OK) {
			opserr << "WARNING invalid Breadth" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;
		if (Tcl_GetDouble (interp, argv[count], &Height) != TCL_OK) {
			opserr << "WARNING invalid Height" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
	   }
    SIFMaterial* theSIFMaterial = theSIFDomain->getSIFMaterial(SIFSectionMatTag);

	theSIFSection = new SIFSection(SIFSectionTag,SIFSectionTypeTag,theSIFMaterial);

	SectionPars.resize(2);

	SectionPars(0) =  Breadth;
	SectionPars(1) = Height;

	if(theSIFSection->AssignSectionPars(SectionPars)!=0){
		opserr << "WARNING AssignSectionPars failed " << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
	}

	}
	else if(strcmp(argv[count],"ISection") == 0 ||strcmp(argv[count],"ISec") == 0 ) {
		SIFSectionTypeTag = 2; 
		double xLoc =0; double yLoc=0;
		double d = 0;			//nominal depth;
		double tw = 0;			//web thickness;
		double bf = 0;			//flange width;
		double tf = 0;			//flange thickness;
		double coat = 0;        //coating thickness;
		count++;

		if (Tcl_GetInt(interp, argv[count], &SIFSectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFSection Tag recieved
		count++;

		if (Tcl_GetInt(interp, argv[count], &SIFSectionMatTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
	    count++;

	   if(argc-count<4){
			opserr << "WARNING:Dimension should be speciefied for SIFSection "<<SIFSectionTag<<endln;
	   }
	   else{
	    if (Tcl_GetDouble (interp, argv[count], &d) != TCL_OK) {
			opserr << "WARNING invalid nominal depth" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &bf) != TCL_OK) {
			opserr << "WARNING invalid flange width" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &tw) != TCL_OK) {
			opserr << "WARNING invalid web thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &tf) != TCL_OK) {
			opserr << "WARNING invalid flange thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;
		if(argc>count){
			if(strcmp(argv[count],"protected") == 0 ||strcmp(argv[count],"-protected") == 0 ){
				SIFSectionTypeTag = 22; 
				count++;
				if (Tcl_GetDouble (interp, argv[count], &coat) != TCL_OK) {
				opserr << "WARNING invalid coating thickness" << endln;
				opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
				return TCL_ERROR;
				}
			}
		}

	 }
    SIFMaterial* theSIFMaterial = theSIFDomain->getSIFMaterial(SIFSectionMatTag);

	theSIFSection = new SIFSection(SIFSectionTag,SIFSectionTypeTag,theSIFMaterial);

	SectionPars.resize(4);

	SectionPars(0) = d;
	SectionPars(1) = tw;
	SectionPars(2) = bf;
	SectionPars(3) = tf;
	if(SIFSectionTypeTag==22){
		SectionPars.resize(5);
		SectionPars(0) = d;
		SectionPars(1) = tw;
		SectionPars(2) = bf;
		SectionPars(3) = tf;
		SectionPars(4) = coat;
	}

	if(theSIFSection->AssignSectionPars(SectionPars)!=0){
		opserr << "WARNING AssignSectionPars failed " << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
	}
	}
	//End of aading ISection
	else if(strcmp(argv[count],"Composite") == 0 ||strcmp(argv[count],"Comp")==0||strcmp(argv[count],"composite") == 0 ) {
		SIFSectionTypeTag = 3; 
		double xLoc =0; double yLoc=0;
		double d = 0;			//nominal depth;
		double tw = 0;			//web thickness;
		double bf = 0;			//flange width;
		double tf = 0;			//flange thickness;
		double slabW = 0;       //Slab width
		double slabT =0;        //slab thickness;
		count++;

		if (Tcl_GetInt(interp, argv[count], &SIFSectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFSection Tag recieved
		count++;

		if (Tcl_GetInt(interp, argv[count], &SIFSectionMatTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
	    count++;

		if(strcmp(argv[count],"-SecMat") == 0 ||strcmp(argv[count],"-secmat")==0){
			count++;
			if (Tcl_GetInt(interp, argv[count], &SIFSectionSecMatTag) != TCL_OK) {
				opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[count] << "\n";	
				return TCL_ERROR;
			}
			count++;
		 }
			
	   if(argc-count<6){
			opserr << "WARNING:Dimension should be speciefied for SIFSection "<<SIFSectionTag<<endln;
	   }
	   else{
	    if (Tcl_GetDouble (interp, argv[count], &d) != TCL_OK) {
			opserr << "WARNING invalid nominal depth" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &bf) != TCL_OK) {
			opserr << "WARNING invalid flange width" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &tw) != TCL_OK) {
			opserr << "WARNING invalid web thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &tf) != TCL_OK) {
			opserr << "WARNING invalid flange thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;
		if (Tcl_GetDouble (interp, argv[count], &slabW) != TCL_OK) {
			opserr << "WARNING invalid slab width" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &slabT) != TCL_OK) {
			opserr << "WARNING invalid slab thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}

	   }
    SIFMaterial* theSIFMaterial = theSIFDomain->getSIFMaterial(SIFSectionMatTag);

	theSIFSection = new SIFSection(SIFSectionTag,SIFSectionTypeTag,theSIFMaterial);
	if(SIFSectionSecMatTag!=0){
		 SIFMaterial* theSecSIFMaterial = theSIFDomain->getSIFMaterial(SIFSectionMatTag);
		 theSIFSection->assignSecondMaterial(theSecSIFMaterial);
		
		}

	SectionPars.resize(6);

	SectionPars(0) = d;
	SectionPars(1) = tw;
	SectionPars(2) = bf;
	SectionPars(3) = tf;
	SectionPars(4) = slabW;
	SectionPars(5) = slabT;

	if(theSIFSection->AssignSectionPars(SectionPars)!=0){
		opserr << "WARNING AssignSectionPars failed " << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
	}
	}
	//end of adding composite section
	else if (strcmp(argv[count],"SlabSection") == 0 ||strcmp(argv[count],"slabSec") == 0||strcmp(argv[count],"slab") == 0) {
		SIFSectionTypeTag = 10; 
		double xLoc =0;
		double Thickness = 0;
		double Width = 0;
		
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFSectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFSection Tag recieved
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFSectionMatTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFSection: " << argv[2] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
	    
		count++;

	   if(argc-count<1){
			opserr << "WARNING:Dimension should be speciefied for SIFSection "<<SIFSectionTag<<endln;
	   }
	   else if(argc-count==1){
	    if (Tcl_GetDouble (interp, argv[count], &Thickness) != TCL_OK) {
			opserr << "WARNING invalid mositure" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
		}
		count++;
	   }
	   else if(argc-count==2){
	    if (Tcl_GetDouble (interp, argv[count], &Thickness) != TCL_OK) {
			opserr << "WARNING invalid slab thickness" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		count++;
		  if (Tcl_GetDouble (interp, argv[count], &Width) != TCL_OK) {
			opserr << "WARNING invalid slab width" << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
			}
		  count++;
     }//IF INPUT ==8

    SIFMaterial* theSIFMaterial = theSIFDomain->getSIFMaterial(SIFSectionMatTag);

	theSIFSection = new SIFSection(SIFSectionTag,SIFSectionTypeTag,theSIFMaterial);

	SectionPars.resize(2);

	SectionPars(0) =  Thickness;
	SectionPars(1) =  Width;

	if(theSIFSection->AssignSectionPars(SectionPars)!=0){
		opserr << "WARNING AssignSectionPars failed " << endln;
			opserr << " for SIFBuilder to add SIFSection: " << argv[2] << endln;	    
			return TCL_ERROR;
	}

	}
	//end of adding slab section
	if(theSIFSection!=0){
		theSIFDomain->addSIFSection(theSIFSection);
	}
	else 
		opserr<<"WARNING: TclHTModule fail to add SIFSection: "<<argv[2]<<endln;

	return TCL_OK;

}

//----------------Assign SIFSection----------------
int 
TclSIFBuilderCommand_assignSection(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }
	
	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		//for SIFxbeam mesh;
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before assigning sections"<<endln;
   }
	SIFXBeam *theSIFXBeam;
	SIFXBeamIter &theSIFXBeams = theSIFDomain->getSIFXBeams();

	SIFYBeam *theSIFYBeam;
	SIFYBeamIter &theSIFYBeams = theSIFDomain->getSIFYBeams();
	int SIFSectionTag;
	int count = 1; 
	SIFSection* theSection;

	if (Tcl_GetInt(interp, argv[count], &SIFSectionTag) != TCL_OK) {
		opserr << "WARNING:: invalid Section tag for assigning section: " << argv[1] << "\n";
		return TCL_ERROR;
	}
	theSection = theSIFDomain->getSIFSection(SIFSectionTag);
	count++;

    //AddMaterial steel 1 -type EC3 fy E0;
	if (strcmp(argv[count],"Beams") == 0 || strcmp(argv[count], "AllBeams") == 0 ||strcmp(argv[count],"beams") == 0 ) {
		count++;
		while((theSIFXBeam = theSIFXBeams())!=0){
			theSIFXBeam->assignSection(theSection);
		}
  
		
		while((theSIFYBeam = theSIFYBeams())!=0){
			theSIFYBeam->assignSection(theSection);
		}

	}
	else if (strcmp(argv[count],"XBeam") == 0 ||strcmp(argv[count],"Xbeam") == 0 ) {
		count++;

		int XBayTag = 0; int ZBayTag = 0; int XBeamTag = 0;
	
		//specify bays
		if(argc>count){
			if (strcmp(argv[count],"-XBay") == 0 ||strcmp(argv[count],"XBay") == 0 ) {
				count++;
				if (Tcl_GetInt(interp, argv[count], &XBayTag) != TCL_OK) {
					opserr << "WARNING:: invalid Section tag for assigning section: " << argv[1] << "\n";	
					return TCL_ERROR;
				}
				count++;
			}
			else if (strcmp(argv[count], "-ZBay") == 0 || strcmp(argv[count], "ZBay") == 0) {
				count++;
				if (Tcl_GetInt(interp, argv[count], &ZBayTag) != TCL_OK) {
					opserr << "WARNING:: invalid Section tag for assigning section: " << argv[1] << "\n";
					return TCL_ERROR;
				}
				count++;
			}
			else if (strcmp(argv[count], "-all") == 0 || strcmp(argv[count], "-All") == 0) {
				count++;
				XBayTag = 0;
				XBeamTag = 0;
			}
			else if (Tcl_GetInt(interp, argv[count], &XBeamTag) != TCL_OK) {
				opserr << "WARNING:: invalid XBeam tag to assign section: " << argv[1] << "\n";
				return TCL_ERROR;
			}
		}
		//end of xBay

		if (XBeamTag != 0) {
			theSIFXBeam = theSIFDomain->getSIFXBeam(XBeamTag);
			theSIFXBeam->assignSection(theSection);
		}
		else {
			while ((theSIFXBeam = theSIFXBeams()) != 0) {
				if (XBayTag != 0) {
					if (theSIFXBeam->getMemberInfo()(0) == XBayTag) {
						theSIFXBeam->assignSection(theSection);
					}
				}
				else if (ZBayTag != 0) {
					if (theSIFXBeam->getMemberInfo()(1) == ZBayTag) {
						theSIFXBeam->assignSection(theSection);
					}
				}
				else
					theSIFXBeam->assignSection(theSection);

			}
		}
	}
	else if (strcmp(argv[count],"-ZBeam") == 0 ||strcmp(argv[count],"Zbeam") == 0 ) {
		count++;

		int XBayTag = 0; int ZBayTag = 0; int ZBeamTag = 0;
	
		//specify bays
		if (argc>count) {
			if (strcmp(argv[count], "-ZBay") == 0 || strcmp(argv[count], "ZBay") == 0) {
				count++;
				if (Tcl_GetInt(interp, argv[count], &ZBayTag) != TCL_OK) {
					opserr << "WARNING:: invalid Section tag for assigning section: " << argv[1] << "\n";
					return TCL_ERROR;
				}
				count++;
			}
			else if (strcmp(argv[count], "-XBay") == 0 || strcmp(argv[count], "XBay") == 0) {
				count++;
				if (Tcl_GetInt(interp, argv[count], &XBayTag) != TCL_OK) {
					opserr << "WARNING:: invalid Section tag for assigning section: " << argv[1] << "\n";
					return TCL_ERROR;
				}
				count++;
			}
			else if (strcmp(argv[count], "-all") == 0 || strcmp(argv[count], "-All") == 0) {
				count++;
				ZBayTag = 0;
				ZBeamTag = 0;
			}
			else if (Tcl_GetInt(interp, argv[count], &ZBeamTag) != TCL_OK) {
				opserr << "WARNING:: invalid ZBeam tag to assign section: " << argv[1] << "\n";
				return TCL_ERROR;
			}
		}
		//end of xBay
		
		if (ZBeamTag != 0) {
			theSIFYBeam = theSIFDomain->getSIFYBeam(ZBeamTag);
			theSIFYBeam->assignSection(theSection);
		}
		else {
			while ((theSIFYBeam = theSIFYBeams()) != 0) {
				if (ZBayTag != 0) {
					if (theSIFYBeam->getMemberInfo()(1) == ZBayTag) {
						theSIFYBeam->assignSection(theSection);
					}
				}
				else if (XBayTag != 0) {
					if (theSIFYBeam->getMemberInfo()(1) == XBayTag) {
						theSIFYBeam->assignSection(theSection);
					}
				}
				else  
					theSIFYBeam->assignSection(theSection);

			}
		}
	}
	else if (strcmp(argv[count],"Column") == 0 ||strcmp(argv[count],"column") == 0 ) {
		SIFColumn *theSIFColumn;
		SIFColumnIter &theSIFColumns = theSIFDomain->getSIFColumns();

		count++;

		int ColumnTag = 0;
		
		//specify bays
		if (argc>count) {
			if (strcmp(argv[count], "-all") == 0 || strcmp(argv[count], "-All") == 0) {
				count++;
				ColumnTag = 0;
			}
			else if (Tcl_GetInt(interp, argv[count], &ColumnTag) != TCL_OK) {
				opserr << "WARNING:: invalid XBeam tag to assign section: " << argv[1] << "\n";
				return TCL_ERROR;
			}
		}
		//end of column

		if (ColumnTag != 0) {
			theSIFColumn = theSIFDomain->getSIFColumn(ColumnTag);
			theSIFColumn->assignSection(theSection);
		}
		else {
			while ((theSIFColumn = theSIFColumns()) != 0) {

				theSIFColumn->assignSection(theSection);

			}
		}

	}
	else if (strcmp(argv[count],"slab") == 0 ||strcmp(argv[count],"Slab") == 0 ) {
  
		SIFSlab *theSIFSlab;
		SIFSlabIter &theSIFSlabs  = theSIFDomain->getSIFSlabs();

		count++;

		int SlabTag = 0;

		//specify bays
		if (argc>count) {
			if (strcmp(argv[count], "-all") == 0 || strcmp(argv[count], "-All") == 0) {
				count++;
				SlabTag = 0;
			}
			else if (Tcl_GetInt(interp, argv[count], &SlabTag) != TCL_OK) {
				opserr << "WARNING:: invalid XBeam tag to assign section: " << argv[1] << "\n";
				return TCL_ERROR;
			}
		}
		//end of column

		if (SlabTag != 0) {
			theSIFSlab = theSIFDomain->getSIFSlab(SlabTag);
			theSIFSlab->assignSection(theSection);
		}
		else {
			while ((theSIFSlab = theSIFSlabs()) != 0) {

				theSIFSlab->assignSection(theSection);

			}
		}
		
	}

return 0;
}

//-----------------Add SIFMaterial----------------
int
TclSIFBuilderCommand_addMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	/*
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}
	*/

	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - HTMaterial\n";
    return TCL_ERROR;
    }
	
	SIFMaterial* theSIFMaterial=0;
	int SIFMaterialTag = 0;
	int SIFMaterialTypeTag = 0;

	//if(argc<3){
	//opserr << "WARNING:: Creating HeatTransfer material requires at least 3 arguments." << "\n";	
	//return TCL_ERROR;
	//}
    int count = 1; 
    //AddMaterial steel 1 -type EC3 fy E0;
	if (strcmp(argv[count],"steel") == 0 ||strcmp(argv[count],"Steel") == 0 ) {
		
		double fy = 0;
		double E0 = 0;
		
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFMaterialTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFMaterial: " << argv[1] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
		count++;
	
	   if (strcmp(argv[count],"-type") == 0) {

		count++;

		if(strcmp(argv[count], "EC3") ==0){
			SIFMaterialTypeTag = 130;
		}else if(strcmp(argv[count], "EC2Nh") ==0 ||strcmp(argv[count], "EC2NH") ==0) {
			SIFMaterialTypeTag = 121;
		}else if (strcmp(argv[count], "EC2NC") ==0 ||strcmp(argv[count], "EC2Nc") ==0) {
			SIFMaterialTypeTag = 122;
		}else if (strcmp(argv[count], "EC2X") ==0 ||strcmp(argv[count], "EC2x") ==0) {
			SIFMaterialTypeTag = 123;
		}else{
			opserr << "WARNING:invalid material type tag in SIFBuilder::AddMaterial "<<SIFMaterialTag<<endln;
		    return TCL_ERROR;
		}
		count++;
	    }
	   //end of recieving type tag
	   if(argc-count<2){
			opserr << "WARNING:fy and E0 should be speciefied for SIFBuilder::AddMaterial "<<SIFMaterialTag<<endln;
	   }
	   else{
	    if (Tcl_GetDouble (interp, argv[count], &fy) != TCL_OK) {
			opserr << "WARNING invalid fy" << endln;
			opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;	    
			return TCL_ERROR;
			}
		count++;
		if (Tcl_GetDouble (interp, argv[count], &E0) != TCL_OK) {
			opserr << "WARNING invalid E0" << endln;
			opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;	    
			return TCL_ERROR;
			}
	   }

	theSIFMaterial = new SIFMaterial(SIFMaterialTag,SIFMaterialTypeTag,fy, E0);
	}
	//End of aading steel material
	if (strcmp(argv[count], "stainlessSteel") == 0 || strcmp(argv[count], "StainlessSteel") == 0) {
		opserr << "adding stainless steel" << endln;
		int gradeTag = 0;
		double fy = 0;
		double E0 = 0;
		double fu = 0;

		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFMaterialTag) != TCL_OK) {
			opserr << "WARNING:: invalid material tag for Adding SIFMaterial: " << argv[1] << "\n";
			return TCL_ERROR;
		}
		//SIFMaterial Tag recieved
		count++;

		if (strcmp(argv[count], "-type") == 0) {

			count++;

			if (strcmp(argv[count], "Grade14571") == 0) {
				SIFMaterialTypeTag = 301;
			}
			
			count++;
		}
		//end of recieving type tag
		if (argc - count<2) {
			opserr << "WARNING:fy and E0 should be speciefied for SIFBuilder::AddMaterial " << SIFMaterialTag << endln;
		}
		else {
			if (Tcl_GetDouble(interp, argv[count], &fy) != TCL_OK) {
				opserr << "WARNING invalid fy" << endln;
				opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;
				return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &fu) != TCL_OK) {
				opserr << "WARNING invalid fu" << endln;
				opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;
				return TCL_ERROR;
			}
			count++;
			if (Tcl_GetDouble(interp, argv[count], &E0) != TCL_OK) {
				opserr << "WARNING invalid E0" << endln;
				opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;
				return TCL_ERROR;
			}
			
			
		}

		theSIFMaterial = new SIFMaterial(SIFMaterialTag, SIFMaterialTypeTag, fy, fu, E0);
		opserr << "stainless steel added" << endln;
	}
	//end of adding stainless steel
	else if (strcmp(argv[count],"concrete") == 0 ||strcmp(argv[count],"Concrete") == 0 ) {
		
		double moist = 0;
		double fc = 0;
		
		count++;
		if (Tcl_GetInt(interp, argv[count], &SIFMaterialTag) != TCL_OK) {
		opserr << "WARNING:: invalid material tag for Adding SIFMaterial: " << argv[1] << "\n";	
		return TCL_ERROR;
		}
	    //SIFMaterial Tag recieved
		count++;
	
	   if (strcmp(argv[count],"-type") == 0) {

		count++;

		if(strcmp(argv[count], "EC2") ==0){
			SIFMaterialTypeTag = 220;
		}else{
			opserr << "WARNING:invalid material type tag in SIFBuilder::AddMaterial "<<SIFMaterialTag<<endln;
		    return TCL_ERROR;
		}
		count++;
	    }
	   //end of recieving type tag
	   if(argc-count<2){
			opserr << "WARNING:fy and E0 should be speciefied for SIFBuilder::AddMaterial "<<SIFMaterialTag<<endln;
	   }
	   else{
	    if (Tcl_GetDouble (interp, argv[count], &moist) != TCL_OK) {
			opserr << "WARNING invalid moisture" << endln;
			opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;	    
			return TCL_ERROR;
			}
		count++;
		if (Tcl_GetDouble (interp, argv[count], &fc) != TCL_OK) {
			opserr << "WARNING invalid fc" << endln;
			opserr << " for SIFBuilder to add material: " << SIFMaterialTag << endln;	    
			return TCL_ERROR;
			}
	   }
//SIFMaterial::SIFMaterial(int tag, int MaterialTypeTag, 
						 //double fc, double epsc0, double fcu, double epscu, double rat, double ft,
						// double Ets, double moisture)
	   fc = (-1)*fc;
	   double epsc0 =-0.0025; double fcu=0.2*fc; double epscu=-0.02; 
	   double  rat=0.1; double ft=-0.1*fc; double Ets=ft/0.004;
	theSIFMaterial = new SIFMaterial(SIFMaterialTag,SIFMaterialTypeTag, fc, epsc0, fcu, epscu, rat, ft, Ets, moist);
	}
	//End of aading concrete material

	if(theSIFMaterial!=0){
		theSIFDomain->addSIFMaterial(theSIFMaterial);
	}
	else 
		opserr<<"WARNING: TclHTModule fail to add HeatTransfer Material: "<<argv[1]<<endln;

	return TCL_OK;

}


//------------------------AddXBay-------------------------------------------------
int
TclSIFBuilderCommand_addXBay(ClientData clientData, Tcl_Interp *interp, int argc,
			 TCL_Char **argv)
{
	/*
	if (theTclHTModule == 0) {
	opserr << "WARNING current HeatTransfer Module has been destroyed - HTMaterial\n";    
    return TCL_ERROR;
	}
	*/
  
  double originLocX = 0;
  
  int numXBays = argc; //The number of xBays specifies the bounds of Bays, so +1;
  XBayVec = new Vector(numXBays);
  
  double data=0;
  
  (*XBayVec)(0) = originLocX;  //The origin could be adjusted if necessary;
  
  double currentLoc = originLocX;
  
  for (int i=1;i<numXBays; i++) {
    if(Tcl_GetDouble(interp, argv[i], &data)==TCL_OK){
      currentLoc = currentLoc+ data;
      (*XBayVec)(i) = currentLoc;
    }
    else
    {
      opserr << "WARNING:: invalid data "<< argv[i] << " for defining XBay: " << "\n";
      return TCL_ERROR;
    }
      
  }
  

	return TCL_OK;

}

//------------------------------------AddYBay
int
TclSIFBuilderCommand_addYBay(ClientData clientData, Tcl_Interp *interp, int argc,
			 TCL_Char **argv)
{
	
  double originLocY = 0;
  
  int numYBays = argc; //The number of YBays specifies the bounds of Bays, so +1;
  
  YBayVec = new Vector(numYBays);
  
  double data=0;
  
  (*YBayVec)(0) = originLocY;  //The origin could be adjusted if necessary;
  
  double currentLoc = originLocY;
  
  for (int i=1;i<numYBays; i++) {
    if(Tcl_GetDouble(interp, argv[i], &data)==TCL_OK){
      currentLoc = currentLoc+ data;
      (*YBayVec)(i) = currentLoc;
    }
    else
    {
      opserr << "WARNING:: invalid data "<< argv[i] << " for defining YBay: " << "\n";
      return TCL_ERROR;
    }
    
  }
  
  

	return TCL_OK;

}


int
TclSIFBuilderCommand_addSkew_Angle(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Skew_Angle\n";
    return TCL_ERROR;
    }
	
	return TCL_OK;

}

int
TclSIFBuilderCommand_addStorey(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	
  double originLocZ = 0;
  
  int numStoreys = argc;
  
  StoreyVec = new Vector(numStoreys);
  
  double data=0;
  
  (*StoreyVec)(0) = originLocZ;
  double currentLoc = originLocZ;
  
  for (int i=1;i<numStoreys; i++) {
    if(Tcl_GetDouble(interp, argv[i], &data)==TCL_OK){
      currentLoc = currentLoc+ data;
      (*StoreyVec)(i) = currentLoc;
    }
    else
    {
      opserr << "WARNING:: invalid data "<< argv[i] << " for defining Storey: " << "\n";
      return TCL_ERROR;
    }
    
  }
  
  return TCL_OK;

}



int
TclSIFBuilderCommand_SetBC(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	double MinusInf = -65555;
	ID theJointsID = ID();
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }
	
	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		//for SIFxbeam mesh;
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before defining boundary conditions"<<endln;
   }


	int BCTypeTag;
	double ylocLB = MinusInf;
	double ylocUB = MinusInf;
	double xlocLB = MinusInf;
	double xlocUB = MinusInf;
	double zlocLB = MinusInf;
	double zlocUB = MinusInf;
	int count = 1;

	
	if (strcmp(argv[count],"pinnedJoint") == 0 ||strcmp(argv[count],"-pinnedJoint") == 0||strcmp(argv[count],"-Pinned") == 0 ) {
		count++;
		BCTypeTag =1;
	}
	else if (strcmp(argv[count],"fixedJoint") == 0 ||strcmp(argv[count],"-fixedJoint") == 0||strcmp(argv[count],"-Fixed") == 0 ) {
		count++;
		BCTypeTag =2;
	}
	else {
		opserr<<"WARNING:: TclSIFBuilderCommand_SetBC can not recognise the BCType "<< argv[count]<<endln;
	}
		
	if(argc-count<=0){
         opserr<< "WARNING:: TclSIFBuilderCommand_SetBC wants locx, or locy , or locz   "<<endln;
		 return TCL_ERROR;
	}

if(argc-count>0){
  if(strcmp(argv[count],"-Locx") == 0||strcmp(argv[count],"Locx") == 0||strcmp(argv[count],"locx") == 0){
		count++;
		//-----for geting uncertain number of double data.
    
		if (Tcl_GetDouble (interp, argv[count], &xlocLB) != TCL_OK) {
			opserr << "WARNING invalid xloc" << endln;
			opserr << " for finding SIFJoint: " << endln;
			return TCL_ERROR;
		}
		count++;
		if(argc-count>0){
			if (Tcl_GetDouble (interp, argv[count], &xlocUB) == TCL_OK) {
			count++;
			}
			else
			xlocUB=xlocLB;
		}
		else
			xlocUB=xlocLB;

		  if (xlocUB<xlocLB) {
			 opserr << "WARNING::TclSIFBuilder getting xlocUB " << xlocUB <<" not greater than xlocLb "<<xlocLB;
			 return TCL_ERROR;
    }
  }
}
 // end of locx
  if(argc-count>0){
	if(strcmp(argv[count],"-Locy") == 0||strcmp(argv[count],"Locy") == 0||strcmp(argv[count],"locy") == 0){
			count++;
	//-----for geting uncertain number of double data.
    
		if (Tcl_GetDouble (interp, argv[count], &ylocLB) != TCL_OK) {
			opserr << "WARNING invalid yloc" << endln;
			opserr << " for finding SIFJoint: " << endln;
			return TCL_ERROR;
		}
		count++;
		if(argc-count>0){
			if (Tcl_GetDouble (interp, argv[count], &ylocUB) == TCL_OK) {
			count++;
			}
			else
			ylocUB=ylocLB;
		}
		else
			ylocUB=ylocLB;

		 if (ylocUB<ylocLB) {
			 opserr << "WARNING::TclSIFBuilder getting ylocUB " << ylocUB <<" not greater than ylocLb "<<ylocLB;
			 return TCL_ERROR;
		}
	}
  }
 // end of locy
  if(argc-count>0){
	 if(strcmp(argv[count],"-Locz") == 0||strcmp(argv[count],"Locz") == 0||strcmp(argv[count],"locz") == 0){
		count++;
		//-----for geting uncertain number of double data.
    
		if (Tcl_GetDouble (interp, argv[count], &zlocUB) != TCL_OK) {
			opserr << "WARNING invalid zloc" << endln;
			opserr << " for finding SIFJoint: " << endln;
			return TCL_ERROR;
		}
		count++;
		if(argc-count>0){
			if (Tcl_GetDouble (interp, argv[count], &zlocLB) == TCL_OK) {
			count++;
			}
			else
			zlocLB=zlocUB;
		}
		else
			zlocLB=zlocUB;

		if (zlocUB<zlocLB) {
			 opserr << "WARNING::TclSIFBuilder getting zlocUB " << zlocUB <<" not greater than zlocLb "<<zlocLB;
			 return TCL_ERROR;
		}
	}
  }
	
 // end of locx
  if((xlocUB == MinusInf)&&(ylocUB == MinusInf)&&(zlocUB == MinusInf)){
		 opserr<<"WARNING:: SIFBUilder : setBC failed to find joints!"<<endln;
	 return TCL_ERROR;
	  }

    
	SIFJoint *theSIFJoint;
	SIFJointIter &theSIFJoints  = theSIFDomain->getSIFJoints();

	while((theSIFJoint = theSIFJoints())!=0){
		Vector theJtCrds = theSIFJoint->getCrds();
        //may have multiple location requirements

		int countIn =0;
		if(xlocLB-MinusInf>1e-5){
			//if Locx LB specified
			if(xlocLB-(1e-5)<theJtCrds(0)&&theJtCrds(0)<xlocUB+(1e-5))
				countIn++;
		} else 
			countIn++;
		if(ylocLB-MinusInf>1e-5){
			if(ylocLB-(1e-5)<theJtCrds(1)&&theJtCrds(1)<ylocUB+(1e-5))
				countIn++;
		} else 
			countIn++;
		if(zlocLB-MinusInf>1e-5){
			if(zlocLB-(1e-5)<theJtCrds(2)&&theJtCrds(2)<zlocUB+(1e-5))
				countIn++;
		}
		else 
			countIn++;
		
		if(countIn ==3){
			theSIFJoint->setBCtype(BCTypeTag);
		#ifdef _DEBUG 
			//opserr<< "TclSIFBuilder::A joint with tag "<< theSIFJoint->getTag()<< " has been selected"<<endln;
		#endif
		}


	}
	//End of determining selected joints
	
	
	return TCL_OK;
    //This function will define basic fire model type with tag, its linked compartment id, 


}




///////////////////////////////////----------PartDamage///////////////////////

int
TclSIFBuilderCommand_AddPartDamage(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{

//SIFfireAction -type standard -compartment 1 - duration 3600; 
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }


	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		//for SIFxbeam mesh;
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before adding partial damage to coatings"<<endln;
   }

	int count = 1;
	ID theColumnID = ID();
	double DamageLoc1=0;
	double DamageLoc2=0;

	
	if(strcmp(argv[count],"Column") == 0 ||strcmp(argv[count],"-Column") == 0||strcmp(argv[count],"-column") == 0 ) {
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
	    theColumnID.resize(ArgEnd-ArgStart);
		if (ArgStart != ArgEnd) {
			for (int i=ArgStart; i<ArgEnd; i++) {
			Tcl_GetInt(interp, argv[i], &data);
			theColumnID(i-ArgStart) = data;
			}
		}
		
	}
	if(strcmp(argv[count],"Damage") == 0 ||strcmp(argv[count],"-damage") == 0||strcmp(argv[count],"-damage") == 0 ) {
		count++;
		
		if (Tcl_GetDouble (interp, argv[count], &DamageLoc1) != TCL_OK) {
        opserr << "WARNING invalid damage location 1" << endln;
        opserr << " for SIFBuilder to add partial damage: "  << endln;
        return TCL_ERROR;
		}
		count++;

		if (Tcl_GetDouble (interp, argv[count], &DamageLoc2) != TCL_OK) {
        opserr << "WARNING invalid damage location 2" << endln;
        opserr << " for SIFBuilder to add partial damage: "  << endln;
        return TCL_ERROR;
		}
		count++;

	}


	//Adding the load info to slab
	 if(theColumnID.Size()!=0){
		int theSIFColumnTag;  SIFColumn *theSIFColumn;
		for(int i = 0; i<theColumnID.Size();i++){
			theSIFColumnTag = theColumnID(i);
			theSIFColumn  = theSIFDomain->getSIFColumn(theSIFColumnTag);
			if(theSIFColumn==0){
				opserr<<"WARNING::TclSIFBuider failed to find SIFColumn with tag "<<theSIFColumnTag<<" for adding Partial damage"<<endln;
				return TCL_ERROR;
			}
			else{
				//setPartialDamage(double lowerLoc, double upperLoc);
				Vector damageVec = Vector(2);
				damageVec(0)= DamageLoc1; damageVec(1)= DamageLoc2;
				theSIFColumn->setPartialDamage(damageVec) ;
			}

		}

	 }


}

///////////////////////////////////----------PartDamage///////////////////////




int
TclSIFBuilderCommand_AddLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{

//SIFfireAction -type standard -compartment 1 - duration 3600; 
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }


	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		//for SIFxbeam mesh;
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before adding loads"<<endln;
   }

	int count = 1;
	ID theJointID = ID();
	ID theXBeamID = ID();
	ID theYBeamID = ID();
	ID theSlabID = ID();
	Vector theLoadVec = Vector(6);  //fo 3D model, this may be changed to 3 for 2D model

	if (strcmp(argv[count],"Joint") == 0 ||strcmp(argv[count],"-Joint") == 0||strcmp(argv[count],"-joint") == 0 ) {
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
	    theJointID.resize(ArgEnd-ArgStart);
		if (ArgStart != ArgEnd) {
			for (int i=ArgStart; i<ArgEnd; i++) {
			Tcl_GetInt(interp, argv[i], &data);
			theJointID(i-ArgStart) = data;
			}
		}
	
      // for geting uncertain number of doubel values
	
	}
	//end of obtaining joint ID
	else if(strcmp(argv[count],"Member") == 0 ||strcmp(argv[count],"-Member") == 0||strcmp(argv[count],"-member") == 0 ) {
		count++;
		
		if(strcmp(argv[count],"Slabs") == 0 ||strcmp(argv[count],"AllSlabs") == 0||strcmp(argv[count],"allslabs") == 0 ) {
			count++;
			SIFSlab *theSIFSlab;
			SIFSlabIter &theSIFSlabs  = theSIFDomain->getSIFSlabs();
			int size = 0;
			while((theSIFSlab = theSIFSlabs())!=0){
				size++;
				theSlabID.resize(size);
				theSlabID(size-1)=theSIFSlab->getTag();
			}

		}
		else if(strcmp(argv[count],"Beams") == 0 ||strcmp(argv[count],"AllBeams") == 0||strcmp(argv[count],"allbeams") == 0 ) {
			count++;
			SIFXBeam *theSIFXBeam;
			SIFXBeamIter &theSIFXBeams  = theSIFDomain->getSIFXBeams();
			int size = 0;
			while((theSIFXBeam = theSIFXBeams())!=0){
				size++;
				theXBeamID.resize(size);
				theXBeamID(size-1)=theSIFXBeam->getTag();
			}
			SIFYBeam *theSIFYBeam;
			SIFYBeamIter &theSIFYBeams  = theSIFDomain->getSIFYBeams();
			size = 0;
			while((theSIFYBeam = theSIFYBeams())!=0){
				size++;
				theYBeamID.resize(size);
				theYBeamID(size-1)=theSIFYBeam->getTag();
			}

		}

	}
	 
	//end of obtaining member ID

	 if (strcmp(argv[count],"-load") == 0 ||strcmp(argv[count],"-Load") == 0||strcmp(argv[count],"Load") == 0 ) {
		count++;
		
		int ArgStart = count;
		int ArgEnd = 0;
		double data;
		Vector LoadVec = Vector();
     
		while (count < argc && ArgEnd == 0) {
			if (Tcl_GetDouble(interp, argv[count], &data) != TCL_OK)
			ArgEnd = count;
			else
			count++;
		}
     
		if(count ==argc)
			ArgEnd = argc;

      //~ detecting the remianing number of input
	    LoadVec.resize(ArgEnd-ArgStart);
		if (ArgStart != ArgEnd) {
			for (int i=ArgStart; i<ArgEnd; i++) {
			Tcl_GetDouble(interp, argv[i], &data);
			LoadVec(i-ArgStart) = data;
			}
		}
     //if input less than 6
     if (LoadVec.Size()<6) {
       for (int i= 0; i<6; i++) {
		   if(i< LoadVec.Size())
			   theLoadVec(i)= LoadVec(i);
		   else
				theLoadVec(i) = 0;
   }
     }
#ifdef _DEBUG
		opserr<<"TCLSIFBuilder::AddLoad "<<theLoadVec<<endln;
#endif
	 }
	 //end of input for load

	//Adding the load info to the joint
	 if(theJointID.Size()!=0){
		 //if the size of theJointID is not zero, load info is to be stored in these joints
		int theSIFJointTag;  SIFJoint *theSIFJoint;
		for(int i = 0; i<theJointID.Size();i++){
			theSIFJointTag = theJointID(i);
			theSIFJoint  = theSIFDomain->getSIFJoint(theSIFJointTag);
			if(theSIFJoint==0){
				opserr<<"WARNING::TclSIFBuider failed to find SIFJoint with tag "<<theSIFJointTag<<" for adding Load"<<endln;
				return TCL_ERROR;
			}
			else
				theSIFJoint->addLoad(theLoadVec);
		}
	 }
   //Adding the load info to Xbeam
	 if(theXBeamID.Size()!=0){
		int theSIFXBeamTag;  SIFXBeam *theSIFXBeam;
		for(int i = 0; i<theXBeamID.Size();i++){
			theSIFXBeamTag = theXBeamID(i);
			theSIFXBeam  = theSIFDomain->getSIFXBeam(theSIFXBeamTag);
			if(theSIFXBeam==0){
				opserr<<"WARNING::TclSIFBuider failed to find SIFJoint with tag "<<theSIFXBeamTag<<" for adding Load"<<endln;
				return TCL_ERROR;
			}
			else{
				theSIFXBeam->addLoad(theLoadVec);
			}

		}

	 }
	 //adding the load for ybeam
	 if(theYBeamID.Size()!=0){
		int theSIFYBeamTag;  SIFYBeam *theSIFYBeam;
		for(int i = 0; i<theYBeamID.Size();i++){
			theSIFYBeamTag = theYBeamID(i);
			theSIFYBeam  = theSIFDomain->getSIFYBeam(theSIFYBeamTag);
			if(theSIFYBeam==0){
				opserr<<"WARNING::TclSIFBuider failed to find SIFJoint with tag "<<theSIFYBeamTag<<" for adding Load"<<endln;
				return TCL_ERROR;
			}
			else{
				theSIFYBeam->addLoad(theLoadVec);
			}

		}

	 }
   //Adding the load info to slab
	 if(theSlabID.Size()!=0){
		int theSIFSlabTag;  SIFSlab *theSIFSlab;
		for(int i = 0; i<theSlabID.Size();i++){
			theSIFSlabTag = theSlabID(i);
			theSIFSlab  = theSIFDomain->getSIFSlab(theSIFSlabTag);
			if(theSIFSlab==0){
				opserr<<"WARNING::TclSIFBuider failed to find SIFJoint with tag "<<theSIFSlabTag<<" for adding Load"<<endln;
				return TCL_ERROR;
			}
			else{
				theSIFSlab->addLoad(theLoadVec);
			}

		}

	 }

	 return TCL_OK;
	
}




int
TclSIFBuilderCommand_addFireAction(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	
	int FireTypeTag;
	ID CompartmentID = 0;
    double fireStart =0;
  double fireOriginx =0;
  double fireOriginy =0;
  double fireOriginz=0;
  double fireHRR = 1e6;
  double fireDia = 0.5;
  
  double thi = 0;
  double avent = 0; double hvent = 0;
  double atotal = 0; double afire = 0;
  double qfire = 0; double tlimt = 0;

	
	//SIFfireAction -type standard -compartment 1 - duration 3600; 
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }

	if(SIFModelStatus==0){
		if( TclSIFBuilderCommand_BuildSIFModel()==TCL_OK){
		//for SIFxbeam mesh;
		SIFModelStatus++;
		}
		else
			opserr<<"WARNING::SIFModel has not been built before adding fire action"<<endln;
   }

	int count = 1;

	if (strcmp(argv[count],"Compartment") == 0 ||strcmp(argv[count],"-Compartment") == 0||strcmp(argv[count],"-compartment") == 0 ) {
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
	    CompartmentID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          CompartmentID(i-ArgStart) = data;
        }
      }
      // for geting uncertain number of doubel values
	}
	//end of getting compartment

	if (strcmp(argv[count],"-type") == 0 ||strcmp(argv[count],"-Type") == 0||strcmp(argv[count],"Type") == 0 ) {
		count++;
		if (strcmp(argv[count],"standard") == 0 ||strcmp(argv[count],"Standard") == 0 ) {
			FireTypeTag=1;
		}
		else if (strcmp(argv[count],"Parametric") == 0 ||strcmp(argv[count],"parametric") == 0) {
			FireTypeTag=2;
		}
		else if (strcmp(argv[count], "Exponent") == 0 || strcmp(argv[count], "exp") == 0) {
			FireTypeTag = 21;
		}
		else if (strcmp(argv[count],"EC1Localised") == 0 ||strcmp(argv[count],"EC1Local") == 0) {
			FireTypeTag=3;
		}

		else{
			opserr<<"WARNING::TclSIFBuilderCommand_addFireAction wants fire type as "<<endln
				  <<"Standard, or Parametric, or HydroCarbon, or Localised" <<endln;
		}
		count++;
	}

    //if more arguments are detected, check if it is start time
	if(argc-count>0){
		if (strcmp(argv[count],"-start") == 0 ||strcmp(argv[count],"Start") == 0||strcmp(argv[count],"start") == 0 ) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &fireStart) != TCL_OK) {
				opserr << "WARNING invalid fire start time" << endln;
				opserr << " for SIFBuilder to add fire action: " << argv[count-2] << endln;	    
				return TCL_ERROR;
			}
			count++;
		}
	}
  
  //More arguments are detected
  
	if (argc - count > 0) {
		//for localised fire
		if (FireTypeTag == 3) {
			//information of fire origin
			if (strcmp(argv[count], "-origin") == 0 || strcmp(argv[count], "-Origin") == 0 || strcmp(argv[count], "Origin") == 0) {
				count++;
				if (Tcl_GetDouble(interp, argv[count], &fireOriginx) != TCL_OK) {
					opserr << "WARNING invalid fireOriginx" << endln;
					opserr << " for SIFBuilder to add fire action: " << endln;
					return TCL_ERROR;
				}
				count++;

				if (Tcl_GetDouble(interp, argv[count], &fireOriginy) != TCL_OK) {
					opserr << "WARNING invalid fireOriginy" << endln;
					opserr << " for SIFBuilder to add fire action: " << endln;
					return TCL_ERROR;
				}
				count++;

				if (Tcl_GetDouble(interp, argv[count], &fireOriginz) != TCL_OK) {
					opserr << "WARNING invalid fireOriginz" << endln;
					opserr << " for SIFBuilder to add fire action: " << endln;
					return TCL_ERROR;
				}
				count++;
			}
			//end of origin
			if (argc - count > 0) {
				if (strcmp(argv[count], "-HRR") == 0 || strcmp(argv[count], "HRR") == 0) {
					count++;
					if (Tcl_GetDouble(interp, argv[count], &fireHRR) != TCL_OK) {
						opserr << "WARNING invalid HRR" << endln;
						opserr << " for SIFBuilder to add fire action: " << endln;
						return TCL_ERROR;
					}
					count++;
				}
			}
			//end of HRR

			if (argc - count > 0) {
				if (strcmp(argv[count], "-diameter") == 0 || strcmp(argv[count], "-Diameter") == 0 || strcmp(argv[count], "-Dia") == 0) {
					count++;
					if (Tcl_GetDouble(interp, argv[count], &fireDia) != TCL_OK) {
						opserr << "WARNING invalid diameter" << endln;
						opserr << " for SIFBuilder to add fire action: " << endln;
						return TCL_ERROR;
					}
					count++;
				}
			}
			//end of Diameter
		}
		//end of localised fire
		else if (FireTypeTag == 2) {

		}
	}
    //end of recieving parameters
  

//------------------------for considering members affected by fires in adjacent compartments
	if(CompartmentID.Size()!=0){
		for (int i =0; i<CompartmentID.Size();i++){
		  int CompartmentTag = CompartmentID(i);

		  SIFCompartment* theCompartment = theSIFDomain->getSIFCompartment(CompartmentTag);

		  if(theCompartment==0){
			  opserr<<"WARNING:: invalid compartment tag "<<CompartmentTag<<" when defining fire"<<endln;
				return TCL_ERROR;
		  }
		  
      SIFfireAction* theSIFfireAction = new SIFfireAction(CompartmentTag, FireTypeTag, CompartmentTag);
	  //set the start time for fire;
	  if (fireStart != 0) {
		  theSIFfireAction->SetStartTime(fireStart);
	  }

      if(FireTypeTag==3){
        Vector fireOrigin = Vector(3);
        fireOrigin(0)= fireOriginx; fireOrigin(1)=fireOriginy;fireOrigin(2)=fireOriginz;
        theSIFfireAction->SetFireOrigin(fireOrigin);
        
        if (fireHRR!=0) {
          theSIFfireAction->SetFireHRR(fireHRR);
        }
        else{
          opserr << "WARNING invalid fireHRR 0" << endln;
          opserr << " for SIFBuilder to add fire action: "  << endln;
          return TCL_ERROR;
        }
        
        if (fireDia!=0) {
          theSIFfireAction->SetFireDia(fireDia);
        }
        else{
          opserr << "WARNING invalid fireDia 0" << endln;
          opserr << " for SIFBuilder to add fire action: "  << endln;
          return TCL_ERROR;
        }

      }
	

//-----------------------Here send the fire action to each members in the compartment--------------------	
		if(theSIFfireAction!=0){
			ID theSecXBeams = theCompartment->getConnectedSecXBeams();
			if(theSecXBeams!=0){
				int numSecXBeams = theSecXBeams.Size();
				for(int k=0;k<numSecXBeams;k++){
					int SecxBeamID= theSecXBeams(k);
					SIFXBeamSec* theSecXBeam = theSIFDomain->getSIFXBeamSec(SecxBeamID);
					theSecXBeam->AddFireAction(CompartmentTag);
				}
			}

			ID theXBeams = theCompartment->getConnectedXBeams();
			if(theXBeams!=0){
				int numXBeams = theXBeams.Size();
				for(int k=0;k<numXBeams;k++){
					int xBeamID= theXBeams(k);
					SIFXBeam* theXBeam = theSIFDomain->getSIFXBeam(xBeamID);
					theXBeam->AddFireAction(CompartmentTag);
				}
			}

			ID theYBeams = theCompartment->getConnectedYBeams();
			if(theYBeams!=0){
				int numYBeams = theYBeams.Size();
				for(int k=0;k<numYBeams;k++){
					int yBeamID= theYBeams(k);
					SIFYBeam* theYBeam = theSIFDomain->getSIFYBeam(yBeamID);
					theYBeam->AddFireAction(CompartmentTag);
				}
			}

			ID theColumns = theCompartment->getConnectedColumns();
			if(theColumns!=0){
				int numColumns = theColumns.Size();
				for(int k=0;k<numColumns;k++){
					int ColumnID= theColumns(k);
					SIFColumn* theColumn = theSIFDomain->getSIFColumn(ColumnID);
					theColumn->AddFireAction(CompartmentTag);
				}
			}
			ID theSlabs = theCompartment->getConnectedSlabs();
			if(theSlabs!=0){
				int numSlabs = theSlabs.Size();
				for(int k=0;k<numSlabs;k++){
					int SlabID= theSlabs(k);
					SIFSlab* theSlab = theSIFDomain->getSIFSlab(SlabID);
					theSlab->AddFireAction(CompartmentTag);
				}
			}
			
			theSIFDomain->addSIFfireAction(theSIFfireAction);
			theSIFfireAction->setSIFDomain(theSIFDomain);
			}
		  else {
			opserr<<"WARNING: TclSIFBuilderModule fail to add SIFfireAction: "<<argv[1]<<endln;
		    delete theSIFfireAction;
		  }
		}
	}
	
	return TCL_OK;
    //This function will define basic fire model type with tag, its linked compartment id, 
}


int
TclSIFBuilderCommand_addFirePars(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }
	
	//Here we will get fire duration and timestep
		
	
	return TCL_OK;

}



int
TclSIFBuilderCommand_applyAndanalyze(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
	
	if (theSIFDomain == 0) {
    opserr << "WARNING no active SIFBuilder Domain - Storey\n";
    return TCL_ERROR;
    }

    theSIFDomain->defineBC();
	int count = 1;
	int thePatternTag = 0;

	const char* pwd = getInterpPWD(interp);
    theSIFDomain->setInterpPWD(pwd);

	//apply gravity if "gravity" tag detected
	if (strcmp(argv[count],"SelfWeight") == 0 ||strcmp(argv[count],"-SelfWeight") == 0||strcmp(argv[count],"-selfweight") == 0 ||strcmp(argv[count],"-selfWeight") == 0 ) {
		count++;
		thePatternTag++;
		double dt = 1.0;
		if (strcmp(argv[count],"-dt") == 0 ||strcmp(argv[count],"dt") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &dt) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying gravity load in SIFBuilder: " << endln;
			return TCL_ERROR;
			}
			count++;
		}
				
		theSIFDomain->applySelfWeight(thePatternTag,dt);

	}
	
	if(argc-count<=0)
		return TCL_OK;

	if (strcmp(argv[count],"Load") == 0 ||strcmp(argv[count],"-load") == 0||strcmp(argv[count],"Miscload") == 0 ||strcmp(argv[count],"miscload") == 0 ) {
		count++;
		thePatternTag++;
		double dt = 1.0;
		if (strcmp(argv[count],"-dt") == 0 ||strcmp(argv[count],"dt") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &dt) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying load in SIFBuilder: " << endln;
			return TCL_ERROR;
			}
		count++;
		}
		
		
		theSIFDomain->applyMiscLoad(thePatternTag,dt);

	}

	if(argc-count<=0)
		return TCL_OK;

	if (strcmp(argv[count],"Fire") == 0 ||strcmp(argv[count],"-fire") == 0||strcmp(argv[count],"-Fire") == 0 ) {
		count++;
		thePatternTag++;
		double dt = 10; double duration = 3600;
		if (strcmp(argv[count],"-dt") == 0 ||strcmp(argv[count],"dt") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &dt) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying fire action in SIFBuilder: " << endln;
			return TCL_ERROR;
			}

			count++;
//detecting duration
			if (strcmp(argv[count],"-duration") == 0 ||strcmp(argv[count],"duration") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &duration) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying fire action in SIFBuilder: " << endln;
			return TCL_ERROR;
			}

			count++;
		   }
		}

		else if (strcmp(argv[count],"-duration") == 0 ||strcmp(argv[count],"duration") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &duration) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying fire action in SIFBuilder: " << endln;
			return TCL_ERROR;
			}

			count++;
//detecting dt
			if (strcmp(argv[count],"-dt") == 0 ||strcmp(argv[count],"dt") == 0) {
			count++;
			if (Tcl_GetDouble (interp, argv[count], &dt) != TCL_OK) {
			opserr << "WARNING invalid dt" << endln;
			opserr << " for applying fire action in SIFBuilder: " << endln;
			return TCL_ERROR;
			}

			count++;
			}
		}

  const char* theDir =0;
		if(argc-count>0){
			if (strcmp(argv[count],"-output") == 0 ||strcmp(argv[count],"output") == 0) {
			count++;
			if(argc-count>0){
				theDir = argv[count];
				theSIFDomain->setHTdir(theDir);
			}
			count++;
			}

		}
		
		
		theSIFDomain->applyFireAction(thePatternTag,dt, duration);
		

	}


	return TCL_OK;

}


int TclSIFBuilderCommand_AddSIFRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	ID SIFJointID = 0;
	ID SIFXBeamID = 0;
	ID SIFZBeamID = 0;
	ID SIFSecXBeamID = 0;
	ID SIFYBeamID =0;
	ID SIFSlabID  =0;
	ID SIFColumnID = 0;
	ID* recNodes = 0;
	ID theDofs =0;

	Recorder *theRecorder =0;
	TCL_Char *fileName =0;
	OPS_Stream *theOutputStream = 0;
	int RecorderType =0; //1:Node-Disp; 2: Member_Disp; 3: Member-EleEndForce
	int count = 1;
	if (strcmp(argv[count],"Joint") == 0) {
		count++;
		RecorderType =1;

	}
	else if(strcmp(argv[count],"Member") == 0) {
		count++;
		RecorderType =2;
	}
	//end of Joints


	if (strcmp(argv[count],"-file") == 0 ||strcmp(argv[count],"file") == 0) {
		count++;
		fileName = argv[count];
		const char *pwd = getInterpPWD(interp);
		simulationInfo.addOutputFile(fileName, pwd);
		theOutputStream = new DataFileStream(fileName);
		count++;
	}

	if (strcmp(argv[count],"-joint") == 0) {
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
	    SIFJointID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          SIFJointID(i-ArgStart) = data;
        }
      }
	}
	else if (strcmp(argv[count],"-xBeam") == 0) {
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
	    SIFXBeamID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          SIFXBeamID(i-ArgStart) = data;
        }
      }
	}
	else if (strcmp(argv[count], "-zBeam") == 0) {
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
		SIFZBeamID.resize(ArgEnd - ArgStart);
		if (ArgStart != ArgEnd) {
			for (int i = ArgStart; i<ArgEnd; i++) {
				Tcl_GetInt(interp, argv[i], &data);
				SIFZBeamID(i - ArgStart) = data;
			}
		}
	}
	else if (strcmp(argv[count], "-Column") == 0) {
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
		SIFColumnID.resize(ArgEnd - ArgStart);
		if (ArgStart != ArgEnd) {
			for (int i = ArgStart; i<ArgEnd; i++) {
				Tcl_GetInt(interp, argv[i], &data);
				SIFColumnID(i - ArgStart) = data;
			}
		}
	}
	else if (strcmp(argv[count],"-SecXBeam") == 0) {
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
	    SIFSecXBeamID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          SIFSecXBeamID(i-ArgStart) = data;
        }
      }
	}
	else if (strcmp(argv[count],"-slab") == 0) {
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
	    SIFSlabID.resize(ArgEnd-ArgStart);
      if (ArgStart != ArgEnd) {
        for (int i=ArgStart; i<ArgEnd; i++) {
          Tcl_GetInt(interp, argv[i], &data);
          SIFSlabID(i-ArgStart) = data;
        }
      }
	}

      // for geting uncertain number of doubel values
	if (strcmp(argv[count],"-disp") == 0 ||strcmp(argv[count],"disp") == 0) {
		theDofs = ID(3);
		theDofs(0) = 0;theDofs(1) = 1;theDofs(2) = 2;
		count++;
		if(RecorderType==1){
			recNodes = new ID(SIFJointID.Size());
			for(int i=0; i<SIFJointID.Size();i++){
				SIFJoint* theJt = theSIFDomain->getSIFJoint(SIFJointID(i));
				(*recNodes)(i) = theJt->getNodeTag()(0);
			}
		}
	}
	else if (strcmp(argv[count],"-Mideflect") == 0 ||strcmp(argv[count],"Mideflect") == 0) {
		theDofs = ID(3);
		theDofs(0) = 0; theDofs(1) = 1; theDofs(2) = 2;
			count++;
			RecorderType =2;
			if(SIFXBeamID.Size()!=0){
				recNodes = new ID(SIFXBeamID.Size());
				for(int i=0; i<SIFXBeamID.Size();i++){
					SIFXBeam* theXBeam = theSIFDomain->getSIFXBeam(SIFXBeamID(i));
					if(theXBeam==0){
						opserr<<"WARNING:: TclSIFBuilder detected XBeam "<<SIFXBeamID(i)<< " not existing for defining SIFRecorders"<<endln;
						return TCL_ERROR;
					}
					ID theIntNodes = theXBeam->getIntNodeTags();
					int MidTag = (theIntNodes.Size())/2;
					(*recNodes)(i) = theIntNodes(MidTag);
				}
			}
			if (SIFZBeamID.Size() != 0) {
				recNodes = new ID(SIFZBeamID.Size());
				for (int i = 0; i<SIFZBeamID.Size();i++) {
					SIFYBeam* theZBeam = theSIFDomain->getSIFYBeam(SIFZBeamID(i));
					if (theZBeam == 0) {
						opserr << "WARNING:: TclSIFBuilder detected XBeam " << SIFZBeamID(i) << " not existing for defining SIFRecorders" << endln;
						return TCL_ERROR;
					}
					ID theIntNodes = theZBeam->getIntNodeTags();
					int MidTag = (theIntNodes.Size()) / 2;
					(*recNodes)(i) = theIntNodes(MidTag);
				}
			}

			if (SIFColumnID.Size() != 0) {
				recNodes = new ID(SIFColumnID.Size());
				for (int i = 0; i<SIFColumnID.Size(); i++) {
					SIFColumn* theColumn = theSIFDomain->getSIFColumn(SIFColumnID(i));
					if (theColumn == 0) {
						opserr << "WARNING:: TclSIFBuilder detected Column " << SIFColumnID(i) << " not existing for defining SIFRecorders" << endln;
						return TCL_ERROR;
					}
					ID theIntNodes = theColumn->getIntNodeTags();
					int MidTag = (theIntNodes.Size()) / 2;
					(*recNodes)(i) = theIntNodes(MidTag);
				}
			}

			if(SIFSecXBeamID.Size()!=0){
				recNodes = new ID(SIFSecXBeamID.Size());
				for(int i=0; i<SIFSecXBeamID.Size();i++){
					SIFCompartment* theComp = theSIFDomain->getSIFCompartment(SIFSecXBeamID(i));
					int SecBeamTag = theComp->getConnectedSecXBeams()(0);
					SIFXBeamSec* theXBeamSec = theSIFDomain->getSIFXBeamSec(SecBeamTag);
					if(theXBeamSec==0){
						opserr<<"WARNING:: TclSIFBuilder detected SecXBeam "<<SIFSecXBeamID(i)<< " not existing for defining SIFRecorders"<<endln;
						return TCL_ERROR;
					}
					ID theIntNodes = theXBeamSec ->getIntNodeTags();
					int MidTag = (theIntNodes.Size())/2;
					(*recNodes)(i) = theIntNodes(MidTag);
				}
			}
			

			if(SIFSlabID.Size()!=0){
				recNodes = new ID(SIFSlabID.Size());
				for(int i=0; i<SIFSlabID.Size();i++){
					SIFSlab* theSlab = theSIFDomain->getSIFSlab(SIFSlabID(i));
					if(theSlab==0){
						opserr<<"WARNING:: TclSIFBuilder detected Slab "<<SIFSlabID(i)<< " not existing for when defining SIFRecorders"<<endln;
						return 0;
					}
					ID theIntNodes = theSlab->getIntNodeTags();
					int MidTag = (theIntNodes.Size())/2;
					(*recNodes)(i) = theIntNodes(MidTag);
				}
			}
		
			
	}
	else if (strcmp(argv[count],"-deflect") == 0 ||strcmp(argv[count],"deflect") == 0) {
		theDofs = ID(1);
		theDofs(0) = 1;
		count++;
		RecorderType =2;
		
		int SlabId = SIFSlabID(0);
		SIFSlab* theSlab = theSIFDomain->getSIFSlab(SlabId);
		if(theSlab==0){
			opserr<<"WARNING:: TclSIFBuilder detected Slab "<<SlabId<< " not existing for when defining SIFRecorders"<<endln;
			return 0;
		}
		ID theIntNodes = theSlab->getIntNodeTags();
		recNodes = new ID(theIntNodes.Size());
		*recNodes = theIntNodes;
	}
	else if (strcmp(argv[count],"-force") == 0 ||strcmp(argv[count],"force") == 0) {
			count++;
			if(RecorderType ==2)
				RecorderType =3;
	}



	 int sensitivity = 0;
	 TCL_Char *responseID = 0;
	 double dT = 0.0;
	 bool echoTimeFlag = true;
	 TimeSeries **theTimeSeries = 0;
	 if(recNodes ==0){
		opserr<<"WARNING:: no node defined for SIFRecorder"<<endln;
		 return TCL_ERROR;
	 }
#ifdef _DEBUG
	opserr<<"TCLSIFbuilder RecNode"<<*recNodes<<endln;
#endif
	theRecorder = new NodeRecorder(theDofs, 
					   recNodes, 
					   sensitivity,
					   responseID, 
					   *theDomain, 
					   *theOutputStream, 
					   dT, 
					   echoTimeFlag,
					   theTimeSeries);
	
	if (theRecorder == 0) {
      //sprintf(interp->result,"-1");
      return TCL_ERROR;
	}

  if ((theDomain->addRecorder(*theRecorder)) < 0) {
    opserr << "WARNING could not add to domain - recorder " << argv[1]<< endln;
    delete theRecorder;
    return TCL_ERROR;
  } 


	return 0;

}