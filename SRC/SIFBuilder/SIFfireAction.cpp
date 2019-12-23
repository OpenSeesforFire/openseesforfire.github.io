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
// This file constructs the SIFCompartment iterator.
// Created by Liming Jiang @UoE

//#SIFBUILDER PROJECT#//
//Created by Liming Jiang, UoE,2014
//-----------------------------------
//SIFfireAction will be created for compartment fire
//The function SIFfireAction::Apply will look for all the members linked to the comparment,and then define the  //
//boundary conditions, building up the heat transfer model. The SIFHTforMember will generate a PathtimeSeries to//
//define elemental thermal action, used in thermo-mechanical analysis.

//This class will probably work as tclHTmodelbuilder, taking the responsibilities to store the constants and environment parameters
//------------------------------------



#include <SIFfireAction.h>
#include <SIFCompartment.h>
#include <SIFBuilderDomain.h>
#include <SIFMember.h>
#include <NorminalFireEC1.h>
#include <ParametricFireEC1.h>
#include <LocalizedFireEC1.h>
#include <PathTimeSeriesThermal.h>
#include <Beam3dThermalAction.h>
#include <ThermalActionWrapper.h>
#include <NodalThermalAction.h>
#include <ShellThermalAction.h>
#include <Matrix.h>
#include <math.h>


//extern int theEleLoadTag;
//extern int theNodalLoadTag;
static int thePathSeriesTag =0;
SIFfireAction::SIFfireAction(int tag, int fireModelType, int compartmentID):TaggedObject(tag),
FireModelType(fireModelType), FireHRR(0), FireDia(0),StartTime(0),
theFireModel(0),FireDuration(0), TimeStep(0),FireOrigin(0),theHTDomain(0),theSIFDomain(0)
{
	
	if(compartmentID!=0){
		CompartmentID = compartmentID;
		
		//Here it is necessary to identify whether there are neighboring comparments 
		
	}else{
		opserr<<"WARNING:SIFfireAction "<< tag<< " has no compartment assigned "<<endln;
	}
  
  if (theHTDomain==0) {
	  theHTDomain = new HeatTransferDomain(); 
  } else {
    theHTDomain->clearAll();
  }

  FireOrigin = new Vector(3);
  FireOrigin->Zero();
  
  
	
}

SIFfireAction::~SIFfireAction()
{
  if (theHTDomain!=0) {
    theHTDomain->clearAll();
    delete theHTDomain;
  }
}

void
SIFfireAction::setSIFDomain(SIFBuilderDomain *theSifDomain){
	 theSIFDomain = theSifDomain;
	
}


int 
SIFfireAction::SetFireOrigin(const Vector& fireOrigin)
{
    (*FireOrigin) = fireOrigin;
	return 0;
}

int
SIFfireAction::SetFireDia(double fireDia)
{
  FireDia=fireDia;
  return 0;
}

int
SIFfireAction::SetFireHRR(double fireHRR)
{
  FireHRR = fireHRR;
  return 0;
}

	
int 
SIFfireAction::SetFirePath(const Matrix& FirePath)
{
	
	return 0;
	
	
}

int 
SIFfireAction::SetStartTime(double startTime)
{
	
	if (startTime!=0) {
    StartTime=startTime;
  }
  else {
    StartTime =0;
  }
	return 0;
}

int 
SIFfireAction::SetFireDuration(double fireDuration)
{
	
	if (fireDuration!=0) {
    FireDuration=fireDuration;
  }
  else {
    opserr<<"WARNING::SIFfireAction can't apply fire action with 0s duration"<<endln;
    return -1;
  }
	
	
}


//////////////////////////----------------------Update Fire Model-----------------------------//////////
int 
SIFfireAction::UpdateFireModel(int MemberType)
{
  //define and set up the fire model
  if(FireModelType==1){
	  if(theFireModel==0)
		theFireModel= new NorminalFireEC1(this->getTag(),1,StartTime);
      return 0;
	 }
  else if(FireModelType==2){
	  // ParametricFireEC1(int tag, double I, double Av, double H, double At, double Af, double Qf, double tlim, double startTime = 0);
	  if(theFireModel==0)
		 theFireModel= new ParametricFireEC1(this->getTag(), 1159, 28.13, 3.5,752,240, 420, 1200,StartTime);

	  return 0;
	 }
  else if (FireModelType == 21) {
	  // ParametricFireEC1(int tag, double I, double Av, double H, double At, double Af, double Qf, double tlim, double startTime = 0);
	  if (theFireModel == 0)
		  theFireModel = new NorminalFireEC1(this->getTag(), 5, StartTime);

	  return 0;
  }
  else if(FireModelType==3){
      int centreLineTag;
	  double fireLocx, fireLocy, fireLocz, D, Q;
      if(MemberType ==1){
		 centreLineTag=2;
		fireLocx = -(*FireOrigin)(2);
		fireLocy = (*FireOrigin)(1);
		fireLocz = (*FireOrigin)(0);
	  }
	  else if(MemberType ==21){
		  centreLineTag=2;
		fireLocx = -(*FireOrigin)(2);
		fireLocy = (*FireOrigin)(1);
		fireLocz = (*FireOrigin)(0);
	  }
      else if(MemberType ==2){
		  centreLineTag=2;
		fireLocx = (*FireOrigin)(0);
		fireLocy = (*FireOrigin)(1);
		fireLocz = (*FireOrigin)(2);
	  }
	  else if(MemberType ==3){
		  centreLineTag=3;
		fireLocx = -(*FireOrigin)(2);
		fireLocy = -(*FireOrigin)(0);
		fireLocz = (*FireOrigin)(1);
	  }
	  else if(MemberType ==10){
		  centreLineTag=3;
		fireLocx = (*FireOrigin)(2);
		fireLocy = (*FireOrigin)(0);
		fireLocz = (*FireOrigin)(1);
	  }

	  else 
		opserr<<"WARNING::SIFfireAction can not update fire model for MemberType: "<<MemberType;
     
#ifdef _DEBUG
	  opserr<<"sIFFIREAction::Updated fireOrigin "<<fireLocx <<", "<<fireLocy<<", "<<fireLocz<<endln;
#endif

    if (FireDia!=0) {
      D = FireDia;
    }
    else
    {
		opserr<<"WARNING::SIFfireAction failed to find the diameter for localised  fire "<<this->getTag()<<endln;
      return -1;
    }
    
    if (FireHRR!=0) {
      Q = FireHRR;
    }
    else
    {
      opserr<<"WARNING::SIFfireAction failed to find the fire load for localised fire "<<this->getTag()<<endln;
     return -1;
    }
	SIFCompartment* theComp = theSIFDomain->getSIFCompartment(CompartmentID);
	SIFColumn* theColumn = theSIFDomain->getSIFColumn(theComp->getConnectedColumns()(0));
	SIFJoint* theJt1 =theSIFDomain->getSIFJoint((theColumn->getConnectedJoints()(0)));
	SIFJoint* theJt2 =theSIFDomain->getSIFJoint((theColumn->getConnectedJoints()(1)));

	SIFSlab* theSlab = 0;
	double SlabT;
	double ColumnLen = theJt2->getCrds()(1)-theJt1->getCrds()(1);
	if(theComp->getConnectedSlabs()==0)
		SlabT=0.1;
	else{
		theSlab = theSIFDomain->getSIFSlab(theComp->getConnectedSlabs()(0));
		SlabT =  (theSlab->getSIFSectionPtr())->getSectionPars()(0);
	}
		
    double storeyH = ColumnLen- SlabT;

	if(storeyH<=0){
		opserr<<"WARNING::invalid storey height defined for localised fire model"<<endln;
		storeyH = 3.0;
	}

	if(theFireModel!=0){
		delete theFireModel;
		theFireModel=0;
	}

	theFireModel = new LocalizedFireEC1(this->getTag(),fireLocx, fireLocy, fireLocz, D , Q, storeyH, centreLineTag);

  }
  else{
     opserr<<"WARNING::SIFfireAction failed to initialise fire model for compartment "<<CompartmentID<<endln;
   }
  
  return 0;
  
}

//////////////////////////----------------------Apply fire action-----------------------------//////////
	
int 
SIFfireAction::Apply(int LoadPatternTag, double timeStep, double fireDuration )
{
  
//This function is to apply fire action which is associated with the department
//Looped operation will be carried out over the members within the comparment

  Domain* theDomain = theSIFDomain->getStructureDomain();
  int thePatternTag = LoadPatternTag;
  
  //check the fireDuration is not 0;
  if (fireDuration!=0) {
    FireDuration=fireDuration;
  }
  else {
    opserr<<"WARNING::SIFfireAction can't apply fire action with 0s duration"<<endln;
    return -1;
  }
  
  //check the timeStep is not 0;
  if (timeStep!=0) {
    TimeStep=timeStep;
  }
  else {
    opserr<<"WARNING::SIFfireAction can't apply fire action with 0s timeStep"<<endln;
    return -1;
  }
  
  
  SIFCompartment* theCompartment = theSIFDomain->getSIFCompartment(CompartmentID);
  if (theCompartment==0) {
    opserr<<"WARNING::SIFfireAction failed to allocate the pointer to the compartment "<<CompartmentID<< endln;
    return -1;
  }

 
  ID SecXBeams = theCompartment->getConnectedSecXBeams();
  ID XBeams = theCompartment->getConnectedXBeams();
  ID YBeams = theCompartment->getConnectedYBeams();
  ID Columns = theCompartment->getConnectedColumns();
  ID Slabs = theCompartment->getConnectedSlabs();
  int theMemberID;

  int dislabel = 0;
  
  //Then loop over all the SecXBeams in the compartment
  int NumSecXBeams = SecXBeams.Size();
  if(NumSecXBeams>0)
		opserr<<endln<<"Secondary XBeam  ";
  if (dislabel == 7)
	  NumSecXBeams = 0;
  for(int i=0; i<NumSecXBeams; i++){
    //Get the member pointer from SIFDomain;
    theMemberID = SecXBeams(i);
    SIFXBeamSec* theMember =  theSIFDomain->getSIFXBeamSec(theMemberID);
    opserr<<" "<<theMember->getTag();
    this->RunHTforMember(theMember, thePatternTag);
  }

  //Then loop over all the XBeams in the compartment
	

  int NumXBeams = XBeams.Size();
  if(NumXBeams>0)
		opserr<<endln<<"XBeam  ";
  if (dislabel == 1|| dislabel == 123)
	  NumXBeams = 0;
  for(int i=0; i<NumXBeams; i++){
    //Get the member pointer from SIFDomain;
    theMemberID = XBeams(i);
    SIFXBeam* theMember =  theSIFDomain->getSIFXBeam(theMemberID);
    opserr<<" "<<theMember->getTag();
    this->RunHTforMember(theMember, thePatternTag);
  }
  //end of loop over the XBeams
  

  //Then loop over all the YBeams in the compartmens
  int NumYBeams = YBeams.Size();
   if(NumYBeams>0)
		opserr<<endln<<"YBeam  ";
   if (dislabel == 2 || dislabel == 123)
	   NumYBeams = 0;
  for(int i=0; i<NumYBeams; i++){
    theMemberID = YBeams(i);
    SIFYBeam* theMember =  theSIFDomain->getSIFYBeam(theMemberID);
	opserr<<" "<<theMember->getTag();
    this->RunHTforMember(theMember, thePatternTag);
  }
  //end of loop over the YBeams
  
  
    //Then loop over all the Columns in the compartmens
 int NumColumns = Columns.Size();
 if (dislabel == 3 || dislabel == 123)
	 NumColumns = 0;
  if(NumColumns>0)
		opserr<<endln<<"Column  ";
  for(int i=0; i<NumColumns; i++){
    theMemberID = Columns(i);
    SIFColumn* theMember =  theSIFDomain->getSIFColumn(theMemberID);
	opserr<<" "<<theMember->getTag();
	    if(theMember->getTag()==3301)
		this->RunHTforMember(theMember, thePatternTag);
  }
  //end of loop over the Columns
  

   int NumSlabs = Slabs.Size();
   if (dislabel == 9)
	   NumSlabs = 0;
	 if(NumSlabs>0)
		opserr<<endln<<"Slab ";
  for(int i=0; i<NumSlabs; i++){
    theMemberID = Slabs(i);
    SIFSlab* theMember =  theSIFDomain->getSIFSlab(theMemberID);
	opserr<<" "<<theMember->getTag();
   this->RunHTforMember(theMember, thePatternTag);
  }
  //end of loop over the Slabs
 
	//Finnaly all the members will have thermalActions defined directly with pathtimeseriesThermal.
	return 0;
}
	


int 
SIFfireAction::RunHTforMember(SIFMember* theMember,int LoadPatternTag)
{
	 //double PartialDamage =0.5;   //Partial damage length for column


	 ID theModelInfo = theSIFDomain->getSIFBuilderInfo();
	int ModelType = theModelInfo(0);
    ID appliedFireActions = theMember->getAppliedFireActions();
	if(appliedFireActions==0)
		return 0;

	Domain* theDomain = theSIFDomain->getStructureDomain();
	Vector damageVec = 0;
	
    //Create a SIFHTforMember;
    //SIFHTforMember(int tag, SIFMember* theMember,HeatTransferDomain* theSIFDomain, double fireDuration, double timeSteps);
    SIFHTforMember* theHTforMember = new SIFHTforMember(theMember->getTag(),theMember, theSIFDomain,FireDuration,TimeStep);
    int MemberTypeTag = theMember->getMemberTypeTag();
	int theSectionType = (theMember->getSIFSectionPtr())->getSectionTypeTag();

    if(theHTforMember->SetHTDomain(theHTDomain)!=0){
      opserr<<"WARNING:SIFfireAction failed to set theHTDomain for SIFHTforMember "<<theMember->getTag()<<endln;
      return -1;
    }
    ////////////////////////////////////////////////MULTI//////FIRE////////////////////////////////////////
	if(appliedFireActions.Size()>1||ModelType>10){
		//If more than one fireActions are attached with  this member
		int theFireActionID = appliedFireActions(0);
		if(theFireActionID!=this->getTag()){
			opserr<<"WARNING::SIFfireAction is over defined for member "<< this->getTag();
			return -6;
		}
		ID theCompInfo = (theSIFDomain->getSIFCompartment(CompartmentID))->getCompartmentInfo();
		if(theCompInfo.Size()!=3){
			opserr<<"WARNING::SIFfireAction receives incorrect size of the compartment info "<< this->getTag();
			return -5;
		}
		ID theMemberInfo = theMember->getMemberInfo();
		
		ID* expoFaceID = new ID();
		ID AmbFaceID = ID(0);
		if(MemberTypeTag==1){
			//XBeam
			expoFaceID->resize(7);
				(*expoFaceID)(0)=1; 
				(*expoFaceID)(1)=4; (*expoFaceID)(2)=5; 
				(*expoFaceID)(3)=6; (*expoFaceID)(4)=7; (*expoFaceID)(5)=8; (*expoFaceID)(6)=9; 
			//end of yBayInfo == ComYBayInfo+1;
		}
		//end of XBeam
		else if(MemberTypeTag==21){
			
			expoFaceID->resize(9);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=5; 
				(*expoFaceID)(3)=6; (*expoFaceID)(4)=7; (*expoFaceID)(5)=8; (*expoFaceID)(6)=9; 
				(*expoFaceID)(7)=12; (*expoFaceID)(8)=13; 
			AmbFaceID.resize(1);
				AmbFaceID(0)=16;
		}
		//end of Secondary XBeam
		else if(MemberTypeTag==2){
				expoFaceID->resize(7);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=5; 
				(*expoFaceID)(3)=6; (*expoFaceID)(4)=7; (*expoFaceID)(5)=8; (*expoFaceID)(6)=9; 
		}
		//END OF yBeam
		else if(MemberTypeTag==3){
			expoFaceID->resize(8);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=5; 
				(*expoFaceID)(3)=6; (*expoFaceID)(4)=7; (*expoFaceID)(5)=8; (*expoFaceID)(6)=9; (*expoFaceID)(7)=12;
		}
		//END OF Column
		else if (MemberTypeTag == 10 || MemberTypeTag == 11) {
			expoFaceID->resize(1);
			(*expoFaceID)(0) = 1;
			AmbFaceID.resize(1);
			AmbFaceID(0) = 2;
		}
		//end of Slab
		this->UpdateFireModel(MemberTypeTag);
		theHTforMember->SetFireExposureCons(theFireModel,expoFaceID);
		theHTforMember->SetAmbExposureCons(AmbFaceID);
		theMember->WipeFireAction();
		
	}
	   ////////////////////////////////////////////////SINGLE//////FIRE////////////////////////////////////////
	else{
		//if only one fireAction is attached with this member
		int theFireActionID = appliedFireActions(0);
		if(theFireActionID!=this->getTag()){
			opserr<<"WARNING::SIFfireAction is over defined for member "<< this->getTag();
			return -6;
		}
		ID theCompInfo = (theSIFDomain->getSIFCompartment(CompartmentID))->getCompartmentInfo();
		if(theCompInfo.Size()!=3){
			opserr<<"WARNING::SIFfireAction receives incorrect size of the compartment info "<< this->getTag();
			return -5;
		}
		ID theMemberInfo = theMember->getMemberInfo();
		ID* expoFaceID = new ID();
		ID AmbFaceID = ID(0);
		if(MemberTypeTag==1){
			//XBeam
			int yBayInfo = theMemberInfo(1);
			int CompYBayInfo = theCompInfo(1);
			if(yBayInfo==CompYBayInfo){
			//Beam locates at the left side, wchich means exposure at right side
			if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=3; 
			}
			else if(theSectionType==2){
				expoFaceID->resize(4);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=6; (*expoFaceID)(3)=8; 
				AmbFaceID.resize(3);
				AmbFaceID(0)=5;AmbFaceID(1)=7;AmbFaceID(2)=9;;
			}
			}
			//end of yBayInfo==CompYBayInfo;
			else if(yBayInfo==CompYBayInfo+1){
			//Beam locates at the right side, wchich means exposure at left side
			if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=2; 
			}
			else if(theSectionType==2){
				expoFaceID->resize(4);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=5; (*expoFaceID)(2)=7; (*expoFaceID)(3)=9; 
				AmbFaceID.resize(3);
				AmbFaceID(0)=4;AmbFaceID(1)=6;AmbFaceID(2)=8;;
			}
			}

			//end of yBayInfo == ComYBayInfo+1;
		}
		//end of XBeam
		else if(MemberTypeTag==21){
			//SecXBeam
			int yBayInfo = theMemberInfo(1);
			int CompYBayInfo = theCompInfo(1);
			
			//Beam locates at the left side, wchich means exposure at right side
			if(theSectionType==1){
				expoFaceID->resize(3);
				(*expoFaceID)(0)=1;(*expoFaceID)(1)=2;  (*expoFaceID)(2)=3; 
			}
			else if(theSectionType==2){
				expoFaceID->resize(9);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4;(*expoFaceID)(2)=5; (*expoFaceID)(3)=6;
				(*expoFaceID)(4)=7; (*expoFaceID)(5)=8; (*expoFaceID)(6)=9;(*expoFaceID)(7)=12; (*expoFaceID)(8)=13;
				AmbFaceID.resize(1);
				AmbFaceID(0)=16;
			}

		}
		//end of XBeam
		else if(MemberTypeTag==2){
			//YBeam
			int xBayInfo = theMemberInfo(0);
			int CompXBayInfo = theCompInfo(0);
			if(xBayInfo==CompXBayInfo){
			//Beam locates at the left side, wchich means exposure at right side
			if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=3; 
			}
			else if(theSectionType==2){
				expoFaceID->resize(4);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=5; (*expoFaceID)(2)=7; (*expoFaceID)(3)=9; 
				AmbFaceID.resize(3);
				AmbFaceID(0)=4;AmbFaceID(1)=6;AmbFaceID(2)=8;;    
			}
			}
			//end of yBayInfo==CompYBayInfo;
			else if(xBayInfo==CompXBayInfo+1){
			//Beam locates at the right side, wchich means exposure at left side
			if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=2; 
			}
			else if(theSectionType==2){
				expoFaceID->resize(4);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=6; (*expoFaceID)(3)=8;  
				AmbFaceID.resize(3);
				AmbFaceID(0)=5;AmbFaceID(1)=7;AmbFaceID(2)=9;;
			}
			}
			//end of yBayInfo == ComYBayInfo+1;
			
		}
		//END OF yBeam
		else if(MemberTypeTag==3){
			//Column
			int xBayInfo = theMemberInfo(0);
			int yBayInfo = theMemberInfo(1);
			int CompXBayInfo = theCompInfo(0);
			int CompYBayInfo = theCompInfo(1);

			if(xBayInfo==CompXBayInfo&&yBayInfo==CompYBayInfo){
				//Column locates at the left top corner in the compartment, wchich means exposure at right down side
				if(xBayInfo==1){
					if(yBayInfo==1){
						//Corner column
					}
					else{
						//Side Columns
					}
				}
				//.....
				//Here we firstly consider all the columns locates as center aligned
				if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=3; 
				}
				else if(theSectionType==2||theSectionType==22){
				expoFaceID->resize(3);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=4; (*expoFaceID)(2)=6;
				AmbFaceID.resize(5);
				AmbFaceID(0)=5;AmbFaceID(1)=7;AmbFaceID(2)=8;AmbFaceID(3)=9;AmbFaceID(4)=12;
				}
			}
			//end of left top
			else if(xBayInfo==CompXBayInfo+1&&yBayInfo==CompYBayInfo){
			//Column locates at the right top corner in the compartment
				if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=2; 
				}
				else if(theSectionType==2||theSectionType==22){
				expoFaceID->resize(3);
				(*expoFaceID)(0)=12; (*expoFaceID)(1)=8; (*expoFaceID)(2)=6;   
				AmbFaceID.resize(5);
				AmbFaceID(0)=1;AmbFaceID(1)=4;AmbFaceID(2)=5;AmbFaceID(3)=7;AmbFaceID(4)=9;    
				}
			}
			//end of right top
			else if(xBayInfo==CompXBayInfo&&yBayInfo==CompYBayInfo+1){
			//Column locates at the left bottom corner in the compartment
				if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=2; 
				}
				else if(theSectionType==2||theSectionType==22){
				expoFaceID->resize(3);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=5; (*expoFaceID)(2)=7;   
				AmbFaceID.resize(5);
				AmbFaceID(0)=4;AmbFaceID(1)=6;AmbFaceID(2)=8;AmbFaceID(3)=9;AmbFaceID(4)=12;
				}
			}
			//end of left bottom
			else if(xBayInfo==CompXBayInfo+1&&yBayInfo==CompYBayInfo+1){
			//Column locates at the right bottom corner in the compartment
				if(theSectionType==1){
				expoFaceID->resize(2);
				(*expoFaceID)(0)=1; (*expoFaceID)(1)=2; 
				}
				else if(theSectionType==2||theSectionType==22){
				expoFaceID->resize(3);
				(*expoFaceID)(0)=12; (*expoFaceID)(1)=9; (*expoFaceID)(2)=7;   
				AmbFaceID.resize(5);
				AmbFaceID(0)=1;AmbFaceID(1)=4;AmbFaceID(2)=6;AmbFaceID(3)=8;AmbFaceID(4)=9;    
				}
			}
			//end of right bottom 
		}
		//END OF Column
		else if(MemberTypeTag==10||MemberTypeTag==11){
			expoFaceID->resize(1);
			(*expoFaceID)(0)=1;  
			AmbFaceID.resize(1);
			AmbFaceID(0)=2;

		}
		//end of Slab
		this->UpdateFireModel(MemberTypeTag);
		theHTforMember->SetFireExposureCons(theFireModel,expoFaceID);
		theHTforMember->SetAmbExposureCons(AmbFaceID);
	}
	//end of appliedFireActions;
	
    //----------------------Data transaction---------------------
    int NumofSeries;

    if(FireModelType ==1 ||FireModelType ==2 ||FireModelType == 21)
    {
      NumofSeries = 1;
    }
    else if(FireModelType ==3)
    {
      NumofSeries = 4;//may be 5
    }
    //Partial damage applied to central column;
	if(MemberTypeTag==3&&theMember->getPartialDamage(damageVec)>0 ){
		 NumofSeries = 2;
	}
  
    PathTimeSeriesThermal** thePathTimeSeries;
    thePathTimeSeries = new PathTimeSeriesThermal*[NumofSeries];
    
	for(int k=0;k<NumofSeries;k++){
		if(MemberTypeTag!=10 & MemberTypeTag!=11)
			thePathTimeSeries[k] = new PathTimeSeriesThermal(thePathSeriesTag++,15);
		else
			thePathTimeSeries[k] = new PathTimeSeriesThermal(thePathSeriesTag++,9);
	}
    theHTforMember->SetPathTimeSeries(thePathTimeSeries);
    
	//Vector CompOriginCrds = Vector(3);
	//CompOriginCrds = (theSIFDomain->getSIFCompartment(CompartmentID))->getCompartmentOrigin();
    
	//---------------********-------------------Key executing------------*******-------------------- 
	Matrix* CrdMat = new Matrix(NumofSeries,4);
	theHTforMember->applyFire(*FireOrigin, CrdMat);
	//---------------********-------------------Key executing------------*******-------------------- 
	Vector locs= Vector();
	locs =   theHTforMember->getRecLocations();
#ifdef _DEBUG
	opserr<<"SIFfireAction->getLocs "<<locs<<endln;
#endif
	//for uniform fire
	if(NumofSeries ==1){
		ID theMemberEleTags = theMember->getIntEleTags();
		if (theMemberEleTags ==0 ||theMemberEleTags.Size()==0){
			opserr<<"WARNING SIFfireAction can not find the element tags for SIFMember "<<theMember->getTag()<<endln;
		}
		ElementalLoad* theThermalAction=0;
		for(int i=0; i<theMemberEleTags.Size(); i++){
			if (MemberTypeTag != 10) {
				if (MemberTypeTag == 11) {
					theThermalAction = new Beam3dThermalAction(theSIFDomain->getEleLoadTag() + 1, locs, thePathTimeSeries[0], theMemberEleTags(i));
				}
				else
					theThermalAction = new Beam3dThermalAction(theSIFDomain->getEleLoadTag() + 1, locs(0), locs(1), locs(2), locs(3), thePathTimeSeries[0], theMemberEleTags(i));

			}
			else
				theThermalAction = new ShellThermalAction(theSIFDomain->getEleLoadTag()+1, locs(0),locs(8), thePathTimeSeries[0], theMemberEleTags(i));
			
			theSIFDomain->incrEleLoadTag();

			int loadPatternTag = LoadPatternTag;
			// add the load to the domain
			if (theDomain->addElementalLoad(theThermalAction, loadPatternTag) == false) {
			opserr << "WARNING eleLoad - could not add following load to domain:\n ";
			opserr << theThermalAction;
			delete theThermalAction;
			return -1;
			}
			theSIFDomain->incrEleLoadTag();
		}
	}	
	//non-uniform fire action
	else if(NumofSeries >1){
		ID theMemberEleTags = theMember->getIntEleTags();
		if (theMemberEleTags ==0 ||theMemberEleTags.Size()==0){
			opserr<<"WARNING SIFfireAction can not find the element tags for SIFMember "<<theMember->getTag()<<endln;
		}

		SIFJoint* theJT0 = theSIFDomain->getSIFJoint(theMember->getConnectedJoints()(0));
		int NodeEndTag0 = theJT0->getNodeTag()(0);
		int NodeEndTag1;
		if (MemberTypeTag != 10) {
			SIFJoint* theJT1 = theSIFDomain->getSIFJoint(theMember->getConnectedJoints()(1));
			NodeEndTag1 = theJT1->getNodeTag()(0);
		}
		else
		{
			SIFJoint* theJT1 = theSIFDomain->getSIFJoint(theMember->getConnectedJoints()(3));
			NodeEndTag1 = theJT1->getNodeTag()(0);
		}
		
		
		ID intNodes = theMember->getIntNodeTags();
		int MidNodeTag = intNodes((intNodes.Size()-1)*0.5);  //0.5!!!!!!

		int MidNodeTag1 = intNodes((intNodes.Size()-1)*0.3);  //0.3
		int MidNodeTag2 = intNodes((intNodes.Size()-1)*0.7);  //0.7

		NodalThermalAction** theNodalTA = new NodalThermalAction*[NumofSeries];
		//NodalThermalAction(int tag, int theNodeTag, const Vector& locy, TimeSeries* theSeries,Vector* crds ) 
		
		for(int i = 0; i<NumofSeries;i++){
			int NodeTag;
			if(i== 0)
				NodeTag = NodeEndTag0;
			else if(i==1){
				if(NumofSeries ==4)
					NodeTag = MidNodeTag1;
				else if(NumofSeries ==2)
					NodeTag = NodeEndTag1;
				else if(NumofSeries ==3)
					NodeTag = MidNodeTag;
				else
					opserr<<"SIFfireAction::incorrect number of series"<<endln;
			}
			else if(i==2){
				if(NumofSeries ==4)
					NodeTag = MidNodeTag2;
				else if(NumofSeries ==3)
					NodeTag = NodeEndTag1;
			}
			else if(i==3)
				NodeTag = NodeEndTag1;
			
			Vector* theCrds = new Vector(3);
			(*theCrds)(0) = (*CrdMat)(i,0);(*theCrds)(1) = (*CrdMat)(i,1);(*theCrds)(2) = (*CrdMat)(i,2);
#ifdef _DEBUG
			opserr<<"SIFfireAction "<<(*theCrds);
#endif		
			if(MemberTypeTag!=10)
				theNodalTA[i] = new NodalThermalAction(theSIFDomain->getNodalLoadTag()+1, NodeTag,locs(0),locs(1),locs(2),locs(3),thePathTimeSeries[i],theCrds);
			else 
				theNodalTA[i] = new NodalThermalAction(theSIFDomain->getNodalLoadTag()+1, NodeTag,locs,thePathTimeSeries[i],theCrds);
			delete theCrds;
			if (theDomain->addNodalLoad(theNodalTA[i], LoadPatternTag) == false) {
				 opserr << "WARNING TclModelBuilder - could not add load to domain\n";
				//printCommand(argc, argv);
				//delete theLoad;
			}
			theSIFDomain->incrNodalLoadTag();
		}

	
		for(int i=0; i<theMemberEleTags.Size(); i++){
			// ThermalActionWrapper(int tag, int EleTag , NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2, NodalThermalAction* theNodalTA3);
			ThermalActionWrapper* theThermalAction =0;
			// add the load to the domain
			Vector ratios = Vector(4);
			//Vector damageVec = 0;
			if(MemberTypeTag==3&&theMember->getPartialDamage(damageVec)>0 ){
				theThermalAction = new ThermalActionWrapper(theSIFDomain->getEleLoadTag()+1,theMemberEleTags(i),theNodalTA[0],theNodalTA[1]);
				ratios.resize(2);
				SIFJoint* theJT1 = theSIFDomain->getSIFJoint(theMember->getConnectedJoints()(0));
				SIFJoint* theJT2 = theSIFDomain->getSIFJoint(theMember->getConnectedJoints()(1));
				Vector crds1 = theJT1->getCrds();
				Vector crds2 = theJT2->getCrds();
				int crdSize = crds1.Size();
				double CrdInt = 0;
				for(int m = 0; m<crdSize;  m++){
					CrdInt += (crds1(m)- crds2(m))*(crds1(m)- crds2(m));
				}

				double theMemLength = sqrt(CrdInt);
				
				ratios(0) = (damageVec(0)-0.5)/theMemLength; ratios(1) = damageVec(0)/theMemLength; //Location of Partial damage!!!!!
			}
			else {
				if(NumofSeries==4){
				ratios.resize(4);
				theThermalAction = new ThermalActionWrapper(theSIFDomain->getEleLoadTag()+1,theMemberEleTags(i),theNodalTA[0],theNodalTA[1],theNodalTA[2],theNodalTA[3]);
				ratios(0) =0.0; ratios(1) =0.3;ratios(2) =0.7; ratios(3) =1.0;
				}
				else if(NumofSeries==3){
				ratios.resize(3);
				theThermalAction = new ThermalActionWrapper(theSIFDomain->getEleLoadTag()+1,theMemberEleTags(i),theNodalTA[0],theNodalTA[1],theNodalTA[2]);
				ratios(0) =0.0; ratios(1) =0.5; ratios(3) =1.0;
				}
				else 
					opserr<<"SIFfireAction::invalid number of series for dimensional reduced heat transfer analysis"<<endln;
			}
			theThermalAction->setRatios(ratios);
			if (theDomain->addElementalLoad(theThermalAction, LoadPatternTag) == false) {
			opserr << "WARNING eleLoad - could not add following load to domain:\n ";
			opserr << theThermalAction;
			delete theThermalAction;
			return -1;
			}
			theSIFDomain->incrEleLoadTag();
		}
	}

  delete thePathTimeSeries;
  return 0;
}




