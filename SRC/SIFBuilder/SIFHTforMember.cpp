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
// This file constructs the class SIFHTforMember which processes the heat 
// transfer inside a structure member.
// Created by Liming Jiang @UoE

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>

using std::ios;
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#include <SIFHTforMember.h>
#include <Simple_Boundary.h>
#include <FireImposedPattern.h>
#include <HTNodeRecorder.h>

#include <HT_AnalysisModel.h>
#include <CTestNormTempIncr.h>
#include <HT_SolutionAlgorithm.h>
#include <HT_TransientIntegrator.h>
#include <TemperatureBCHandler.h>
#include <RCM.h>
#include <HT_DOF_Numberer.h>
#include <BandGenLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <HT_TransientAnalysis.h>
#include <NewtonMethod.h>
#include <BackwardDifference.h>
#include <PenaltyBC_Handler.h>

 #include <StandardStream.h>
 #include <DataFileStream.h>

 #include <SimulationInformation.h>

#include <SIFXBeamSec.h>

#include <SimpleMaterial.h>

 extern SimulationInformation simulationInfo;

  static HT_AnalysisModel*        theHTModel = 0;
  static CTestNormTempIncr*       theHTTest= 0 ;
  static HT_SolutionAlgorithm*    theHTSolnAlgo= 0 ;
  static HT_TransientIntegrator*  theHTIntegrator=0;
  static TemperatureBCHandler*    theHTHandler = 0;
  static RCM*                     theHTRCM = 0;
  static HT_DOF_Numberer*         theHTNumberer= 0 ;
  static BandGenLinLapackSolver*  theHTSolver= 0;
  static BandGenLinSOE*             theHTSOE= 0;
  static HT_TransientAnalysis*     theHTAnalysis =0;



SIFHTforMember::SIFHTforMember(int tag, SIFMember* theMember, SIFBuilderDomain* theSifDomain, 
							   double fireDuration, double timeStep):
    TaggedObject(tag),theMember(theMember),theSection(0),theSection1(0), 
		theHTMaterial(0),theHTMaterial1(0),theSlabT(0),
    FireDuration(fireDuration),TimeStep(timeStep), FireActionIndex(0),
	theHTDomain(0),thePathSeries(0),RecLocations(0),theSIFDomain(theSifDomain)
{
	
	//Create Heat Transfer Materials for HT analyses. Material should be reset after the analysis.
	//HeatTransfer analysis need section type, dimension,material,fire model, exposure condition;
  
  if (theSifDomain==0) {
    opserr<<"WARNING:SIFHTforMember can't find theSIFBuilderDomain"<<endln;
  } else {
    theSIFDomain = theSifDomain;
  }
	
  theSection = theMember->getSIFSectionPtr();
  
  SIFMaterial* theMaterial = 0;
  SIFMaterial* theMaterial1 = 0;
  
  //theMember should be identified as composite beam or not
  //if yes, theMaterial2 should be obtained
  theMaterial = theSection->getSIFMaterialPtr();
  theHTMaterial=theMaterial->getHeatTransferMaterial();
  
  int theSlabTag = theMember->getConnectedSlab();
  SIFSlab* theSlab =0;
  if(theSlabTag!=0)
	theSlab = theSIFDomain->getSIFSlab(theSlabTag);
  
  if (theSlab!=0) {
    theSection1 = theSlab->getSIFSectionPtr();
    theMaterial1 = theSection1->getSIFMaterialPtr();
    theHTMaterial1 = theMaterial1->getHeatTransferMaterial();
    //thickness of the slab, material of the slab
  }
  
	theFireExpFacesID = new ID*[4];
	theFireModel=new FireModel*[4];

  if(theAmbExpFacesID!=0){
	theAmbExpFacesID=0;
  }
  
  theMesh = 0;
  theMesh1 = 0;
  
  theEntity = 0;
  theEntity1 = 0;
  
  
  
  //Analysis operators
  
  
}


SIFHTforMember::~SIFHTforMember()
{
	
  //------------------Heat Transfer Aanalysis--------------------------//
  if (theHTModel!=0) {
    delete theHTModel;
    theHTModel = 0;
  }
  if (theHTTest !=0) {
    delete theHTTest ;
    theHTTest  = 0;
  }
  if (theHTSolnAlgo !=0) {
    delete theHTSolnAlgo ;
    theHTSolnAlgo =0;
  }
  if (theHTIntegrator  !=0) {
    delete theHTIntegrator ;
    theHTIntegrator   = 0;
  }
  if (theHTHandler!=0) {
    delete theHTHandler ;
    theHTHandler  = 0;
  }
  if (theHTRCM!=0) {
    delete theHTRCM ;
    theHTRCM = 0;
  }
  if (theHTNumberer!=0) {
    delete theHTNumberer ;
    theHTNumberer = 0;
  }
  if (theHTSolver!=0) {
    delete theHTSolver ;
    theHTSolver = 0;
  }
  if (theHTSOE!=0) {
    delete theHTSOE  ;
    theHTSOE  = 0;
  }
  if (theHTAnalysis!=0) {
	delete theHTAnalysis;
	theHTAnalysis= 0;
  }
}



int
SIFHTforMember::SetHTDomain(HeatTransferDomain* theHTdomain)
{
  if (theHTdomain==0) {
    opserr<<"WARNING:SIFHTforMember can't set an empty Heat transfer Domain"<<endln;
    return -1;
  } else {
    theHTDomain = theHTdomain;
  }
  return 0;

}

/*
int 
SIFHTforMember::SetFireModel(FireModel* thefireModel)
{
	
  if (thefireModel==0) {
    opserr<<"WARNING:SIFHTforMember can't set an empty firemodel "<<thefireModel->getTag()<<endln;
    return -1;
  } else {
    theFireModel = thefireModel;
  }
  return 0;
}
*/
int
SIFHTforMember::SetFireExposureCons(FireModel* thefireModel, ID* theFireExposedfaces)
{
  
  
  if (thefireModel==0) {
    opserr<<"WARNING:SIFHTforMember defined empty fire model "<<thefireModel->getTag()<<endln;
    return -1;
  } 
  if (theFireExposedfaces==0) {
    opserr<<"WARNING:SIFHTforMember defined zero exposure to fire "<<thefireModel->getTag()<<endln;
    return -1;
  } 
  
  //  FireModel** theFireModel;
  // ID** theFireExpFacesID;
   theFireModel[FireActionIndex] = thefireModel;
   theFireExpFacesID[ FireActionIndex] = theFireExposedfaces;
 
  FireActionIndex++;
  return 0;
}

int
SIFHTforMember::SetAmbExposureCons(const ID&theAmbExposedfaces)
{
  
  if (theAmbExposedfaces==0) {
    theAmbExpFacesID=0;
  } else {
	theAmbExpFacesID = new ID();
    (*theAmbExpFacesID) = theAmbExposedfaces;
  }
  return 0;
}


//////////////------------Applying fire to member------------------///////
int
SIFHTforMember::applyFire(const Vector& fireOrigin, Matrix* crdMat)
{
   //HeatTransfer Model will be initialised here
	//Composite beam should be identified
	int NumSections = 1;
	int FireType = theFireModel[0]->getFireTypeTag();
	const Vector FireOrigin = fireOrigin;
	Vector damageVec = 0;
	if( (FireType==1||FireType==2)&&theMember->getPartialDamage(damageVec)==0){
		//Build 2d model if it is the uniform fire
	Vector SecCentroids(2);
	SecCentroids(0)=0;  SecCentroids(1) =0 ;
	//opserr<<"memberType ---"<< theMember->getMemberTypeTag();
	if(theMember->getMemberTypeTag()==1||theMember->getMemberTypeTag()==2||theMember->getMemberTypeTag()==3||theMember->getMemberTypeTag()==21){
		this->BuildHTModel2D(SecCentroids);
	}
	else if(theMember->getMemberTypeTag()==10||theMember->getMemberTypeTag()==11){
		SecCentroids.resize(1); SecCentroids(0)=0;
		this->BuildHTModel1D(SecCentroids);
	}
    this->RunHTAnalysis(FireDuration,TimeStep);
	}
	
	
	else if(FireType==3||FireType==4||theMember->getPartialDamage(damageVec)>0){
    //Then it can be series 2D section or 3d ;
	NumSections = 4;
	//For considering partial damage!!!!!!!!!!!!!!
	if(theMember->getMemberTypeTag()==3&&theMember->getPartialDamage(damageVec)>0)
		NumSections = 2;

    Vector SecCentroids(3);
	ID joints = theMember->getConnectedJoints();
	SIFJoint* theJoint1=0;
	SIFJoint* theJoint2=0;
	SIFJoint* theJoint3=0;
	SIFJoint* theJoint4=0;
	Vector** jtCrds = new Vector* [4];
	if(theMember->getMemberTypeTag()==1||theMember->getMemberTypeTag()==21||theMember->getMemberTypeTag()==2||theMember->getMemberTypeTag()==3||theMember->getMemberTypeTag()==11){
		theJoint1 = theSIFDomain->getSIFJoint(joints(0));
		theJoint2 = theSIFDomain->getSIFJoint(joints(1));
		jtCrds[0] = new Vector();jtCrds[1] = new Vector();
		*(jtCrds[0]) = theJoint1->getCrds();
		*(jtCrds[1]) = theJoint2->getCrds();
	}
	else if(theMember->getMemberTypeTag()==10){
		theJoint1 = theSIFDomain->getSIFJoint(joints(0));
		theJoint2 = theSIFDomain->getSIFJoint(joints(1));
		theJoint3 = theSIFDomain->getSIFJoint(joints(2));
		theJoint4 = theSIFDomain->getSIFJoint(joints(3));
		jtCrds[0] = new Vector();jtCrds[1] = new Vector();jtCrds[2] = new Vector();jtCrds[3] = new Vector();
		*(jtCrds[0]) = theJoint1->getCrds();
		*(jtCrds[1]) = theJoint2->getCrds();
		*(jtCrds[2]) = theJoint3->getCrds();
		*(jtCrds[3]) = theJoint4->getCrds();
	}
	for(int i=0; i<NumSections;i++){
		if(theMember->getMemberTypeTag()==1){
			//XBeam, x=-z(g), y=y(g), z =x(g)
			SecCentroids(0) =-(*(jtCrds[0]))(2)+(-(*(jtCrds[1]))(2)+(*(jtCrds[0]))(2))*i/(NumSections-1);
			SecCentroids(1) = (*(jtCrds[0]))(1)+((*(jtCrds[1]))(1)-(*(jtCrds[0]))(1))*i/(NumSections-1);
			SecCentroids(2) = (*(jtCrds[0]))(0)+((*(jtCrds[1]))(0)-(*(jtCrds[0]))(0))*i/(NumSections-1); 
			this->BuildHTModel2D(SecCentroids,i+1);
			(*crdMat)(i,0) = SecCentroids(2);
			(*crdMat)(i,1) = SecCentroids(1);
			(*crdMat)(i,2) = -SecCentroids(0);

		}
		else if(theMember->getMemberTypeTag()==21){
			//XBeam, x=-z(g), y=y(g), z =x(g)
			SecCentroids(0) =-(*(jtCrds[0]))(2)+(-(*(jtCrds[1]))(2)+(*(jtCrds[0]))(2))*i/(NumSections-1);
			SecCentroids(1) = (*(jtCrds[0]))(1)+((*(jtCrds[1]))(1)-(*(jtCrds[0]))(1))*i/(NumSections-1);
			SecCentroids(2) = (*(jtCrds[0]))(0)+((*(jtCrds[1]))(0)-(*(jtCrds[0]))(0))*i/(NumSections-1); 
			this->BuildHTModel2D(SecCentroids,i+1);
			(*crdMat)(i,0) = SecCentroids(2);
			(*crdMat)(i,1) = SecCentroids(1);
			(*crdMat)(i,2) = -SecCentroids(0);

		}
		else if(theMember->getMemberTypeTag()==2){
			//ZBeam
			SecCentroids(0) =(*(jtCrds[0]))(0)+((*(jtCrds[1]))(0)-(*(jtCrds[0]))(0))*i/(NumSections-1);  
			SecCentroids(1) = (*(jtCrds[0]))(1)+((*(jtCrds[1]))(1)-(*(jtCrds[0]))(1))*i/(NumSections-1); 
			SecCentroids(2) = (*(jtCrds[0]))(2)+((*(jtCrds[1]))(2)-(*(jtCrds[0]))(2))*i/(NumSections-1); 
			this->BuildHTModel2D(SecCentroids,i+1);
			(*crdMat)(i,0) = SecCentroids(0);
			(*crdMat)(i,1) = SecCentroids(1);
			(*crdMat)(i,2) = SecCentroids(2);
		}
		else if(theMember->getMemberTypeTag()==3){
			//Column
			SecCentroids(0) =-(*(jtCrds[0]))(2)+(-(*(jtCrds[1]))(2)+(*(jtCrds[0]))(2))*i/(NumSections-1);  
			SecCentroids(1) = -(*(jtCrds[0]))(0)+(-(*(jtCrds[1]))(0)+(*(jtCrds[0]))(0))*i/(NumSections-1); 
			SecCentroids(2) = (*(jtCrds[0]))(1)+((*(jtCrds[1]))(1)-(*(jtCrds[0]))(1))*i/(NumSections-1); 

			this->BuildHTModel2D(SecCentroids,i+1);

			(*crdMat)(i,0) = -SecCentroids(1);
			(*crdMat)(i,1) = SecCentroids(2);
			(*crdMat)(i,2) = -SecCentroids(0);
		}
		else if(theMember->getMemberTypeTag()==10|| theMember->getMemberTypeTag() == 11){
			//Slab
			//opserr << FireOrigin << endln;
			double dist =0; double MaxDist =0; int MaxDistJtTag =0;
			for(int k=0;k<4;k++){
				//opserr << *(jtCrds[k]) << endln;
				double dist = ((*(jtCrds[k]))(0)-FireOrigin(0))*((*(jtCrds[k]))(0)-FireOrigin(0));
				dist += ((*(jtCrds[k]))(2)-FireOrigin(2))*((*(jtCrds[k]))(2)-FireOrigin(2));
				dist = sqrt(dist);
				if(dist>MaxDist){
					MaxDist = dist;
					MaxDistJtTag =k;
				}
			}
			
			if(i==0){
				SecCentroids(0) =FireOrigin(2);  
				SecCentroids(1) = FireOrigin(0); 
				SecCentroids(2) = (*(jtCrds[MaxDistJtTag]))(1);
			}
			else if(i==1){
				SecCentroids(0) =(FireOrigin(2)+(*(jtCrds[MaxDistJtTag]))(2))*0.3;  
				SecCentroids(1) = (FireOrigin(0)+(*(jtCrds[MaxDistJtTag]))(0))*0.3;  
				SecCentroids(2) = (*(jtCrds[MaxDistJtTag]))(1);  
			
			}
			else if(i==2){
				SecCentroids(0) =(FireOrigin(2)+(*(jtCrds[MaxDistJtTag]))(2))*0.7;  
				SecCentroids(1) = (FireOrigin(0)+(*(jtCrds[MaxDistJtTag]))(0))*0.7;  
				SecCentroids(2) = (*(jtCrds[MaxDistJtTag]))(1);  
			
			}
			else if(i==3){
				//Slab corner node which is farthest to fire origin
				SecCentroids(0) =(*(jtCrds[MaxDistJtTag]))(2);  
				SecCentroids(1) = (*(jtCrds[MaxDistJtTag]))(0);
				SecCentroids(2) = (*(jtCrds[MaxDistJtTag]))(1);

			}
			(*crdMat)(i,0) = SecCentroids(1);
			(*crdMat)(i,1) = SecCentroids(2);
			(*crdMat)(i,2) = SecCentroids(0);
			this->BuildHTModel1D(SecCentroids,i+1);
		}

#ifdef _DEBUG
		opserr<<"SIFHTforMember:"<<this->getTag()<<" SecCentroids "<< SecCentroids<<endln<<(*(jtCrds[0]))<<endln<<(*(jtCrds[1]))<<endln;
#endif
		
		this->RunHTAnalysis(FireDuration,TimeStep);
	}
  }


  return 0;

}




//set up a 2d heat transfer model
int 
SIFHTforMember::BuildHTModel2D(const Vector& SecCentroids, int SecTag)
{
	Vector OriginLocs = SecCentroids;  // section centriod cordinates
    Vector sectionLoc = 0;
  //check centroid coordinates
  if (OriginLocs.Size()<2) {
    opserr<<"WARNING::SIFHTforMember is expecting to have at least 2 arguments as section cordinates in BuildHTModel2D"<<this->getTag()<<endln;
    return -2;
  }
  else if(OriginLocs.Size()==2)
  {
    sectionLoc = 0;
  }
  else
  {
    sectionLoc =  Vector(1);
    sectionLoc(0)= OriginLocs(2); //Section location
    
  }
	
  ID PCIndexes; // phase change indexes;
	
  int sectionType = theSection->getSectionTypeTag();
  
  
  //-----------------------***Geometrical mesh***------------------------------
  
  //theHTMeshes would store the pointers to SimpleMesh
  //theHTMeshes = new ArrayOfTaggedObjects(10);
  
  
  if (theMember->getTypeTag()==1||theMember->getTypeTag()==3) {
	  //--**********---single beam or column(without slab)-------********-------------
    if (sectionType == 1) {
      //Rectangular section
      Vector SectionPars (2);
      SectionPars = theSection->getSectionPars();
    

      theEntity = new Simple_Block(1,OriginLocs(0),OriginLocs(1), SectionPars(0), SectionPars(1));
      Vector MeshCtrls = Vector(2);
      //Using default settings(20eles) for meshCtrl, This can be refined later
      MeshCtrls(0) = SectionPars(0)/24;    //Breadth of block
      MeshCtrls(1) = SectionPars(1)/24;    //Height of block
	
	  RecLocations = Vector(9);   //to be improved
	  //double Orilocy = OriginLocs(1); 
	  double secd = SectionPars(1);
	  for (int k =0;k<9;k++){  
		 RecLocations(k) = -secd/2+secd/8*k;
	  }


      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls);
    
      PCIndexes.resize(1);
      PCIndexes(0)= 1;
    }
    else if(sectionType == 2)
    {
      //HeatTransfer model for I beam;
      Vector SectionPars (4);
      SectionPars = theSection->getSectionPars();
	  double d = SectionPars(0);
	  double tw = SectionPars(1);
		double bf = SectionPars(2);
		double tf = SectionPars(3);
		double wd = d - 2*tf;				//web depth;
    
      //Simple_Isection(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw);
      theEntity = new Simple_Isection(1,OriginLocs(0),OriginLocs(1),bf, tf ,tw, wd);
      
	   RecLocations = Vector(4);  
	   RecLocations(0)= -wd/2; RecLocations(1)= wd/2;
	    RecLocations(2)= -bf/2; RecLocations(3)= bf/2;

      
      Vector MeshCtrls = Vector(4);
      MeshCtrls(0)=(bf-tw)/20;
      MeshCtrls(1)=tf/4;
      MeshCtrls(2)= tw/4;
      MeshCtrls(3)= wd/20;
    
      
      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls);
    
      PCIndexes.resize(1);
      PCIndexes(0)= 0;
    }

	 else if(sectionType == 22)
    {
	Vector damageVec=0;
	if(theMember->getMemberTypeTag()==3&&theMember->getPartialDamage(damageVec)>0&&SecTag==2)
	{
		//HeatTransfer model for I beam;
      Vector SectionPars (4);
      SectionPars = theSection->getSectionPars();
	  double d = SectionPars(0);
	  double tw = SectionPars(1);
		double bf = SectionPars(2);
		double tf = SectionPars(3);
		double wd = d - 2*tf;				//web depth;
    
      //Simple_Isection(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw);
      theEntity = new Simple_Isection(1,OriginLocs(0),OriginLocs(1),bf, tf ,tw, wd);
      
	   RecLocations = Vector(4);  
	   RecLocations(0)= -wd/2; RecLocations(1)= wd/2;
	    RecLocations(2)= -bf/2; RecLocations(3)= bf/2;

      
      Vector MeshCtrls = Vector(4);
      MeshCtrls(0)=(bf-tw)/20;
      MeshCtrls(1)=tf/4;
      MeshCtrls(2)= tw/4;
      MeshCtrls(3)= wd/20;
    
      
      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls);
    
      PCIndexes.resize(1);
      PCIndexes(0)= 0;
	}
	else {
      //HeatTransfer model for I section protected;
      Vector SectionPars (5);
      SectionPars = theSection->getSectionPars();
	  double d = SectionPars(0);
	  double tw = SectionPars(1);
		double bf = SectionPars(2);
		double tf = SectionPars(3);
		double wd = d - 2*tf;				//web depth;
		double coat = SectionPars(4);
    
      //Simple_Isection(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw);
      theEntity = new Simple_IsecProtected(1,OriginLocs(0),OriginLocs(1),bf, tf ,tw, wd,coat);
      
	   RecLocations = Vector(4);  
	   RecLocations(0)= -wd/2; RecLocations(1)= wd/2;
	   RecLocations(2)= -bf/2; RecLocations(3)= bf/2;
      
      Vector MeshCtrls = Vector(5);
      MeshCtrls(0)=(bf-tw-2*coat)/20;
      MeshCtrls(1)=tf/4;
      MeshCtrls(2)= tw/4;
      MeshCtrls(3)= (wd-2*coat)/20;
	  MeshCtrls(4)= coat/5;

      //  Heat transfer Material for SFRM coating
	  //SimpleMaterial(int tag, double my_rho, double my_cp, double kc)
	  HeatTransferMaterial* theHTMaterialSFRM = new SimpleMaterial(10000,350,1100,0.05);
      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls,theHTMaterialSFRM);
    
      PCIndexes.resize(2);
      PCIndexes(0)= 0;
	  PCIndexes(1)= 0;

	 }
    }

    else
    {
      opserr<<"WARNING::SIFHTforMember recieves an invalid section type "<<sectionType<<endln;
    }
	//end of section type
  }
  //--Composite beam(with slab)-------
  else if(theMember->getTypeTag()==2){
    
     if (sectionType == 1) {
    //Rectangular section
    // here it is going to be reinforced concrete beam incluidng slab(under development)
    
    
    //Do the mesh
    }
    else if(sectionType == 2) {
      //HeatTransfer model for I beam;
      //Rectangular section
      Vector SectionPars (4);
      SectionPars = theSection->getSectionPars();

	  double d = SectionPars(0);
	  double tw = SectionPars(1);
		double bf = SectionPars(2);
		double tf = SectionPars(3);
		double wd = d - 2*tf;				//web depth;

      double slabThickness = (theSection1->getSectionPars())(0);
	  SIFMaterial* theSlabMat = theSection1->getSIFMaterialPtr();
	  HeatTransferMaterial* theSlabHTMat = theSlabMat->getHeatTransferMaterial();
      double slabBreadth = SectionPars(0)+0.6; //Add 600mm to beam breadth
      if (slabThickness <=0){
        opserr<<"WARNING::SIFHTforMember can't define slab without thickness"<<endln;
      }

	  RecLocations = Vector(4);  
	   RecLocations(0)= -wd/2; RecLocations(1)= wd/2;
	    RecLocations(2)= -bf/2; RecLocations(3)= bf/2;
                                       
      //Simple_Composite2D(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw,double Slab_W, double Slab_H);
                                       
      theEntity = new Simple_Composite2D(1, OriginLocs(0),OriginLocs(1),bf, tf,tw, wd,slabBreadth, slabThickness);

      Vector MeshCtrls = Vector(6);
      
      MeshCtrls(0)=(bf-tw)/20;
      MeshCtrls(1)=tf/4;
      MeshCtrls(2)= tw/4;
      MeshCtrls(3)= wd/20;
      MeshCtrls(4)= (slabBreadth-slabThickness)/20;
      MeshCtrls(5)= slabThickness/16;
      
      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls,theSlabHTMat);
    
    }
    else{
      opserr<<"WARNING::SIFHTforMember recieves an invalid section type "<<sectionType<<endln;
    }
    //end of section type
    
  }
  //end of member type 
  
  
  //Once the meshes are reday, start to mesh
  if (theMesh!=0) {
    theMesh->GeneratingNodes(sectionLoc);
    theMesh->GeneratingEles(PCIndexes);
  }
  
  if (theMesh1!=0) {
    theMesh1->GeneratingNodes(sectionLoc);
    theMesh1->GeneratingEles(PCIndexes);
  }

  
  //--------------------***Definition of Fire exposure***------------------------//
  
  //set the intial temperature of the structural member;
  theHTDomain->setInitial(293.15);
  
  
  //Creating a boudary operator to apply boundary conditions
  Simple_Boundary* Boundary1 = new Simple_Boundary(1, theHTDomain);
  
  
  //@@Ambient boundary conditions
  BoundaryPattern* thePattern = new BoundaryPattern(1);
  theHTDomain->addBoundaryPattern(thePattern);
  
  //This can be improved later(under development)
  double AmbientHFConstants[6] = {4.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
  Vector AmbientHeatFluxConstants(AmbientHFConstants,6);

   int NumAmbExpFaces;

  if(theAmbExpFacesID==0)
	NumAmbExpFaces =0;
  else 
	NumAmbExpFaces = theAmbExpFacesID->Size();
  if(NumAmbExpFaces>0){
    //if there is ambient exposed faces existed
   for(int i=0; i<NumAmbExpFaces; i++ ){
		int FaceID = (*theAmbExpFacesID)(i);
		ID AmbExpEles = ID();int facetag;
		theMesh->SelectingElesbyFace(AmbExpEles,FaceID,facetag);
      
		Boundary1->GeneratingHeatFluxBC(AmbExpEles,facetag,1,1,AmbientHeatFluxConstants);
		Boundary1->GeneratingHeatFluxBC(AmbExpEles,facetag,2,1,AmbientHeatFluxConstants);
      
    }
    
  }

//-----------------------***Definition of Fire exposure***--------------------------//
  if(FireActionIndex>0){
    //if there is ambient exposed faces existed
	  for(int k=0; k<FireActionIndex;k++){
		//@@Fire exposure boundary conditions
		FireModel* TheFireModel = theFireModel[k];
		int FireModelType = TheFireModel->getFireTypeTag();
		

		FireImposedPattern* FirePattern = new FireImposedPattern(k+2);
		FirePattern->setFireModel(TheFireModel);
		theHTDomain->addBoundaryPattern(FirePattern);
		
		ID* theFireExpoID = theFireExpFacesID[k];
		int NumFireExpFaces = theFireExpoID->Size();

		for(int i=0; i<NumFireExpFaces; i++ ){
		int FaceID = (*theFireExpoID)(i);
		ID FireExpEles = ID();int facetag;
		theMesh->SelectingElesbyFace(FireExpEles,FaceID,facetag);

		if(FireModelType==1||FireModelType==2){
			//This can be improved later(under development)
			//This should be checked for different nominal fire
			double FireHFConstants[6] = {25.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
			Vector FireHeatFluxConstants(FireHFConstants,6);

			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,1,k+2,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,2,k+2,FireHeatFluxConstants);
		}
		//Uniform fire exposrure
		if(FireModelType==3||FireModelType==5||FireModelType==6||FireModelType==7){
			//This can be improved later(under development)
			double FireHFConstants[6] = {35.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
			Vector FireHeatFluxConstants(FireHFConstants,6);
			
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,3,k+2,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,1,1,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,2,1,FireHeatFluxConstants);
		}
		//localised fire
		}
		//for dffernt faces 
    }
  //for fire action index
	}
  
    OPS_Stream *theOutputStream = 0;
	TCL_Char *fileName = 0;
	ID* theNodes =0;
	HTRecorder* theHTRecorder =0;
	ID* MonitorNodes;

	const char* theHTdir = theSIFDomain->getHTdir();
  
  //----------*****Adding recorder for internal transition or storing in a file***---------------//
  //First determine what kind of recorder is needed
  
	  //Rect
	  if(theSection->getSectionTypeTag()==1){
		MonitorNodes = new ID();
		//GetNodesForRecorder(ID& NodesRange, const Vector& MonitorLocX, double PlaneLoc)
		Vector locx = Vector(1);
		locx(0)=0.0;
		theMesh->GetNodesForRecorder(*MonitorNodes,locx);
	  }
	  //Isec
	  if(theSection->getSectionTypeTag()==2||theSection->getSectionTypeTag()==22){
		MonitorNodes = new ID();
		//GetNodesForRecorder(ID& NodesRange, const Vector& MonitorLocX, double PlaneLoc)
		theMesh->GetNodesForRecorder(*MonitorNodes,3);
	  }
    
		//HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes,
		//HeatTransferDomain &theDom,PathTimeSeriesThermal* thePathSeries)
	  PathTimeSeriesThermal* thePathtimeSeries;
	  if(SecTag==0)
		 thePathtimeSeries  = thePathSeries[SecTag];
	  else if(SecTag>0)
		thePathtimeSeries  = thePathSeries[SecTag-1];
	  else 
		  opserr<<"SIFHTforMember::receiving invalid SecTag "<<SecTag<<endln;

		HTNodeRecorder* theRecorderBeam;
		
		if(theHTdir==0) 
			theRecorderBeam= new HTNodeRecorder(theMember->getTag()+SecTag, MonitorNodes,*theHTDomain, thePathtimeSeries);
		else {
			std::string thePrefix = theHTdir; 

			std::string recorder;
			int SecBeamTag;
			if(theMember->getMemberTypeTag()==1)
				recorder = "/HT_SIFXBeam";
			else if(theMember->getMemberTypeTag()==2)
				 recorder = "/HT_SIFYBeam";
			else if(theMember->getMemberTypeTag()==3)
				 recorder = "/HT_SIFColumn";
			else if(theMember->getMemberTypeTag()==21)
				 recorder = "/HT_SIFSecXBeam";


			std::string recorderSuffix = ".out";
			int Rtag = theMember->getTag();

			if(theMember->getMemberTypeTag()==21){
				ID theCompInfo = theMember->getMemberInfo();
				Rtag = theCompInfo(0)*1000+theCompInfo(1)*100+theCompInfo(2);
			}

			std::string section = "-Sec";

			std::stringstream f;
			if(SecTag==0)
				f <<thePrefix<< recorder << Rtag << recorderSuffix;
			else{
				if(theMember->getMemberTypeTag()==21){
				ID theCompInfo = theMember->getMemberInfo();
				SecTag = theCompInfo(3)*100+ SecTag;
				}
				f <<thePrefix<< recorder << Rtag <<section<<SecTag<< recorderSuffix;
			}

			std::string fileNameStr;
			fileNameStr = f.str();
			fileName = const_cast<char*>(fileNameStr.c_str());

			const char *pwd = theSIFDomain->getInterpPWD();
			simulationInfo.addOutputFile(fileName, pwd);
			theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 0 );

			theRecorderBeam= new HTNodeRecorder(theMember->getTag()+SecTag, MonitorNodes,*theHTDomain,*theOutputStream, thePathtimeSeries);
		}
		theHTDomain->addRecorder(*theRecorderBeam);
	
  
  return 0;
}


//set up a 1d heat transfer model
int 
SIFHTforMember::BuildHTModel1D(const Vector& SecCentroids, int SecTag)
{
	Vector OriginLocs = SecCentroids;  // section centriod cordinates
    Vector sectionLoc = 0;
  //check centroid coordinates
  if (OriginLocs.Size()==1) {
	sectionLoc = 0;
	OriginLocs.resize(3);
	OriginLocs.Zero();
  }
  else if(OriginLocs.Size()==3)
  {
    sectionLoc =  Vector(2);
    sectionLoc(0)= OriginLocs(0); 
	sectionLoc(1)= OriginLocs(1);//Section location
  }
  else
  {
	  opserr<<"WARNING:: SIFHTforMember-Unknown section location for BuildHTModel1D"<<endln;
  }
	
  ID PCIndexes; // phase change indexes;
	
  int sectionType = theSection->getSectionTypeTag();
  
  
  //-----------------------***Geometrical mesh***------------------------------
  
  //theHTMeshes would store the pointers to SimpleMesh
  //theHTMeshes = new ArrayOfTaggedObjects(10);
  
 
	//--**********---single slab-------********-------------
    if (sectionType == 10) {
      //Rectangular section
      Vector SectionPars (2);
      SectionPars = theSection->getSectionPars();
    
      theEntity = new Simple_Line(1,OriginLocs(2), SectionPars(0));
      Vector MeshCtrls = Vector(2);
      //Using default settings(20eles) for meshCtrl, This can be refined later
      MeshCtrls(0) = SectionPars(0)/24;    //thickness of slab
   
	  RecLocations = Vector(9);   //to be improved
	  //double Orilocy = OriginLocs(1); 
	  double secd = SectionPars(0);
	  for (int k =0;k<9;k++){  
		 RecLocations(k) = -secd/2+secd/8*k;
	  }

      theMesh = new Simple_Mesh(1, theEntity, theHTDomain,theHTMaterial,MeshCtrls);
    
      PCIndexes.resize(1);
      PCIndexes(0)= 1;
	}
  
  //Once the meshes are reday, start to mesh
  if (theMesh!=0) {
    theMesh->GeneratingNodes(sectionLoc);
    theMesh->GeneratingEles(PCIndexes);
  }
  
  //--------------------***Definition of Fire exposure***------------------------//
  
  //set the intial temperature of the structural member;
  theHTDomain->setInitial(293.15);
  
  
  //Creating a boudary operator to apply boundary conditions
  Simple_Boundary* Boundary1 = new Simple_Boundary(1, theHTDomain);
  
  
  //@@Ambient boundary conditions
  BoundaryPattern* thePattern = new BoundaryPattern(1);
  theHTDomain->addBoundaryPattern(thePattern);
  
  //This can be improved later(under development)
  double AmbientHFConstants[6] = {4.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
  Vector AmbientHeatFluxConstants(AmbientHFConstants,6);

   int NumAmbExpFaces;

  if(theAmbExpFacesID==0)
	NumAmbExpFaces =0;
  else 
	NumAmbExpFaces = theAmbExpFacesID->Size();
  if(NumAmbExpFaces>0){
    //if there is ambient exposed faces existed
   for(int i=0; i<NumAmbExpFaces; i++ ){
		int FaceID = (*theAmbExpFacesID)(i);
		ID AmbExpEles = ID();int facetag;
		theMesh->SelectingElesbyFace(AmbExpEles,FaceID,facetag);
      
		Boundary1->GeneratingHeatFluxBC(AmbExpEles,facetag,1,1,AmbientHeatFluxConstants);
		Boundary1->GeneratingHeatFluxBC(AmbExpEles,facetag,2,1,AmbientHeatFluxConstants);
      
    }
    
  }

//-----------------------***Definition of Fire exposure***--------------------------//
  if(FireActionIndex>0){
    //if there is ambient exposed faces existed
	  for(int k=0; k<FireActionIndex;k++){
		//@@Fire exposure boundary conditions
		FireModel* TheFireModel = theFireModel[k];
		int FireModelType = TheFireModel->getFireTypeTag();
		
		FireImposedPattern* FirePattern = new FireImposedPattern(k+2);
		FirePattern->setFireModel(TheFireModel);
		theHTDomain->addBoundaryPattern(FirePattern);
		
		ID* theFireExpoID = theFireExpFacesID[k];
		int NumFireExpFaces = theFireExpoID->Size();

		for(int i=0; i<NumFireExpFaces; i++ ){
		int FaceID = (*theFireExpoID)(i);
		ID FireExpEles = ID();int facetag;
		theMesh->SelectingElesbyFace(FireExpEles,FaceID,facetag);

		if(FireModelType==1||FireModelType==2){
			//This can be improved later(under development)
			//This should be checked for different nominal fire
			double FireHFConstants[6] = {25.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
			Vector FireHeatFluxConstants(FireHFConstants,6);

			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,1,k+2,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,2,k+2,FireHeatFluxConstants);
		}
		else if ( FireModelType == 21) {
			//This can be improved later(under development)
			//This should be checked for different nominal fire
			double FireHFConstants[6] = { 35.0, 293.15, 0.7, 5.67e-8,0.7, 418.738 };
			Vector FireHeatFluxConstants(FireHFConstants, 6);

			Boundary1->GeneratingHeatFluxBC(FireExpEles, facetag, 1, k + 2, FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles, facetag, 2, k + 2, FireHeatFluxConstants);
		}
		//Uniform fire exposrure
		else if(FireModelType==3||FireModelType==6||FireModelType==7){
			//This can be improved later(under development)
			double FireHFConstants[6] = {35.0, 293.15, 0.7, 5.67e-8,0.7, 418.738};
			Vector FireHeatFluxConstants(FireHFConstants,6);
			
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,3,k+2,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,1,1,FireHeatFluxConstants);
			Boundary1->GeneratingHeatFluxBC(FireExpEles,facetag,2,1,FireHeatFluxConstants);
		}
		//localised fire
		}
		//for dffernt faces 
    }
  //for fire action index
	}
  
    OPS_Stream *theOutputStream = 0;
	TCL_Char *fileName = 0;
	ID* theNodes =0;
	HTRecorder* theHTRecorder =0;
	ID* MonitorNodes;

	const char* theHTdir = theSIFDomain->getHTdir();
  
  //----------*****Adding recorder for internal transition or storing in a file***---------------//
  //First determine what kind of recorder is needed
  if(theMember->getMemberTypeTag()==10||theMember->getMemberTypeTag()==11)
  {
	  
		ID* MoniNodes = new ID();
		MonitorNodes = new ID(9);
		//GetNodesForRecorder(ID& NodesRange, const Vector& MonitorLocX, double PlaneLoc)
		theMesh->GetNodesForRecorder(*MoniNodes); //Shell local coordinate system is reversed in y axis
		if(MoniNodes->Size()!=9){
			opserr<<"WARNIGN::SIFHTforMember recives invalid nodes for defining recorder in BuildHTModel 1D"<<endln;
		}
		
		for(int i=0;i<9;i++){
		  if(theMember->getMemberTypeTag()==10)
			(*MonitorNodes)(i) = (*MoniNodes)(8-i);
		  if(theMember->getMemberTypeTag()==11)
			(*MonitorNodes)(i) = (*MoniNodes)(i);
		}
		delete MoniNodes;
        //opserr<<*MonitorNodes<<endln;
		//HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes,
		//HeatTransferDomain &theDom,PathTimeSeriesThermal* thePathSeries)
	  PathTimeSeriesThermal* thePathtimeSeries;
	  if(SecTag==0)
		 thePathtimeSeries  = thePathSeries[SecTag];
	  else if(SecTag>0)
		thePathtimeSeries  = thePathSeries[SecTag-1];
	  else 
		  opserr<<"SIFHTforMember::receiving invalid SecTag "<<SecTag<<endln;

		HTNodeRecorder* theRecorderBeam;
		
		if(theHTdir==0) 
			theRecorderBeam= new HTNodeRecorder(theMember->getTag()+SecTag, MonitorNodes,*theHTDomain, thePathtimeSeries);
		else {
			std::string thePrefix = theHTdir; 

			std::string recorder;
			if(theMember->getMemberTypeTag()==10||theMember->getMemberTypeTag()==11)
				recorder = "/HT_SIFSlab";

			std::string recorderSuffix = ".out";
			int Rtag = theMember->getTag();

			std::string section = "-Sec";

			std::stringstream f;
			if(SecTag==0)
				f <<thePrefix<< recorder << Rtag << recorderSuffix;
			else
				f <<thePrefix<< recorder << Rtag <<section<<SecTag<< recorderSuffix;

			std::string fileNameStr;
			fileNameStr = f.str();
			fileName = const_cast<char*>(fileNameStr.c_str());

			const char *pwd = theSIFDomain->getInterpPWD();
			simulationInfo.addOutputFile(fileName, pwd);
			theOutputStream = new DataFileStream(fileName, OVERWRITE, 2, 0 );

			theRecorderBeam= new HTNodeRecorder(theMember->getTag()+SecTag, MonitorNodes,*theHTDomain,*theOutputStream, thePathtimeSeries);
		}
		theHTDomain->addRecorder(*theRecorderBeam);
	
  }
  
  return 0;
}


int
SIFHTforMember::RunHTAnalysis(double fireDuration, double timeStep)
{
  
  
  if (theHTModel==0)
	 theHTModel = new HT_AnalysisModel();

  if (theHTTest ==0) 
	theHTTest  = new CTestNormTempIncr(1e-3, 200,0);

  if (theHTSolnAlgo ==0) 
	theHTSolnAlgo = new NewtonMethod(*theHTTest);

  if (theHTIntegrator  ==0) 
	 theHTIntegrator   = new BackwardDifference();

  if (theHTHandler==0)
	theHTHandler  = new PenaltyBC_Handler(1.0e10,1.0e10);

  if (theHTRCM==0)
    theHTRCM = new RCM();

  if (theHTNumberer==0) 
     theHTNumberer = new HT_DOF_Numberer(*theHTRCM);

  if (theHTSolver==0)
    theHTSolver = new BandGenLinLapackSolver();

  if (theHTSOE==0) 
    theHTSOE  = new BandGenLinSOE(*theHTSolver);

  if (theHTAnalysis!= 0) {
	delete theHTAnalysis;
  }
  
  theHTAnalysis = new HT_TransientAnalysis(*theHTDomain,
							     *theHTHandler,
							     *theHTNumberer,
							     *theHTModel,
							     *theHTSolnAlgo,
							     *theHTSOE,
							     *theHTIntegrator,
							     theHTTest);

  
  // perform the analysis & print out the results for the domain
  //int numSteps = 60;
  
  
  int numSteps = fireDuration/timeStep;
  double time = 0;
  theHTAnalysis->analyze(numSteps, timeStep,time);
  
 
  //if finished
  theHTDomain->clearAll();
  return 0;
}


int
SIFHTforMember::SetPathTimeSeries(PathTimeSeriesThermal** thePathTimeSeries)
{
		if (thePathTimeSeries!=0) {
      thePathSeries = thePathTimeSeries;
    }
  else
  {
    opserr<<"SIFHTforMember::SetPathTimeSeries can't allocate an empty PathTimeSeries"<<endln;
	return -2;
  }
  return 0;
}



const Vector&
SIFHTforMember::getRecLocations()
{
  return RecLocations;
}