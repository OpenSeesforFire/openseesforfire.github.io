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
** Developed by:                                                     **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $

                                                                        
                                                                        
#ifndef SIFBuilderDomain_h
#define SIFBuilderDomain_h

#include <OPS_Stream.h>
#include <SIFMaterial.h>
#include <SIFSection.h>
#include <SIFJoint.h> 
//include sifmember info
#include <SIFMember.h>
#include <SIFXBeamSec.h>
#include <SIFXBeam.h>
#include <SIFYBeam.h>
#include <SIFColumn.h>
#include <SIFXWall.h>
#include <SIFYWall.h>
#include <SIFSlab.h>
#include <SIFColumn.h>

#include <MapOfTaggedObjectsIter.h>
#include <SIFfireActionIter.h>
#include <SIFJointIter.h>
#include <SIFCompartmentIter.h>
#include <SIFSecXBeamIter.h>
#include <SIFXBeamIter.h>
#include <SIFYBeamIter.h>
#include <SIFColumnIter.h>
#include <SIFSlabIter.h>

#include <SIFCompartment.h>
#include <FireModel.h>
#include <SIFfireAction.h>

#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Domain.h>
class OPS_Stream;
class TaggedObjectStorage;
class SIFMaterial;
class SIFSection;
class SIFJoint;
class SIFCompartment;

class SIFXBeamSec;
class SIFXBeam;
class SIFYBeam;
class SIFXWall;
class SIFYWall;
class SIFColumn;
class SIFSlab;

class SIFJointIter;
class SIFCompartmentIter;
class SIFXBeamIter;
class SIFYBeamIter;
class SIFColumnIter;
class SIFSlabIter;
class SIFSecXBeamIter;

class SIFfireActionIter;
//class SIFfireAction;

class SIFBuilderDomain 
{
  public:
    SIFBuilderDomain( );
	virtual ~SIFBuilderDomain();
	
	int SetStructureDomain(Domain* theStructDomain);
	Domain* getStructureDomain();

	int addSIFMaterial(SIFMaterial *theSIFMaterial);
	SIFMaterial* getSIFMaterial(int tag);

	int addSIFSection(SIFSection *theSIFSection);
	SIFSection* getSIFSection(int tag); 

//---------------------------------

	int addSIFJoint(SIFJoint *theJoint);
	SIFJoint* getSIFJoint(int tag);
	SIFJointIter &getSIFJoints();

	int addSIFXBeam(SIFXBeam *theXBeam);
	SIFXBeam* getSIFXBeam(int tag);
	SIFXBeamIter &getSIFXBeams();

	int addSIFXBeamSec(SIFXBeamSec *theXBeamSec);
	SIFXBeamSec* getSIFXBeamSec(int tag);
	SIFSecXBeamIter &getSIFSecXBeams();

	int addSIFYBeam(SIFYBeam *theYBeam);
	SIFYBeam* getSIFYBeam(int tag);
	SIFYBeamIter &getSIFYBeams();

	int addSIFColumn(SIFColumn *theSIFColumn);
	SIFColumn* getSIFColumn(int tag);
	SIFColumnIter &getSIFColumns();

	int addSIFSlab(SIFSlab *theSlab);
	SIFSlab* getSIFSlab(int tag);
	SIFSlabIter &getSIFSlabs();

	int addSIFCompartment(SIFCompartment *theCompartment);
	SIFCompartment* getSIFCompartment(int tag);
	SIFCompartmentIter &getSIFCompartments();
	SIFCompartment* getSIFCompartment(const ID& compInfo);

//------------------------------------------

	int addSIFXWall(SIFXWall *theXWall);
	SIFXWall* getSIFXWall(int tag);
	//const ID& getAdjCompartments(int theCompID) =0;

	int addSIFfireAction(SIFfireAction *theSIFfireAction);
	SIFfireActionIter &getSIFfireActions();
	SIFfireAction* getSIFfireAction (int tag);

	int applySelfWeight(int thePatternTag, double dt);
	int applyMiscLoad(int thePatternTag, double dt);
	int applyFireAction(int thePatternTag, double timeStep,double fireDuration);
	
	int defineBC();

	int GenStructuralModel(const ID& MeshCtrlPars,bool isElasticModel,bool isDynamicModel, bool isPinnedModel);
  
    int MeshSIFMember(SIFMember* theSIFMember, int NumEles);
	int MeshSIFSlab(SIFMember* theSIFMember, int NumEleX,int NumEleZ);

	CrdTransf* getCrdTransf(int MemberType, int TransfType);

	void setInterpPWD(const char* pwd);
	const char* getInterpPWD();

	void setHTdir(const char* HTdir);
	const char* getHTdir();

	int setSIFBuilderInfo(const ID& theBuilderInfo);
	const ID& getSIFBuilderInfo();

	void clearAll(void);

	int getNodalLoadTag();
	void incrNodalLoadTag();

	int getEleLoadTag();
	void incrEleLoadTag();  
  
  private:
	//TaggedObjectStorage* theSIFJoints;

	TaggedObjectStorage* theSIFMaterials;
	TaggedObjectStorage* theSIFSections;
  
	TaggedObjectStorage* theSIFJoints;
	SIFJointIter* theSIFJoint_Iter;

	TaggedObjectStorage* theSIFCompartments;
	SIFCompartmentIter* theSIFCompartment_Iter;
	
	TaggedObjectStorage* theSIFSecXBeams;
	SIFSecXBeamIter* theSIFSecXBeam_Iter;

	TaggedObjectStorage* theSIFXBeams;
	SIFXBeamIter* theSIFXBeam_Iter;

	TaggedObjectStorage* theSIFYBeams;
	SIFYBeamIter* theSIFYBeam_Iter;
    TaggedObjectStorage* theSIFColumns;
	SIFColumnIter* theSIFColumn_Iter;
    TaggedObjectStorage* theSIFSlabs;
	SIFSlabIter* theSIFSlab_Iter;
  
	TaggedObjectStorage* theSIFXWalls;
    TaggedObjectStorage* theSIFYWalls;
  
	TaggedObjectStorage* theSIFfireActions;
	SIFfireActionIter* theSIFfireAction_Iter;

	bool isElastic;
	bool isPinned;
	int LoadApplied;
	int theNodeTag;
	int theEleTag;
	const char *InterpPWD;
	const char *theHTdir;
	ID SIFBuilderInfo;
	double offset;
	Vector fireFloorTag;

	Domain* theDomain ;                 //static pointer to Domain;
int nDoF ;                            //static tag for number of degree of freedom
 int theNodalLoadTag ;           // static tag for NodalLoadTag;
 int theEleLoadTag ;           // static tag for NodalLoadTag;
};


#endif

