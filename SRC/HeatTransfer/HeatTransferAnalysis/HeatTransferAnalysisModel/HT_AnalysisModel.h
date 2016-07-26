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
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

// This class is a modified version of AnalysisModel 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// AnalysisModel
// Written: fmk 
// Created: 9/96
// Revision: A


#ifndef HT_AnalysisModel_h
#define HT_AnalysisModel_h

//#include <MovableObject.h>

class TaggedObjectStorage;
class HeatTransferDomain;
class HT_FE_EleIter;
class HT_DOF_GrpIter;
class Graph;
class HT_FE_Element;
class HT_DOF_Group;
class Vector;
class TemperatureBCHandler;

class HT_AnalysisModel
{
  public:
    HT_AnalysisModel();    
    virtual ~HT_AnalysisModel();    

    // methods to populate/depopulate the AnalysisModel
    virtual bool addHT_FE_Element(HT_FE_Element* theFE_Ele);
    virtual bool addHT_DOF_Group(HT_DOF_Group* theDOF_Grp);
    virtual void clearAll(void);
    
    // methods to access the FE_Elements and DOF_Groups and their numbers
    virtual int getNumDOF_Groups(void) const;		
    virtual HT_DOF_Group* getDOF_GroupPtr(int tag);	
    virtual HT_FE_EleIter& getFEs();
    virtual HT_DOF_GrpIter& getDOFs();

    // method to access the connectivity for SysOfEqn to size itself
    virtual void setNumEqn(int) ;	
    virtual int getNumEqn(void) const ; 
    virtual Graph& getDOFGraph(void);
    virtual Graph& getDOFGroupGraph(void);
    
    // methods to update the response quantities at the DOF_Groups,
    // which in turn set the new nodal trial response quantities.
    virtual void setResponse(const Vector& T, 
			                 const Vector& Tdot);
    virtual void setTemp(const Vector& T);    
    virtual void setTdot(const Vector& Tdot);                  

    virtual void incrTemp(const Vector& T);    
    virtual void incrTdot(const Vector& Tdot);                  
    
    // methods which trigger operations in the Domain
    virtual void setLinks(HeatTransferDomain& theDomain, TemperatureBCHandler& theHandler);
	
    // virtual void   applyLoadDomain(double newTime);
    virtual int    updateDomain(void);
    virtual int    updateDomain(double newTime, double dT);
    virtual int    analysisStep(double dT =0.0);

    virtual int    commitDomain(void);
    virtual int    revertDomainToLastCommit(void);
    virtual double getCurrentDomainTime(void);
    virtual void   setCurrentDomainTime(double newTime);      

    HeatTransferDomain* getDomainPtr(void) const;

  protected:

    
  private:
    HeatTransferDomain* myDomain;
    TemperatureBCHandler* myHandler;

    Graph* myDOFGraph;
    Graph* myGroupGraph;    
    
    int numFE_Ele;             // number of FE_Elements objects added
    int numDOF_Grp;            // number of DOF_Group objects added
    int numEqn;                // numEqn set by the ConstraintHandler typically

    TaggedObjectStorage* theFEs;
    TaggedObjectStorage* theDOFs;
    
    HT_FE_EleIter* theFEiter;     
    HT_DOF_GrpIter* theDOFiter;    
};

#endif
