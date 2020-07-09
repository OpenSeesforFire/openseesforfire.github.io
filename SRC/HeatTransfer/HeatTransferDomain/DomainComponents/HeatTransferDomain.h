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

// This class is a modified version of Domain 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// Domain
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A

#ifndef HeatTransferDomain_h
#define HeatTransferDomain_h

#include <OPS_Stream.h>
#include <Vector.h>
#include <ID.h>
#include <vector>

class HeatTransferElement;
class HeatTransferNode;
class TemperatureBC;
class MP_TemperatureBC;
class HeatFluxBC;
class BoundaryPattern;
class HT_ElementIter;
class HT_NodeIter;
class TemperatureBCIter;
class BoundaryPatternIter;
class TempBCIter;
class AllTempBCIter;
class MP_TemperatureBCIter;
class MP_TempBCIter;
class Graph;
class TaggedObjectStorage;
class HTRecorder;

class HeatTransferDomain
{
  public:
    HeatTransferDomain(); 
    virtual ~HeatTransferDomain();    

    // methods to populate a domain
    virtual  int addElement(HeatTransferElement* );
    virtual  int addNode(HeatTransferNode* );
    virtual  bool addTemperatureBC(TemperatureBC* );
    virtual  bool addBoundaryPattern(BoundaryPattern* ); 
	virtual  bool addMP_TemperatureBC(MP_TemperatureBC* );  // Added by Liming for MP_TempBC
    
    // methods to add components to a LoadPattern object
    virtual  bool addTemperatureBC(TemperatureBC* , int PatternTag); 
    virtual  bool addHeatFluxBC(HeatFluxBC* , int PatternTag);

    
    // methods to remove the components 
    virtual void clearAll(void);	
    virtual HeatTransferElement* removeElement(int tag);
    virtual HeatTransferNode* removeNode(int tag);    
    virtual TemperatureBC* removeTemperatureBC(int tag);
	virtual MP_TemperatureBC* removeMP_TemperatureBC(int tag);

    // virtual int removeTemperatureBC(int nodeTag, int dof, int PatternTag);  // not really useful at the moment
    virtual BoundaryPattern* removeBoundaryPattern(int tag);
    virtual HeatFluxBC* removeHeatFluxBC(int tag, int Pattern);
    virtual TemperatureBC* removeTemperatureBC(int tag, int Pattern);
    
    // methods to access the components of a domain
    virtual  HT_ElementIter& getElements();
    virtual  HT_NodeIter& getNodes();
    virtual  TemperatureBCIter& getTemperatureBCs();
	virtual  MP_TemperatureBCIter& getMPTemperatureBCs();
    virtual  BoundaryPatternIter& getBoundaryPatterns();
    virtual  TemperatureBCIter& getAllTempBCs();   //getDomainAndLoadPatternSPs()
    

    virtual  HeatTransferElement* getElement(int tag);
    virtual  HeatTransferNode* getNode(int tag);
    virtual  TemperatureBC* getTemperatureBC(int tag); 
	virtual  MP_TemperatureBC* getMPTemperatureBC(int tag);
    virtual  BoundaryPattern* getBoundaryPattern(int tag);       

    // methods to query the state of the domain
    virtual double getCurrentTime(void) const;
    virtual int getCommitTag(void) const;    	
    virtual int getNumElements(void) const;
    virtual int getNumNodes(void) const;
    virtual int getNumTempBCs(void) const;
	virtual int getNumMPTempBCs(void) const; 
    virtual int getNumBoundaryPatterns(void) const;    
	virtual const Vector& getPhysicalBounds(void); 
    
    // methods to update the domain
    virtual void setCommitTag(int newTag);    	
    virtual void setCurrentTime(double newTime);    
    virtual void setCommittedTime(double newTime);        
    virtual void applyBCs(double Time);
    virtual void setFluxConstant(void);
	virtual void unsetFluxConstant(void);
    virtual int setInitial(double vaule, int dofTag = 0);
	virtual int setInitial(int nodeTag, double value, int dofTag = 0);
    virtual int SelectingNodes(ID& NodesRange, int crdTag, double MinValue, double MaxValue, double Tolerance = 1e-6);

    virtual  int commit(void);
    virtual  int revertToLastCommit(void);
    virtual  int update(void);
    virtual  int update(double newTime, double dT);
    
    virtual  int analysisStep(double dT);

    virtual int hasDomainChanged(void);
    virtual bool getDomainChangeFlag(void);    
    virtual void domainChange(void);    
    virtual void setDomainChangeStamp(int newStamp);

	virtual int  addRecorder(HTRecorder &theRecorder);    	
	virtual int  removeRecorders(void);
	virtual int  removeRecorder(int tag);

  protected:     

  private:
    double currentTime;               // current pseudo time
    double committedTime;             // the committed pseudo time
    double dT;                        // difference between committed and current time
    int	currentGeoTag;             // an integer used to mark if domain has changed
    bool hasDomainChangedFlag;      // a bool flag used to indicate if GeoTag needs to be ++
    //int theDbTag;                   // the Domains unique database tag == 0
    
    TaggedObjectStorage* theElements;
    TaggedObjectStorage* theNodes;
    TaggedObjectStorage* theTempBCs; 
	TaggedObjectStorage* theMPTempBCs;
    TaggedObjectStorage* thePatterns;             

    HT_ElementIter* theEleIter;
    HT_NodeIter* theNodIter;
    TempBCIter* theTempBC_Iter;
	MP_TempBCIter* theMP_TempBCIter;
    BoundaryPatternIter* thePatternIter;        
    AllTempBCIter* allTempBC_Iter; 
    int commitTag;
	Vector theBounds;

	HTRecorder** theRecorders;
	int numRecorders;    
};

#endif
