                                                                                                                                            
#ifndef BoundaryPattern_h
#define BoundaryPattern_h

#include <HeatTransferDomainComponent.h>
#include <Vector.h>

class NodalLoad;
class TimeSeries;
class ElementalLoad;
class SP_Constraint;
class NodalLoadIter;
class ElementalLoadIter;
class SingleDomSP_Iter;
class SP_ConstraintIter;
class TaggedObjectStorage;
class GroundMotion;

class LoadPattern : public HeatTransferDomainComponent    
{
  public:
    // constructors
    LoadPattern(int tag);    
    LoadPattern();                      // for FEM_ObjectBroker
    LoadPattern(int tag, int classTag); // for subclasses
    
    // destructor
    virtual ~LoadPattern();

    // method to set the associated TimeSeries and Domain
    virtual void setTimeSeries(TimeSeries *theSeries);
    virtual void setDomain(Domain *theDomain);

    // methods to add loads
    virtual bool addSP_Constraint(SP_Constraint *);
    virtual bool addNodalLoad(NodalLoad *);
    virtual bool addElementalLoad(ElementalLoad *);
    virtual NodalLoadIter     &getNodalLoads(void);
    virtual ElementalLoadIter &getElementalLoads(void);    
    virtual SP_ConstraintIter &getSPs(void);        
    
    // methods to remove loads
    virtual void clearAll(void);
    virtual NodalLoad *removeNodalLoad(int tag);
    virtual ElementalLoad *removeElementalLoad(int tag);
    virtual SP_Constraint *removeSP_Constraint(int tag);

    // methods to apply loads
    virtual void applyBCs(double time = 0.0);
    virtual void setLoadConstant(void);
	virtual void unsetLoadConstant(void);
    virtual double getLoadFactor(void);
    

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void);

    virtual int addMotion(GroundMotion &theMotion, int tag);    
    virtual GroundMotion *getMotion(int tag);        

  protected:
    int    isConstant;     // to indictae whether setConstant has been called
	
  private:
    double loadFactor;     // current load factor

    TimeSeries *theSeries; // pointer to associated TimeSeries

    int	   currentGeoTag;
    int    lastGeoSendTag;
    int    dbSPs, dbNod, dbEle; // database tags for storing info about components
    
    // storage objects for the loads and constraints
    TaggedObjectStorage  *theNodalLoads;
    TaggedObjectStorage  *theElementalLoads;
    TaggedObjectStorage  *theSPs; 	  

    // iterator objects for the objects added to the storage objects
    NodalLoadIter       *theNodIter;
    ElementalLoadIter   *theEleIter;
    SingleDomSP_Iter    *theSpIter;    

    // AddingSensitivity:BEGIN //////////////////////////////////////
    Vector *randomLoads;
    bool RVisRandomProcessDiscretizer;
    // AddingSensitivity:END ////////////////////////////////////////

    int lastChannel; 
};

#endif







