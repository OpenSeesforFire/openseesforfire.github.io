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

//
// Written by Yaqiang Jiang (y.jiang@ed.ac.uk)
//
// Note: This class was adapted from EarthquakePattern

#ifndef FireImposedPattern_h
#define FireImposedPattern_h

#include <BoundaryPattern.h>
#include <Vector.h>

class TemperatureBC;
//class TimeSeries;
class FireModel;
class HeatFluxBC;
class TemperatureBCIter;
class TempBCIter;
class HeatFluxBCIter;
class HTDomain_Iter;
class TaggedObjectStorage;

class FireImposedPattern : public BoundaryPattern    
{
  public:
    // constructors
    FireImposedPattern(int tag);    
    virtual ~FireImposedPattern();

    // method to set the associated TimeSeries and Domain
    //virtual void setTimeSeries(TimeSeries *theSeries);
	virtual void setFireModel(FireModel* theFireModel);
	//virtual void removeFireModel();
    //virtual void setDomain(HeatTransferDomain *theDomain);

    // methods to add loads
    //virtual bool addTemperatureBC(TemperatureBC* );
    //virtual bool addHeatFluxBC(HeatFluxBC* );
/*    virtual HeatFluxBCIter& getHeatFluxBCs(void);    
    virtual TemperatureBCIter& getTemperatureBCs(void);     */   
    
    // methods to remove loads
    //virtual void clearAll(void);
    //virtual HeatFluxBC* removeHeatFluxBC(int tag);
    //virtual TemperatureBC* removeTemperatureBC(int tag);

    // methods to apply loads
    virtual void applyBCs(double time = 0.0);
 //   virtual void setFluxConstant(void);
	//virtual void unsetFluxConstant(void);
 //   virtual double getFactor(void);
    
    // method to obtain a blank copy of the BoundaryPattern
    //virtual BoundaryPattern* getCopy(void);     

  protected:
    //int isConstant;     // to indictae whether setConstant has been called
	
  private:
    //double factor;     // current load factor

    FireModel* the_firemodel; // pointer to associated TimeSeries
    
    //// storage objects for the fluxes and temperatures
    //TaggedObjectStorage* theHeatFluxBCs;
    //TaggedObjectStorage* theTemperatureBCs; 	  

    // iterator objects for the objects added to the storage objects
    //HeatFluxBCIter* theFluxIter;
    //TempBCIter* theTempIter;    
};

#endif







