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
// Note: This class was adapted from ConstraintHandler   
#ifndef TemperatureBCHandler_h
#define TemperatureBCHandler_h


//#include <MovableObject.h>

class ID;
class HeatTransferDomain;
class HT_AnalysisModel;
class HeatTransferIntegrator;
//class FEM_ObjectBroker;

class TemperatureBCHandler
{
  public:
    TemperatureBCHandler();
    virtual ~TemperatureBCHandler();

    void setLinks(HeatTransferDomain& theDomain, 
		  HT_AnalysisModel& theModel,
		  HeatTransferIntegrator& theIntegrator);

    // pure virtual functions
    virtual int handle() = 0;
    virtual int update(void);
    virtual int applyBCs(void);
    virtual int doneNumberingDOF(void);
    virtual void clearAll(void) =0;    

  protected:
    HeatTransferDomain* getDomainPtr(void) const;
    HT_AnalysisModel* getAnalysisModelPtr(void) const;
    HeatTransferIntegrator* getIntegratorPtr(void) const;
    
  private:
    HeatTransferDomain* theDomainPtr;
    HT_AnalysisModel* theAnalysisModelPtr;
    HeatTransferIntegrator* theIntegratorPtr;
};

#endif
