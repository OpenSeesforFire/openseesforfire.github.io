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

// This class is a modified version of SolutionAlgorithm 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// SolutionAlgorithm
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef HT_SolutionAlgorithm_h
#define HT_SolutionAlgorithm_h

class HeatTransferIntegrator;
class HT_AnalysisModel;
class LinearSOE;
class HT_ConvergenceTest;

class HT_SolutionAlgorithm
{
  public:
    HT_SolutionAlgorithm();
    virtual ~HT_SolutionAlgorithm();

    // public functions defined for subclasses
    virtual void setLinks(HT_AnalysisModel& theModel, 
		                  HeatTransferIntegrator& theIntegrator,
		                  LinearSOE& theSOE,
		                  HT_ConvergenceTest* theTest);
    virtual int domainChanged(void);
    // virtual functions
    virtual int setConvergenceTest(HT_ConvergenceTest* theNewTest);    
    virtual HT_ConvergenceTest* getConvergenceTest(void);       
	virtual int solveCurrentStep(void) = 0;

    virtual int getNumFactorizations(void) {return 0;}
    virtual int getNumIterations(void) {return 0;}
    virtual double getTotalTimeCPU(void)   {return 0.0;}
    virtual double getTotalTimeReal(void)  {return 0.0;}
    virtual double getSolveTimeCPU(void)   {return 0.0;}
    virtual double getSolveTimeReal(void)  {return 0.0;}
    //virtual double getAccelTimeCPU(void)   {return 0.0;}
    //virtual double getAccelTimeReal(void)  {return 0.0;}
 
    // the following are not protected as convergence test
    // may need access to them
    HT_AnalysisModel* getHT_AnalysisModelPtr(void) const;
    HeatTransferIntegrator* getIntegratorPtr(void) const;
    LinearSOE* getLinearSOEptr(void) const;

  protected:
	HT_ConvergenceTest* theTest;

  private:
    HT_AnalysisModel* theModel;
    HeatTransferIntegrator* theIntegrator;
    LinearSOE* theSysOfEqn;
};

#endif
