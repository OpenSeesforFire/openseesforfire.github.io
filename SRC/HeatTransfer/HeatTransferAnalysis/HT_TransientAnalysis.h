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

// This class is a modified version of TransientAnalysis 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// TransientAnalysis
// Written: fmk 
// Created: 11/96
// Revision: A


#ifndef HT_TransientAnalysis_h
#define HT_TransientAnalysis_h

#include <HeatTransferAnalysis.h>

class TemperatureBCHandler;
class HT_DOF_Numberer;
class HT_AnalysisModel;
class HT_TransientIntegrator;
class LinearSOE;
class HT_SolutionAlgorithm;
class HT_ConvergenceTest;


class HT_TransientAnalysis: public HeatTransferAnalysis
{
  public:
    HT_TransientAnalysis(HeatTransferDomain& theDomain, 
			             TemperatureBCHandler& theHandler,
			             HT_DOF_Numberer& theNumberer,
			             HT_AnalysisModel& theModel,
			             HT_SolutionAlgorithm& theSolnAlgo,		   
			             LinearSOE& theSOE,
			             HT_TransientIntegrator& theIntegrator,
			             HT_ConvergenceTest* theTest = 0);

    virtual ~HT_TransientAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps, double dT,double& lasttime, double monitortime = 0 );
    int initialize(void);
    int domainChanged(void);

    int setNumberer(HT_DOF_Numberer& theNumberer);    
    int setAlgorithm(HT_SolutionAlgorithm& theAlgorithm);
    int setIntegrator(HT_TransientIntegrator& theIntegrator);
    int setLinearSOE(LinearSOE& theSOE); 
    int setConvergenceTest(HT_ConvergenceTest& theTest);
    
    int checkDomainChange(void);

    HT_SolutionAlgorithm* getAlgorithm(void);
    HT_TransientIntegrator* getIntegrator(void);
    HT_ConvergenceTest* getConvergenceTest(void); 
    HT_AnalysisModel* getModel(void) ;
    
  protected:
    
  private:
    TemperatureBCHandler*  temp_bc_handler;    
    HT_DOF_Numberer*  dof_numberer;
    HT_AnalysisModel*  analysis_model;
    HT_SolutionAlgorithm*  solution_algorithm;
    LinearSOE*  linear_soe;
    HT_TransientIntegrator*  transient_integrator;
    HT_ConvergenceTest*  convergence_test;

    int domainStamp;
};

#endif
