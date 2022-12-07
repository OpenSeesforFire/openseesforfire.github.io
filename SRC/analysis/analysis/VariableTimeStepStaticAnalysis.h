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
                                                                        

// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/VariableTimeStepStaticAnalysis.h,v $
                                                                        
                                                                        
#ifndef VariableTimeStepStaticAnalysis_h
#define VariableTimeStepStaticAnalysis_h

// Written: fmk 
// Created: 10/00
// Revision: A
//
// Description: This file contains the class definition for 
// VariableTimeStepStaticAnalysis. VariableTimeStepStaticAnalysis 
// is a subclass of DirectIntegrationAnalysis. It is used to perform a 
// dynamic analysis on the FE\_Model using a direct integration scheme.  
//
// What: "@(#) VariableTimeStepStaticAnalysis.h, revA"

#include <StaticAnalysis.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class StaticIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;

class VariableTimeStepStaticAnalysis: public StaticAnalysis
{
  public:
    VariableTimeStepStaticAnalysis(Domain &theDomain,
					      ConstraintHandler &theHandler,
					      DOF_Numberer &theNumberer,
					      AnalysisModel &theModel,
					      EquiSolnAlgo &theSolnAlgo,
					      LinearSOE &theSOE,
                          StaticIntegrator &theIntegrator,
					      ConvergenceTest *theTest =0);
    virtual ~VariableTimeStepStaticAnalysis();

    int analyze(int numSteps, double dT, double dtMin, double dtMax, int Jd);
 

  protected:
    double determineDt(double dT, double dtMin, double dtMax, bool failure,
        ConvergenceTest* theTest);

  private:
      int domainStamp;
};

#endif

