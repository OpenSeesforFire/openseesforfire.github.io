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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-11 21:32:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/VariableTimeStepStaticAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/VariableTimeStepStaticAnalysis.C
// 
// Written: fmk 
// Created: 10/00
// Revision: A
//
// Description: This file contains the implementation of the
// VariableTimeStepStaticAnalysis class.
//
// What: "@(#) VariableTimeStepStaticAnalysis.C, revA"

#include <VariableTimeStepStaticAnalysis.h>
#include <EquiSolnAlgo.h>
#include <StaticIntegrator.h>
#include <LoadControl.h>
#include <Domain.h>
#include <ConvergenceTest.h>
#include <float.h>
#include <AnalysisModel.h>

// Constructor
VariableTimeStepStaticAnalysis::VariableTimeStepStaticAnalysis(
			      Domain &the_Domain,
			      ConstraintHandler &theHandler,
			      DOF_Numberer &theNumberer,
			      AnalysisModel &theModel,
			      EquiSolnAlgo &theSolnAlgo,		   
			      LinearSOE &theLinSOE,
                  StaticIntegrator &theStaticIntegrator,
			      ConvergenceTest *theTest)

:StaticAnalysis(the_Domain, theHandler, theNumberer, theModel, 
			   theSolnAlgo, theLinSOE, theStaticIntegrator, theTest)
{
    domainStamp = 0;
}    

VariableTimeStepStaticAnalysis::~VariableTimeStepStaticAnalysis()
{

}    

int 
VariableTimeStepStaticAnalysis::analyze(int numSteps, double dT, double dtmin, double dtmax, int Jd)
{
  // get some pointers
  Domain *theDom = this->getDomainPtr();
  EquiSolnAlgo *theAlgo = this->getAlgorithm();
  StaticIntegrator*theIntegratr = this->getIntegrator();
  LoadControl* theLoadCtrl = (LoadControl*)theIntegratr;

  ConvergenceTest *theTest = theAlgo->getConvergenceTest();
  AnalysisModel *theModel = this->getModel();


  // set some variables
  int result = 0;  
  double totalTimeIncr = numSteps * dT;
  double currentTimeIncr = 0.0;
  double currentDt = dT;
  double dtMin = dtmin;
  double dtMax = dtmax;
  bool failure = false;

  // loop until analysis has performed the total time incr requested
  while (currentTimeIncr < totalTimeIncr) {
      
      opserr << "Current time:" << theDom->getCurrentTime() << ", dt: " << currentDt << endln;

    if (theModel->analysisStep() < 0) {
      opserr << "DirectIntegrationAnalysis::analyze() - the AnalysisModel failed in newStepDomain";
      opserr << " at time " << theDom->getCurrentTime() << endln;
      theDom->revertToLastCommit();
      return -2;
    }

    int stamp = theDom->hasDomainChanged();


    if (stamp != domainStamp) {
        domainStamp = stamp;

        result = this->domainChanged();

        if (result < 0) {
            opserr << "StaticAnalysis::analyze() - domainChanged failed";
            opserr << " at time step " << theDom->getCurrentTime() << endln;
            return -1;
        }
    }

    //if (this->checkDomainChange() != 0) {
    //  opserr << "VariableTimeStepStaticAnalysis::analyze() - failed checkDomainChange\n";
     // return -1;
    //}

    //
    // do newStep(), solveCurrentStep() and commit() as in regular
    // DirectINtegrationAnalysis - difference is we do not return
    // if a failure - we stop the analysis & resize time step if failure
    //

    theLoadCtrl->setDeltaLambda(currentDt);
    if (theIntegratr->newStep() < 0) {
      result = -2;
    }


    if (result >= 0) {
      result = theAlgo->solveCurrentStep();
      if (result < 0) 
  	    result = -3;
    }    
/*
    // AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY

    TransientIntegrator* theIntegrator = this->getIntegrator();

    if (theIntegrator->shouldComputeAtEachStep()) {
	
      result = theIntegrator->computeSensitivities();
      if (result < 0) {
	opserr << "VariableTimeStepStaticAnalysis::analyze() - the SensitivityAlgorithm failed";
	opserr << " at time ";
	opserr << theDom->getCurrentTime() << endln;
	theDom->revertToLastCommit();
	theIntegrator->revertToLastStep();
	return -5;
      }    
    }

#endif
    // AddingSensitivity:END //////////////////////////////////////
*/
    if (result >= 0) {
      result = theLoadCtrl->commit();
      if (result < 0) 
	result = -4;
    }

    // if the time step was successful increment delta T for the analysis
    // otherwise revert the Domain to last committed state & see if can go on

    if (result >= 0) {
        currentTimeIncr += currentDt;
        failure = false;
        currentDt = this->determineDt(currentDt, dtMin, dtMax, failure, theTest);
    }
      
    else {

      // invoke the revertToLastCommit
      theDom->revertToLastCommit();	    
      theLoadCtrl->revertToLastStep();
      double dt = currentDt;
      failure = true;
      currentDt = this->determineDt(dt, dtMin, dtMax, failure, theTest);
      opserr << "Failure- Current time:" << theDom->getCurrentTime() << ", dt: " << currentDt << endln;

      // if last dT was <= min specified the analysis FAILS - return FAILURE
      if (currentDt <= dtMin) {
	opserr << "VariableTimeStepStaticAnalysis::analyze() - ";
	opserr << " failed at time " << theDom->getCurrentTime() << endln;
	return result;
      }

      
      // if still here reset result for next loop
      result = 0;
    }

    // now we determine a new delta T for next loop
   

  }


  return 0;
}





double 
VariableTimeStepStaticAnalysis::determineDt(double dT,
						       double dtMin, 
						       double dtMax, 
						       bool failure, 
    ConvergenceTest* theTest)
{
  double newDt = dT;
    
  // get the number of trial steps in the last solveCurrentStep()
  double numLastIter = 1.0;
  if (theTest != 0)
    numLastIter = theTest->getNumTests();
  
  
  // determine new dT based on last dT and Jd and #iter of last step
  double factor = 1.0;
  
  if(failure)
      factor = 0.5;
  else {
      if (numLastIter < 10)
          factor = 2.0;
      else if (numLastIter < 100)
          factor = 1.0;
  }
  

  newDt *= factor;
  
  // ensure: dtMin <~~ dT <= dtMax
  if (newDt < dtMin)
    newDt = dtMin;  // to ensure we get out of the analysis 
                               // loop if can't converge on next step
  else if (newDt > dtMax)
    newDt = dtMax;
    
  return newDt;
}


