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
// Note: This class was adapted from Linear 
#include <LinearAlgorithm.h>
#include <HT_AnalysisModel.h>
//#include <StaticAnalysis.h>
#include <HeatTransferIntegrator.h>
#include <LinearSOE.h>
#include <Vector.h>
//#include <Channel.h>
//#include <FEM_ObjectBroker.h>
#include <HT_ConvergenceTest.h>

//#include <Timer.h>
// Constructor
LinearAlgorithm::LinearAlgorithm()
{

}

// Destructor
LinearAlgorithm::~LinearAlgorithm()
{

}

// int run(void)
//    Performs the linear solution algorithm.

int 
LinearAlgorithm::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass

    HT_AnalysisModel* theAnalysisModel = this->getHT_AnalysisModelPtr(); 
    LinearSOE* theSOE = this->getLinearSOEptr();
    HeatTransferIntegrator* theIntegrator = this->getIntegratorPtr(); 

    if ((theAnalysisModel == 0) || (theIntegrator ==0 ) || (theSOE == 0)){
	opserr << "WARNING LinearAlgorithm::solveCurrentStep() -";
	opserr << "setLinks() has not been called.\n";
	return -5;
	}

    if (theIntegrator->formTangent() < 0) {
	opserr << "WARNING LinearAlgorithm::solveCurrentStep() -";
	opserr << "the Integrator failed in formTangent()\n";
	return -1;
    }	

    
    if (theIntegrator->formUnbalance() < 0) {
	opserr << "WARNING LinearAlgorithm::solveCurrentStep() -";
	opserr << "the Integrator failed in formUnbalance()\n";	
	return -2;
    }

    if (theSOE->solve() < 0) {
	opserr << "WARNING LinearAlgorithm::solveCurrentStep() -";
	opserr << "the LinearSOE failed in solve()\n";	
	return -3;
    }

    const Vector& deltaT = theSOE->getX();

    if (theIntegrator->update(deltaT) < 0) {
	opserr << "WARNING LinearAlgorithm::solveCurrentStep() -";
	opserr << "the Integrator failed in update()\n";	
	return -4;
    }

    return 0;
}

int
LinearAlgorithm::setConvergenceTest(HT_ConvergenceTest* theNewTest)
{
    return 0;
}
