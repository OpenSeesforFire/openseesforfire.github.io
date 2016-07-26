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
// Note: This class was adapted from ModifiedNewton  

#include <ModifiedNewtonMethod.h>
#include <HT_AnalysisModel.h>
#include <HeatTransferIntegrator.h>
#include <LinearSOE.h>
#include <HT_ConvergenceTest.h>

// Constructor
ModifiedNewtonMethod::ModifiedNewtonMethod()
{
    // does nothing  
}


ModifiedNewtonMethod::ModifiedNewtonMethod(HT_ConvergenceTest& theT)
{
    theTest = &theT;
}

// Destructor
ModifiedNewtonMethod::~ModifiedNewtonMethod()
{
    // does nothing
}


int 
ModifiedNewtonMethod::solveCurrentStep(void)
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    HT_AnalysisModel*  theModel = this->getHT_AnalysisModelPtr();
    HeatTransferIntegrator*  theIntegrator = this->getIntegratorPtr();
    LinearSOE*  theSOE = this->getLinearSOEptr();

    if ((theModel == 0) || (theIntegrator == 0) || (theSOE == 0)
	|| (theTest == 0)){
		opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() - setLinks() has";
		opserr << " not been called - or no HT_ConvergenceTest has been set\n";
		return -5;
		}	

    if (theIntegrator->formUnbalance() < 0) {
		opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() -";
		opserr << "the Integrator failed in formUnbalance()\n";	
		return -2;
		}	

    if (theIntegrator->formTangent() < 0){
		opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() -";
		opserr << "the Integrator failed in formTangent()\n";
		return -1;
		}		    

    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setHT_SolutionAlgorithm(*this);
    if (theTest->start() < 0) {
		opserr << "ModifiedNewtonMethod::solveCurrentStep() -";
		opserr << "the HT_ConvergenceTest object failed in start()\n";
		return -3;
		}

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    int count = 0;
    do {
		//Timer timer2;
		//timer2.start();
		if (theSOE->solve() < 0) {
			opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() -";
			opserr << "the LinearSOE failed in solve()\n";	
			return -3;
			}	    

		if (theIntegrator->update(theSOE->getX()) < 0) {
			opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() -";
			opserr << "the HeatTransferIntegrator failed in update()\n";	
			return -4;
			}	        

		if (theIntegrator->formUnbalance() < 0) {
			opserr << "WARNING ModifiedNewtonMethod::solveCurrentStep() -";
			opserr << "the HeatTransferIntegrator failed in formUnbalance()\n";	
			return -2;
			}	

		//this->record(count++);
		result = theTest->test();

		} while (result == -1);

    if (result == -2) {
		opserr << "ModifiedNewtonMethod::solveCurrentStep() -";
		opserr << "the HT_ConvergenceTest object failed in test()\n";
		return -3;
		}

	opserr << "ModifiedNewtonMethod::solveCurrentStep() - invoked\n";
    return result;
}
