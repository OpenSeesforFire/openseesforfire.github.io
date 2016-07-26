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
// Note: This class was adapted from SolutionAlgorithm 
#include <HT_SolutionAlgorithm.h>
#include <HT_AnalysisModel.h>
#include <HeatTransferIntegrator.h>
#include <LinearSOE.h>
#include <HT_ConvergenceTest.h>

HT_SolutionAlgorithm::HT_SolutionAlgorithm()
:theModel(0), theIntegrator(0), theSysOfEqn(0),
 theTest(0)
{

}

HT_SolutionAlgorithm::~HT_SolutionAlgorithm()
{

}

void 
HT_SolutionAlgorithm::setLinks(HT_AnalysisModel& theNewModel, 
		                       HeatTransferIntegrator& theNewIntegrator,
		                       LinearSOE& theSOE,
		                       HT_ConvergenceTest* theConvergenceTest)
{
    theModel = &theNewModel;
    theIntegrator = &theNewIntegrator;
    theSysOfEqn = &theSOE;
	theTest = theConvergenceTest;

    this->setConvergenceTest(theConvergenceTest);
}


int 
HT_SolutionAlgorithm::domainChanged()
{
    return 0;
}


int 
HT_SolutionAlgorithm::setConvergenceTest(HT_ConvergenceTest* theConvergenceTest)
{
    theTest = theConvergenceTest;
    return 0;
}

HT_ConvergenceTest*
HT_SolutionAlgorithm::getConvergenceTest(void)
{
    return 0;
}


HT_AnalysisModel*
HT_SolutionAlgorithm::getHT_AnalysisModelPtr(void) const
{
    return theModel;
}



HeatTransferIntegrator*
HT_SolutionAlgorithm::getIntegratorPtr(void) const
{
    return theIntegrator;
}



LinearSOE*
HT_SolutionAlgorithm::getLinearSOEptr(void) const
{
    return theSysOfEqn;
}