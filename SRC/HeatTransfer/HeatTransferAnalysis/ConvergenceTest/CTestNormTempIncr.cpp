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
// Note: This class was adapted from CTestNormDispIncr   
#include <CTestNormTempIncr.h>
#include <Vector.h>
#include <HT_SolutionAlgorithm.h>
#include <LinearSOE.h>


CTestNormTempIncr::CTestNormTempIncr()	    	
:theSOE(0), tol(0), maxNumIter(0), currentIter(0), printFlag(0), 
 norms(25), nType(2)
{
    
}


CTestNormTempIncr::CTestNormTempIncr(double theTol, int maxIter, int printIt, int normType)
:theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
 norms(maxIter), nType(normType)
{
    
}


CTestNormTempIncr::~CTestNormTempIncr()
{
    
}


HT_ConvergenceTest* CTestNormTempIncr::getCopy(int iterations)
{
    CTestNormTempIncr* theCopy ;
    theCopy = new CTestNormTempIncr(this->tol, iterations, this->printFlag, this->nType) ;
    
    theCopy->theSOE = this->theSOE ;
    
    return theCopy ;
}


void CTestNormTempIncr::setTolerance(double newTol)
{
    tol = newTol;
}


int CTestNormTempIncr::setHT_SolutionAlgorithm(HT_SolutionAlgorithm& HT_Algo)
{
    theSOE = HT_Algo.getLinearSOEptr();
    if (theSOE == 0) {
		opserr << "WARNING: CTestNormTempIncr::setHT_SolutionAlgorithm() - no SOE\n";	
		return -1;
		}
	else
		return 0;
}


int CTestNormTempIncr::test(void)
{
    // check to ensure the SOE has been set - this should not happen if the 
    // return from start() is checked
    if (theSOE == 0)
		return -2;
    
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
		opserr << "WARNING: CTestNormTempIncr::test() - start() was never invoked.\n";	
		return -2;
		}
    
    // get the X vector & determine it's norm & save the value in norms vector
    const Vector& x = theSOE->getX();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter) 
		norms(currentIter-1) = norm;
    
    // print the data if required
    if (printFlag == 1) {
        opserr << "CTestNormTempIncr::test() - iteration: " << currentIter;
        opserr << " current Norm: " << norm << " (max: " << tol;
        opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
		} 
    if (printFlag == 4) {
		opserr << "CTestNormTempIncr::test() - iteration: " << currentIter;
		opserr << " current Norm: " << norm << " (max: " << tol << ")\n";
		opserr << "\tNorm deltaX: " << norm << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << endln;
		opserr << "\tdeltaX: " << x << "\tdeltaR: " << theSOE->getB();
		} 

    //
    // check if the algorithm converged
    //
    
    // if converged - print & return ok
    if (norm <= tol) { 
		// do some printing first
		if (printFlag != 0) {
			if (printFlag == 1 || printFlag == 4) 
				opserr << endln;
			else if (printFlag == 2 || printFlag == 6) {
				opserr << "CTestNormTempIncr::test() - iteration: " << currentIter;
				opserr << " current Norm: " << norm << " (max: " << tol;
				opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
				}
			}
		// return the number of times test has been called
		return currentIter;
		}
 
    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
		opserr << "WARNING: CTestNormTempIncr::test() - failed to converge but going on - ";
		opserr << " current Norm: " << norm << " (max: " << tol;
		opserr << ", Norm deltaR: " << theSOE->getB().pNorm(nType) << ")\n";
		return currentIter;
		}

    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter) { // failes to converge
		opserr << "WARNING: CTestNormTempIncr::test() - failed to converge \n";
		opserr << "after: " << currentIter << " iterations\n";	
		currentIter++;    
		return -2;
		} 
    // algorithm not yet converged - increment counter and return -1
    else {
		currentIter++;    
		return -1;
		}
}


int CTestNormTempIncr::start(void)
{
    if (theSOE == 0) {
		opserr << "WARNING: CTestNormTempIncr::test() - no SOE returning true\n";
		return -1;
		}   
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    return 0;
}


int CTestNormTempIncr::getNumTests()
{
    return currentIter;
}


int CTestNormTempIncr::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestNormTempIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestNormTempIncr::getNorms() 
{
    return norms;
}
