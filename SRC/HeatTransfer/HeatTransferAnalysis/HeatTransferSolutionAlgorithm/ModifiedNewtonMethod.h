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

// This class is a modified version of ModifiedNewton 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// ModifiedNewton
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef ModifiedNewtonMethod_h
#define ModifiedNewtonMethod_h

#include <HT_SolutionAlgorithm.h>
//#include <Vector.h>


class ModifiedNewtonMethod: public HT_SolutionAlgorithm
{
  public:
    ModifiedNewtonMethod();    
    ModifiedNewtonMethod(HT_ConvergenceTest& theTest);
    ~ModifiedNewtonMethod();

    //int setConvergenceTest(HT_ConvergenceTest* theNewTest);
	//HT_ConvergenceTest* getConvergenceTest();
    int solveCurrentStep(void);     
    
  protected:

  private:

  
};

#endif
