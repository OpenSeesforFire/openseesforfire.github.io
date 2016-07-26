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

// This class is a modified version of Linear 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// Linear
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 

#ifndef LinearAlgorithm_h
#define LinearAlgorithm_h

#include <HT_SolutionAlgorithm.h>


class LinearAlgorithm: public HT_SolutionAlgorithm
{
  public:
    LinearAlgorithm();
    ~LinearAlgorithm();

    int solveCurrentStep(void);
    int setConvergenceTest(HT_ConvergenceTest* theNewTest);
    
  protected:
    
  private:

};

#endif
