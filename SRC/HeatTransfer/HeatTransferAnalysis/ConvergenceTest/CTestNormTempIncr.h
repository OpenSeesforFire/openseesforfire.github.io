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

// This class is a modified version of CTestNormDispIncr 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// CTestNormDispIncr
// Written: fmk 
// Date: 09/98
// Modified: 05/05 ahs

#ifndef CTestNormTempIncr_h
#define CTestNormTempIncr_h

#include <HT_ConvergenceTest.h>
#include <bool.h>
class HT_SolutionAlgorithm;
class LinearSOE;


class CTestNormTempIncr: public HT_ConvergenceTest
{
  public:
	  // constructors
	  CTestNormTempIncr();	    	
	  CTestNormTempIncr(double tol, int maxNumIter, int printFlag, int normType = 2);

	  // destructor
	  ~CTestNormTempIncr();

	  HT_ConvergenceTest* getCopy(int iterations);

	  void setTolerance(double newTol);
	  int setHT_SolutionAlgorithm(HT_SolutionAlgorithm& HT_Algo);

	  int test(void);
	  int start(void);

	  int getNumTests(void);
	  int getMaxNumTests(void);        
	  double getRatioNumToMax(void);                
	  const Vector& getNorms(void);  
   
  protected:
    
  private:
	  LinearSOE* theSOE;
	  double tol;         // the tol on the norm used to test for convergence

	  int maxNumIter;     // max number of iterations
	  int currentIter;    // number of times test() has been invokes since last start()
	  int printFlag;      // a flag indicating if to print on test
	  int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)

	  Vector norms;       // vector to hold the norms
};

#endif


