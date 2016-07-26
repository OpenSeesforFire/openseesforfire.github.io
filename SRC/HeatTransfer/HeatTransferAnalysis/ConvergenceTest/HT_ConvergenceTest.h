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


// This class is a modified version of ConvergenceTest 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// ConvergenceTest
// Written: fmk 
// Date: 09/98

#ifndef HT_ConvergenceTest_h
#define HT_ConvergenceTest_h

//#include <MovableObject.h>
#include <Vector.h>
#include <bool.h>

class HT_SolutionAlgorithm;


class HT_ConvergenceTest
{
  public:
    // constructors and destructor
    HT_ConvergenceTest();	
    virtual ~HT_ConvergenceTest();

    virtual HT_ConvergenceTest* getCopy( int iterations ) = 0 ;

	virtual int setHT_SolutionAlgorithm(HT_SolutionAlgorithm& HT_Algo) = 0;


    virtual int start(void) =0;
    virtual int test(void) = 0;
  
    virtual int getNumTests(void) =0;    
    virtual int getMaxNumTests(void) =0;        
    virtual double getRatioNumToMax(void) =0;            
    virtual const Vector& getNorms(void) =0;
    
    
  protected:

  private:
};


#endif

