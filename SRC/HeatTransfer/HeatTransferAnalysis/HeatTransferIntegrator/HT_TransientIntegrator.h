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

// This class is a modified version of TransientIntegrator 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// TransientIntegrator
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A

#ifndef TransientIntegrator_h
#define TransientIntegrator_h

#include <HeatTransferIntegrator.h>
class LinearSOE;
class HT_AnalysisModel;
class HT_FE_Element;
class Vector;

class HT_TransientIntegrator: public HeatTransferIntegrator
{
  public:
	  HT_TransientIntegrator();
	  virtual ~HT_TransientIntegrator();
	  virtual int formEleResidual(HT_FE_Element* theFE_Ele);
	  virtual int initialize(void) {return 0;};
	  virtual bool getUpdatingFlag() = 0;
	  virtual double getAlphaDeltat() = 0;

  protected:
    
  private:

};

#endif
