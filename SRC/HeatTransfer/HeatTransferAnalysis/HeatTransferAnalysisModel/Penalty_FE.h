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

// This class is a modified version of PenaltySP_FE 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// PenaltySP_FE
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef Penalty_FE_h
#define Penalty_FE_h

#include <HT_FE_Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class HeatTransferElement;
class HeatTransferIntegrator;
class HT_TransientIntegrator;
class HT_AnalysisModel;
class HeatTransferDomain;
class TemperatureBC;
class HeatTransferNode;

class Penalty_FE: public HT_FE_Element
{
  public:
    Penalty_FE(int tag, 
		       HeatTransferDomain& theDomain, 
		       TemperatureBC& theTempBC, 
		       double pnum = 1.0e8);    

    virtual ~Penalty_FE();    

    // public methods
    virtual int  setID(void);
    virtual const Matrix& getTangent(HeatTransferIntegrator* theIntegrator);
    virtual const Vector& getResidual(HeatTransferIntegrator* theIntegrator);

  protected:
    
  private:
    double penalty_num;
    TemperatureBC* the_temp_bc;
    HeatTransferNode* the_node;
    static Matrix tang;
    static Vector resid;
};

#endif
