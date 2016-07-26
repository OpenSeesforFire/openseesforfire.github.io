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

// This class is a modified version of SP_Constraint 
// for the heat transfer module
// Modified by: Yaqiang Jiang(y.jiang@ed.ac.uk)

// SP_Constraint
// Written: fmk 
// Created: 11/96
// Revision: A

#ifndef TemperatureBC_h
#define TemperatureBC_h

#include <HeatTransferDomainComponent.h>

//class Vector;
class HeatTransferNode;

class TemperatureBC : public HeatTransferDomainComponent
{
  public:
    // constructors      
    TemperatureBC(int nodeTag, int ndof, double value, bool isConstant = true);  // need to determine if keep ndof or not??
    virtual ~TemperatureBC();

    virtual int getNodeTag(void) const;
    virtual int getDOF_Number(void) const;
    virtual int applyTemperatureBC(double factor);    
    virtual double getValue(void);

    virtual void setPatternTag(int patternTag);
    virtual int  getPatternTag(void) const;
	virtual void setNodalValue();

  protected:
    int nodeTag;     // to identify the node in the model
    int dofNumber;   // identifies which of the nodes dof is constrrained 
    double valueR;   // the reference value
    double valueC;   // if constant = the reference value, if not contant =
	                 // the reference value * load factor
    bool isConstant; // flag indicating if constant
    int  patternTag;    
	//HeatTransferNode* theNode;
	//Vector* nodalResponses;
};

#endif
