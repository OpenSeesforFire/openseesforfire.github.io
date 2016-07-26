/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/MP_TemperatureBC.h,v $
                                                                        
                                                                        
#ifndef MP_TemperatureBC_h
#define MP_TemperatureBC_h

// Written: fmk 
// Created: 11/96
// Revision: A
// Modified:
//
// Purpose: This file contains the class definition for MP_TemperatureBC.
// MP_TemperatureBC is a class which stores the information for a multi
// point constraint. A multipoint constraint relates certain dof at 
// a constrained node to be related to certain dof at a retained node: 
//                      {Uc} = [Ccr] {Ur}
//
// The MP_TemperatureBC class assumes time invariant constraints, i.e. the
// constraint matrix does not change over time. All the methods are declared
// as pure virtual, which will allow subclasses for time varying constraints.
//
// What: "@(#) MP_TemperatureBC, revA"
// Modified by Liming for Creating Multipoint TemperatureBC. 2012

#include <HeatTransferDomainComponent.h>

class Matrix;
class ID;

class MP_TemperatureBC : public HeatTransferDomainComponent
{
  public:
    MP_TemperatureBC(int nodeRetain, 
		  int nodeConstr, 
		  Matrix &constrnt,
		  ID &constrainedDOF,
    		  ID &retainedDOF);
  
    MP_TemperatureBC(int nodeRetain,
                   int nodeConstr);

    // destructor    
    virtual ~MP_TemperatureBC();

    // method to get information about the constraint
    virtual int getNodeRetained(void) const;
    virtual int getNodeConstrained(void) const;    
    virtual const ID &getConstrainedDOFs(void) const;        
    virtual const ID &getRetainedDOFs(void) const;            
    virtual int applyConstraint(double pseudoTime);
    virtual bool isTimeVarying(void) const;
    virtual const Matrix &getConstraint(void);    

    

  protected:
    
  private:
    int nodeRetained;        // to identify the retained node
    int nodeConstrained;     // to identify  the constrained node
    Matrix *constraint;      // pointer to the constraint matrix
    ID *constrDOF;           // ID of constrained DOF at constrained node
    ID *retainDOF;           // ID of related DOF at retained node
    
    int dbTag1, dbTag2;      // need a dbTag for the two ID's
};

#endif

