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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/MP_TemperatureBC.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of class MP_TemperatureBC.
//
// The class MP_TemperatureBC interface:
// Modified by Liming for crating multi-point temperature boundary condition;

#include <MP_TemperatureBC.h>

#include <stdlib.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
//#include <FEM_ObjectBroker.h>

static int numMPs = 0;
static int nextTag = 0;
 


// general constructor for ModelBuilder
MP_TemperatureBC::MP_TemperatureBC(int nodeRetain, int nodeConstr, Matrix &constr,
			     ID &constrainedDOF, ID &retainedDOF)
:HeatTransferDomainComponent(nextTag++), 
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0), dbTag1(0), dbTag2(0)
{
  numMPs++;    
  constrDOF = new ID(constrainedDOF);
  retainDOF = new ID(retainedDOF);    
  if (constrDOF == 0 || constrainedDOF.Size() != constrDOF->Size() ||
      retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
    opserr << "MP_TemperatureBC::MP_TemperatureBC - ran out of memory 1\n";
    exit(-1);
  }    
  
  constraint = new Matrix(constr);
  if (constraint == 0 || constr.noCols() != constr.noCols()) { 
    opserr << "MP_TemperatureBC::MP_TemperatureBC - ran out of memory 2\n";
    exit(-1);
  }        
}


MP_TemperatureBC::MP_TemperatureBC(int nodeRetain, int nodeConstr)
:HeatTransferDomainComponent(nextTag++),
nodeRetained(nodeRetain), nodeConstrained(nodeConstr),
constraint(0), constrDOF(0), retainDOF(0), dbTag1(0), dbTag2(0)
{
  numMPs++;
  constrDOF= new ID(1);
  (*constrDOF)(0)=0;
  retainDOF = new ID(1);
  (*retainDOF)(0)=0;
  
  constraint = new Matrix(1,1);
  (*constraint)(0,0)=1.0;      
}



MP_TemperatureBC::~MP_TemperatureBC()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	delete constraint;
    if (constrDOF != 0)
	delete constrDOF;
    if (retainDOF != 0)
	delete retainDOF;    
    
    numMPs--;
    if (numMPs == 0)
      nextTag = 0;
}


int
MP_TemperatureBC::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_TemperatureBC::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_TemperatureBC::getConstrainedDOFs(void) const
{
    if (constrDOF == 0) {
	opserr << "MP_TemperatureBC::getConstrainedDOF - no ID was set, ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";	
	exit(-1);
    }

    // return the ID corresponding to constrained DOF of Ccr
    return *constrDOF;    
}


const ID &
MP_TemperatureBC::getRetainedDOFs(void) const
{
    if (retainDOF == 0) {
	opserr << "MP_TemperatureBC::getRetainedDOFs - no ID was set\n ";
	opserr << "was recvSelf() ever called? or subclass incorrect?\n";		
	exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return *retainDOF;    
}

int 
MP_TemperatureBC::applyConstraint(double timeStamp)
{
    // does nothing MP_TemperatureBC objects are time invariant
    return 0;
}

bool
MP_TemperatureBC::isTimeVarying(void) const
{
    return false;
}


const Matrix &
MP_TemperatureBC::getConstraint(void)
{
    if (constraint == 0) {
	opserr << "MP_TemperatureBC::getConstraint - no Matrix was set\n";
	exit(-1);
    }    

    // return the constraint matrix Ccr
    return *constraint;    
}




