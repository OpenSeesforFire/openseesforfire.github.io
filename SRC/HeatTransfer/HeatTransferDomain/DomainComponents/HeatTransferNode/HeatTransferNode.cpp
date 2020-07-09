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
// Modified by Liming Jiang//
//
// Note: This class was adapted from Node, the DoF should be reduced!!

#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <HeatTransferDomain.h>
#include <HT_DOF_Group.h>
#include <string.h>
#include <Vector.h>
#include <Matrix.h>
#include <stdlib.h>
#include <HT_NodeIter.h>

//#include <OPS_Globals.h>



//  HeatTransferNode(int tag, int ndof, double Crd1, double yCrd);
//	constructor for 2d nodes
HeatTransferNode::HeatTransferNode(int tag, int ndof, double Crd1)
:HeatTransferDomainComponent(tag), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitTemp(0), commitTdot(0), 
 trialTemp(0), trialTdot(0), incrTemp(0),
 incrDeltaTemp(0), T(0), Tdot(0)
{
    Crd = new Vector(2);
    (*Crd)(0) = Crd1;
}

//  HeatTransferNode(int tag, int ndof, double Crd1, double yCrd);
//	constructor for 2d nodes
HeatTransferNode::HeatTransferNode(int tag, int ndof, double Crd1, double Crd2)
:HeatTransferDomainComponent(tag), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitTemp(0), commitTdot(0), 
 trialTemp(0), trialTdot(0), incrTemp(0),
 incrDeltaTemp(0), T(0), Tdot(0)
{
    Crd = new Vector(2);
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
}


//  HeatTransferNode(int tag, int ndof, double Crd1, double Crd2, double zCrd);
//	constructor for 3d nodes
HeatTransferNode::HeatTransferNode(int tag, int ndof, double Crd1, double Crd2, double Crd3)
:HeatTransferDomainComponent(tag), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitTemp(0), commitTdot(0),
 trialTemp(0), trialTdot(0), incrTemp(0),
 incrDeltaTemp(0), T(0), Tdot(0)
{ 
    Crd = new Vector(3);
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
    (*Crd)(2) = Crd3;    
}

// ~Node():
// 	destructor

HeatTransferNode::~HeatTransferNode()
{
    // delete anything that we created with new
    if (Crd != 0)
		delete Crd;

    if (commitTemp != 0)
		delete commitTemp;

    if (commitTdot != 0)
		delete commitTdot;

    if (trialTemp != 0)
		delete trialTemp;

    if (trialTdot != 0)
		delete trialTdot;

    if (incrTemp != 0)
		delete incrTemp;
    
    if (incrDeltaTemp != 0)
		delete incrDeltaTemp;    
    
    if (T != 0)
		delete[] T;

    if (Tdot != 0)
		delete[] Tdot;

    if (theDOF_GroupPtr != 0)
		theDOF_GroupPtr->resetNodePtr();
}


int
HeatTransferNode::getNumberDOF(void) const
{
    // return the number of dof
    return  numberDOF;
}


void
HeatTransferNode::setDOF_GroupPtr(HT_DOF_Group* theDOF_Grp)
{
    // set the DOF_Group pointer
    theDOF_GroupPtr = theDOF_Grp;
}


HT_DOF_Group*
HeatTransferNode::getDOF_GroupPtr(void)
{
    // return the DOF_Group pointer
    return theDOF_GroupPtr;
}


const Vector&
HeatTransferNode::getCrds() const
{
    // return the vector of nodal coordinates
    return *Crd;
}


const Vector&
HeatTransferNode::getTemperature(void) 
{
    // construct memory and Vectors for trial and committed
    // temperature on first call to this method, getTrialTemperature()
    // getTrialTemperature() or incrTrialTemperature()
    if (commitTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::getTemperature() -- ran out of memory\n";
			exit(-1);
			}
		}
    
    // return the committed temperature
    return *commitTemp;
}


const Vector&
HeatTransferNode::getTdot(void) 
{
    // construct memory and Vectors for trial and committed
    // velocity on first call to this method, getTrialTdot()
    // setTrialTdot() or incrTrialTdot()    
    if (commitTdot == 0) {
		if (this->createTdot() < 0) {
			opserr << "FATAL HeatTransferNode::getTdot() -- ran out of memory\n";
			exit(-1);
			}
		}    

    // return the velocity
    return *commitTdot;    
}

/* *********************************************************************
**
**   Methods to return the trial response quantities
**
** *********************************************************************/

const Vector&
HeatTransferNode::getTrialTemperature(void) 
{
    if (trialTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::getTrialTemperature() -- ran out of memory\n";
			exit(-1);
			}
		}    

    return *trialTemp;
}


const Vector&
HeatTransferNode::getTrialTdot(void) 
{
    if (trialTdot == 0) {
		if (this->createTdot() < 0) {
			opserr << "FATAL HeatTransferNode::getTrialTdot() -- ran out of memory\n";
			exit(-1);
			}
		}    

    return *trialTdot;
}


const Vector&
HeatTransferNode::getIncrTemperature(void) 
{
    if (incrTemp == 0) {
	if (this->createTemp() < 0) {
	    opserr << "FATAL HeatTransferNode::getIncrTemperature() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *incrTemp;
}


const Vector&
HeatTransferNode::getIncrDeltaTemperature(void) 
{
    if (incrDeltaTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::getIncrDeltaTemperature() -- ran out of memory\n";
			exit(-1);
			}
		}    

    return *incrDeltaTemp;
}


int
HeatTransferNode::setTrialTemperature(const Vector& newTrialTemp)
{
    // check vector arg is of correct size
    if (newTrialTemp.Size() != numberDOF) {
		opserr << "WARNING HeatTransferNode::setTrialTemperature() - incompatable sizes\n";
		return -2;
		}    

    // construct memory and Vectors for trial and committed
    // temperature on first call to this method, getTrialTemperature(),
    // getTemperature(), or incrTrialTemperature()        
    if (trialTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::setTrialTemperature() - ran out of memory\n";
			exit(-1);
			}    
		}

    // perform the assignment .. we dont't go through Vector interface
    // as we are sure of size and this way is quicker
    for (int i = 0; i < numberDOF; i++) {
		double T_trial = newTrialTemp(i);
		T[i+2*numberDOF] = T_trial - T[i+numberDOF];  // incrTemp
		T[i+3*numberDOF] = T_trial - T[i];  //	incrDeltaTemp
		T[i] = T_trial;  // trialTemp
		}

    return 0;
}


int
HeatTransferNode::setTrialTdot(const Vector& newTrialTot)
{
    // check vector arg is of correct size
    if (newTrialTot.Size() != numberDOF) {
		opserr << "WARNING HeatTransferNode::setTrialTdot() - incompatable sizes\n";
		return -2;
		}    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialTdot(),
    // getTdot(), or incrTrialTdot()        
    if (trialTdot == 0) {
		if (this->createTdot() < 0) {
			opserr << "FATAL HeatTransferNode::setTrialTdot() - ran out of memory\n";
			exit(-1);
			}
		}      
    
    // set the trial quantities
    for (int i = 0; i < numberDOF; i++)
		Tdot[i] = newTrialTot(i);

    return 0;
}


int
HeatTransferNode::setResponse(double value, int dof_num)
{
    if (dof_num < 0 || dof_num >= numberDOF){
		opserr << "FATAL HeatTransferNode::setInitialBC()" 
			   << " - improper dof number assigned\n";
		return -1;
		}
	
	if (trialTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::setInitialBC() -- ran out of memory\n";
			exit(-1);
			}
		}  

	// assign the value to array holding committed and trial values
	T[dof_num] = value;  
	T[dof_num+numberDOF] = value;

	return 0;
}


int
HeatTransferNode::incrTrialTemperature(const Vector& incrTemp)
{
    // check vector arg is of correct size
    if (incrTemp.Size() != numberDOF) {
		opserr << "WARNING HeatTransferNode::incrTrialTemperature() - incompatable sizes\n";
		return -2;
		}    

    // create a copy if no trial exists andd add committed
    if (trialTemp == 0) {
		if (this->createTemp() < 0) {
			opserr << "FATAL HeatTransferNode::incrTrialTemperature() - ran out of memory\n";
			exit(-1);
			}   

		for (int i = 0; i < numberDOF; i++) {
			double incrTempI = incrTemp(i);
			T[i] = incrTempI;
			T[i+2*numberDOF] = incrTempI;
			T[i+3*numberDOF] = incrTempI;
			}

		return 0;
		}

    // otherwise set trial = incr + trial
    for (int i = 0; i < numberDOF; i++) {
		double incrTempI = incrTemp(i);
		T[i] += incrTempI;
		T[i+2*numberDOF] += incrTempI;
		T[i+3*numberDOF] = incrTempI;
		}

    return 0;
}


int
HeatTransferNode::incrTrialTdot(const Vector& incrTdot)
{
    // check vector arg is of correct size
    if (incrTdot.Size() != numberDOF) {
	opserr << "WARNING HeatTransferNode::incrTrialTdot() - incompatable sizes\n";
	return -2;
    }    

    // create Vectors and array if none exist and set trial
    if (trialTdot == 0) {
		if (this->createTdot() < 0) {
			opserr << "FATAL HeatTransferNode::incrTrialTdot - ran out of memory\n";
			exit(-1);
			}    

		for (int i = 0; i < numberDOF; i++)
			Tdot[i] = incrTdot(i);

		return 0;
		}

    // otherwise set trial = incr + trial
    for (int i = 0; i < numberDOF; i++)
		Tdot[i] += incrTdot(i);    

    return 0;
}


int
HeatTransferNode::commitState()
{
    // check disp exists, if does set commit = trial, incr = 0.0
    if (trialTemp != 0) {
		for (int i=0; i < numberDOF; i++) {
			T[i+numberDOF] = T[i];  
			T[i+2*numberDOF] = 0.0;
			T[i+3*numberDOF] = 0.0;
			}
		}		    
    
    // check Tdot exists, if does set commit = trial    
    if (trialTdot != 0) {
		for (int i=0; i < numberDOF; i++)
			Tdot[i+numberDOF] = Tdot[i];
		}

    // if we get here we are done
    return 0;
}


int
HeatTransferNode::revertToLastCommit()
{
    // check T exists, if does set trial = last commit, incr = 0
    if (T != 0) {
		for (int i = 0 ; i < numberDOF; i++) {
			T[i] = T[i+numberDOF];
			T[i+2*numberDOF] = 0.0;
			T[i+3*numberDOF] = 0.0;
			}
		}
    
    // check Tdot exists, if does set trial = last commit
    if (Tdot != 0) {
		for (int i = 0 ; i < numberDOF; i++)
			Tdot[i] = Tdot[numberDOF+i];
		}

    // if we get here we are done
    return 0;
}


int
HeatTransferNode::revertToStart()
{
    // check T exists, if does set all to zero
    if (T != 0) {
		for (int i = 0 ; i < 4*numberDOF; i++)
			T[i] = 0.0;
		}

    // check Tdot exists, if does set all to zero
    if (Tdot != 0) {
		for (int i = 0 ; i < 2*numberDOF; i++)
			Tdot[i] = 0.0;
		}

    // if we get here we are done
    return 0;
}


// createDisp(), createVel() and createAccel():
// private methods to create the arrays to hold the disp, vel and acceleration
// values and the Vector objects for the committed and trial quantaties.

int
HeatTransferNode::createTemp(void)
{
  // trial , committed, incr = (committed-trial)
  T = new double[4*numberDOF];
    
  if (T == 0) {
    opserr << "WARNING - HeatTransferNode::createTemp() ran out of memory for array of size " 
		   << 2*numberDOF << endln;		    
    return -1;
  }

  for (int i=0; i < 4*numberDOF; i++)
	  T[i] = 0.0;
    
  commitTemp = new Vector(&T[numberDOF], numberDOF); 
  trialTemp = new Vector(T, numberDOF);
  incrTemp = new Vector(&T[2*numberDOF], numberDOF);
  incrDeltaTemp = new Vector(&T[3*numberDOF], numberDOF);
  
  if (commitTemp == 0 || trialTemp == 0 || incrTemp == 0 || incrDeltaTemp == 0) {
	  opserr << "WARNING - HeatTransferNode::createTemp() " 
		     << "ran out of memory creating Vectors(double *,int)";
	  return -2;
	  }
    
  return 0;
}


int
HeatTransferNode::createTdot(void)
{
    Tdot = new double[2*numberDOF];
    
    if (Tdot == 0) {
		opserr << "WARNING - HeatTransferNode::createTdot() ran out of memory for array of size " << 2*numberDOF << endln;
		return -1;
		}
    for (int i=0; i < 2*numberDOF; i++)
		Tdot[i] = 0.0;
    
    commitTdot = new Vector(&Tdot[numberDOF], numberDOF); 
    trialTdot = new Vector(Tdot, numberDOF);
    
    if (commitTdot == 0 || trialTdot == 0) {
		opserr << "WARNING - HeatTransferNode::createTdot() %s" 
			   << "ran out of memory creating Vectors(double *,int) \n";
      return -2;
    }
    
    return 0;
}
