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
//
// Note: This class was adapted from DOF_Group 
#include <HT_DOF_Group.h>
#include <stdlib.h>

#include <HeatTransferNode.h>
#include <Vector.h>
#include <Matrix.h>
#include <HeatTransferIntegrator.h>

#define MAX_NUM_DOF 256

Vector HT_DOF_Group::errVect(0);
Vector** HT_DOF_Group::theVectors;
int HT_DOF_Group::numDOFs(0);           // number of objects

//  DOF_Group(Node *);
//	construictor that take the corresponding model node.

HT_DOF_Group::HT_DOF_Group(int tag, HeatTransferNode* node)
:TaggedObject(tag), storage(0), myNode(node), 
 myID(node->getNumberDOF()), 
 numDOF(node->getNumberDOF())
{
    // get number of DOF & verify valid
    int numDOF = node->getNumberDOF();
    if (numDOF <= 0) {
		opserr << "HT_DOF_Group::HT_DOF_Group(Node *) ";
		opserr << " node must have at least 1 dof " << *node;
		exit(-1);
		}	

    // check the ID created is of appropriate size
    if (myID.Size() != numDOF) {
		opserr << "HT_DOF_Group::HT_DOF_Group(Node *) ";
		opserr << " ran out of memory creating ID for node " << *node;
		exit(-1);
		}

    // initially set all the IDs to be -2
    for (int i = 0; i < numDOF; i++)
		myID(i) = -2;

    if (numDOFs == 0) {
		theVectors  = new Vector*[MAX_NUM_DOF+1];

		if (theVectors == 0) {
			opserr << "HT_DOF_Group::DOF_Group(Node *) ";
			opserr << " ran out of memory";	    
			}
		for (int i=0; i < MAX_NUM_DOF; i++) {
			theVectors[i] = 0;
			}
		} 

    // create the storage vector
    if (numDOF <= MAX_NUM_DOF) {
        if (theVectors[numDOF] == 0) {
            theVectors[numDOF] = new Vector(numDOF);
            storage = theVectors[numDOF];
            if (storage == 0 || storage->Size() != numDOF) {  
                    opserr << "HT_DOF_Group::HT_DOF_Group( ) ";
                    opserr << "ran out of memory for vector/Matrix of size :";
                    opserr << numDOF << endln;
                    exit(-1);
                }
            } else {
                storage = theVectors[numDOF];
            }
        } else {
            // create matrices and vectors for each object instance
            storage = new Vector(numDOF);
            if (storage == 0 || storage->Size() ==0) {
                    opserr << "HT_DOF_Group::HT_DOF_Group( ) ";
                    opserr << " ran out of memory for Vector of size :";
                    opserr << numDOF << endln;
                    exit(-1);
                }
        }

    numDOFs++;
}


HT_DOF_Group::~HT_DOF_Group()
{
    numDOFs--;

    // set the pointer in the associated Node to 0, to stop
    // segmentation fault if node tries to use this object after destroyed
    if (myNode != 0) 
		myNode->setDOF_GroupPtr(0);

}    


void
HT_DOF_Group::setID(int index, int value)
{
    if ((index >= 0) && (index < numDOF))
		myID(index) = value;
    else {
		opserr << "WARNING HT_DOF_Group::setID - invalid location ";
		opserr << index << " in ID of size " << numDOF << endln;
    }	
}

// void setID(const ID &);
//	Method to set the ID to be same as that given.

void
HT_DOF_Group::setID(const ID& copy)
{
    myID = copy;
}
 

// const ID &getID(void) const;
//	Method to return the current ID.

const ID&
HT_DOF_Group::getID(void) const
{
    return myID;
}



int
HT_DOF_Group::doneID(void)
{
    return 0;
}



int
HT_DOF_Group::getNumDOF(void) const
{
    return numDOF;
}


int
HT_DOF_Group::getNodeTag(void) const
{
    if (myNode != 0)
    	return myNode->getTag();
    else
    	return -1;
}


int
HT_DOF_Group::getNumFreeDOF(void) const
{
    int numFreeDOF = numDOF;
    for (int i = 0; i < numDOF; i++)
		if (myID(i) == -1 || myID(i) == -4)  // make sure the value is consistent with what
			                                 // being set in TemperatureBCHandler class
	    numFreeDOF--;
    
    return numFreeDOF;
}


int
HT_DOF_Group::getNumConstrainedDOF(void) const
{   
    int numConstr = 0;
    for (int i = 0; i < numDOF; i++)
		if (myID(i) < 0)
			numConstr++;    

    return numConstr;
}    


const Vector&
HT_DOF_Group::getCommittedTemp(void)
{
    if (myNode == 0) {
		opserr << "HT_DOF_Group::getCommittedTemp: no associated Node ";
		opserr << " returning the erroneous value\n";
		return errVect;
		}
    return myNode->getTemperature();
}


const Vector&
HT_DOF_Group::getCommittedTdot(void)
{
    if (myNode == 0) {
	opserr << "HT_DOF_Group::getCommittedTdot: no associated Node ";
	opserr << " returning the erroneous value\n";
	return errVect;	
    }
    return myNode->getTdot();
}


void
HT_DOF_Group::setNodeTemp(const Vector& T)
{
    if (myNode == 0) {
		opserr << "HT_DOF_Group::setNodeDisp: no associated Node\n";
		return;
		}
    
    Vector& Temp = *storage;
	Temp = myNode->getTrialTemperature();
    
    for (int i = 0; i < numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
			Temp(i) = T(loc);  
		}

    myNode->setTrialTemperature(Temp);
}
	
	
// void setNodeVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
HT_DOF_Group::setNodeTot(const Vector& Temp_dot)
{
    if (myNode == 0) {
		opserr << "HT_DOF_Group::setNodeVel: 0 Node Pointer\n";
		return;
		}
    
    Vector& Tdot = *storage;
    Tdot = myNode->getTrialTdot();

    for (int i = 0; i < numDOF; i++) {
	int loc = myID(i);	    	
	if (loc >= 0) 
	    Tdot(i) = Temp_dot(loc);  
    }

    myNode->setTrialTdot(Tdot);
}


void
HT_DOF_Group::incrNodeTemp(const Vector& T_incr)
{
    if (myNode == 0) {
	opserr << "HT_DOF_Group::incrNodeTemp: 0 Node Pointer\n";
	exit(-1);
    }

    Vector& Temp = *storage;;

    if (Temp.Size() == 0) {
		opserr << "HT_DOF_Group::incrNodeTemp - out of space\n";
		return;
		}

    for (int i = 0; i < numDOF; i++) {
		int loc = myID(i);	    			
		if (loc >= 0)
			Temp(i) = T_incr(loc);  
		else Temp(i) = 0.0;  
		}

    myNode->incrTrialTemperature(Temp);
}


void
HT_DOF_Group::incrNodeTdot(const Vector& Tdot_incr)
{

    if (myNode == 0) {
		opserr << "HT_DOF_Group::setNodeVel: 0 Node Pointer\n";
		exit(-1);
		}
    
    Vector& Tdot = *storage;
 
    for (int i = 0; i < numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
			Tdot(i) = Tdot_incr(loc);
		else  Tdot(i) = 0.0;
		}

	myNode->incrTrialTdot(Tdot);
}


void
HT_DOF_Group::resetNodePtr(void)
{
    myNode = 0;
}
