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
                                                                        
// $Revision: 1.8 $
// $Date: 2009-10-13 21:13:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/penalty/HT_PenaltyMP_FE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the HT_PenaltyMP_FE class interface.
//
// the interface:

#include <HT_PenaltyMP_FE.h>
#include <stdlib.h>

#include <HeatTransferElement.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HT_DOF_Group.h>
#include <HeatTransferIntegrator.h>
#include <HT_TransientIntegrator.h>
#include <HT_AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <HeatTransferNode.h>
#include <MP_TemperatureBC.h>
#include <HT_DOF_Group.h>

HT_PenaltyMP_FE::HT_PenaltyMP_FE(int tag, HeatTransferDomain &theDomain, 
			   MP_TemperatureBC &TheMP, double Alpha)
:HT_FE_Element(tag, 2,2),
 theMP(&TheMP), theConstrainedNode(0) , theRetainedNode(0),
 tang(0), resid(0), C(0), alpha(Alpha)
{
    
    int size;
    const ID &id1 = theMP->getConstrainedDOFs();
    size = id1.Size();
    const ID &id2 = theMP->getRetainedDOFs();    
    size += id2.Size();

    tang = new Matrix(2,2);
    resid = new Vector(2);
    C = new Matrix(1,2);

    if (tang == 0 || resid == 0 || C == 0 ||
	tang->noCols() != size || C->noCols() != size || 
	resid->Size() != size) {
	opserr << "FATAL HT_PenaltyMP_FE::HT_PenaltyMP_FE() - out of memory\n";
	exit(-1);
    }
	    
    theRetainedNode = theDomain.getNode(theMP->getNodeRetained());    
    theConstrainedNode = theDomain.getNode(theMP->getNodeConstrained());

    if (theRetainedNode == 0 || theConstrainedNode == 0) {
	opserr << "FATAL HT_PenaltyMP_FE::HT_PenaltyMP_FE() - Constrained or Retained";
	opserr << " Node does not exist in Domain\n";
	opserr << theMP->getNodeRetained() << " " << theMP->getNodeConstrained() << endln;
	exit(-1);
    }	


    // set up the dof groups tags
    HT_DOF_Group *dofGrpPtr = 0;
    dofGrpPtr = theRetainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(0) = dofGrpPtr->getTag();	    
    else 
	opserr << "WARNING HT_PenaltyMP_FE::HT_PenaltyMP_FE() - node no Group yet?\n"; 
    dofGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(1) = dofGrpPtr->getTag();	        
    else
	opserr << "WARNING HT_PenaltyMP_FE::HT_PenaltyMP_FE() - node no Group yet?\n"; 
    
    
    if (theMP->isTimeVarying() == false) {
	this->determineTangent();
	// we can free up the space taken by C as it is no longer needed
	if (C != 0)
	    delete C;
	C = 0;
    }
}

HT_PenaltyMP_FE::~HT_PenaltyMP_FE()
{
    if (tang != 0)
	delete tang;
    if (resid != 0)
	delete resid;
    if (C != 0)
	delete C;
}    

// void setID(int index, int value);
//	Method to set the correMPonding index of the ID to value.
int
HT_PenaltyMP_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrained DOFs, this is obtained from the DOF_Group
    // associated with the constrained node
    HT_DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) 
	{
		opserr << "WARNING HT_PenaltyMP_FE::setID(void)";
		opserr << " - no DOF_Group with Constrained Node\n";
		return -2;
    }    

    const ID &constrainedDOFs = theMP->getConstrainedDOFs();               //Constrained DOF
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    //DOF id of the constrained node
    
    
	int constrained = constrainedDOFs(0);
	if (constrained < 0 || 
	    constrained >= theConstrainedNode->getNumberDOF()) 
	{
	    
	    opserr << "WARNING HT_PenaltyMP_FE::setID(void) - unknown DOF ";
	    opserr << constrained << " at Node\n";
	    myID(0) = -1; // modify so nothing will be added to equations
	    result = -3;
	}    	
	else 
	{
	    if (constrained >= theConstrainedNodesID.Size()) 
		{
			opserr << "WARNING HT_PenaltyMP_FE::setID(void) - ";
			opserr << " Nodes DOF_Group too small\n";
			myID(0) = -1; // modify so nothing will be added to equations
			result = -4;
	    }
	    else
		myID(0) = theConstrainedNodesID(constrained);      //myID(0) = Constrained Nodes ID
	}
    
    
    // now determine the IDs for the retained dof's
    HT_DOF_Group *theRetainedNodesDOFs = theRetainedNode->getDOF_GroupPtr();
    if (theRetainedNodesDOFs == 0) 
	{
		opserr << "WARNING HT_PenaltyMP_FE::setID(void)";
		opserr << " - no DOF_Group with Retained Node\n";
		return -2;
    }    
    
    const ID &RetainedDOFs = theMP->getRetainedDOFs();
    const ID &theRetainedNodesID = theRetainedNodesDOFs->getID();    

	int retained = RetainedDOFs(0);
	if (retained < 0 || retained >= theRetainedNode->getNumberDOF()) 
	{
	    opserr << "WARNING HT_PenaltyMP_FE::setID(void) - unknown DOF ";
	    opserr << retained << " at Node\n";
	    myID(1) = -1; // modify so nothing will be added
	    result = -3;
	}    	
	else 
	{
	    if (retained >= theRetainedNodesID.Size()) 
		{
			opserr << "WARNING HT_PenaltyMP_FE::setID(void) - ";
			opserr << " Nodes DOF_Group too small\n";
			myID(1) = -1; // modify so nothing will be added 
			result = -4;
	    }
	    else
		myID(1) = theRetainedNodesID(retained);
	}
    

    myDOF_Groups(0) = theConstrainedNodesDOFs->getTag();
    myDOF_Groups(1) = theRetainedNodesDOFs->getTag();
   
    return result;
}

const Matrix &
HT_PenaltyMP_FE::getTangent(HeatTransferIntegrator *theNewIntegrator)
{
    if (theMP->isTimeVarying() == true)
	this->determineTangent();    
    return *tang;
}

const Vector &
HT_PenaltyMP_FE::getResidual(HeatTransferIntegrator *theNewIntegrator)
{
    // zero residual, CD = 0
    return *resid;
}



const Vector &
HT_PenaltyMP_FE::getTangForce(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getTangForce() - not yet implemented\n";
 return *resid;
}

const Vector &
HT_PenaltyMP_FE::getK_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getK_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
HT_PenaltyMP_FE::getKi_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getK_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
HT_PenaltyMP_FE::getC_Force(const Vector &disp, double fact)
{
 opserr << "WARNING PenaltySP_FE::getC_Force() - not yet implemented\n";
 return *resid;
}

const Vector &
HT_PenaltyMP_FE::getM_Force(const Vector &disp, double fact)
{
  // opserr << "WARNING PenaltySP_FE::getM_Force() - not yet implemented\n";
 return *resid;
}

void  
HT_PenaltyMP_FE::determineTangent(void)
{
    // first determine [C] = [-I [Ccr]]
    C->Zero();
    const Matrix &constraint = theMP->getConstraint();
    
    
	(*C)(0,0) = -1.0;
    
    (*C)(0,1) = constraint(0,0);       //would be 1.0
    

    // now form the tangent: [K] = alpha * [C]^t[C]
    // *(tang) = (*C)^(*C); 
    // *(tang) *= alpha;

	// THIS IS A WORKAROUND UNTIL WE GET addMatrixTransposeProduct() IN
	// THE Matrix CLASS OR UNROLL THIS COMPUTATION
	int rows = C->noRows();
	int cols = C->noCols();
	Matrix CT(2,1);
	const Matrix &Cref = *C;
	// Fill in the transpose of C
	for (int k = 0; k < cols; k++)
		CT(k,0) = Cref(0,k);
	// Compute alpha*(C^*C)
	tang->addMatrixProduct(0.0, CT, Cref, alpha);
	//opserr <<*tang;
	//opserr<<endln;

}


