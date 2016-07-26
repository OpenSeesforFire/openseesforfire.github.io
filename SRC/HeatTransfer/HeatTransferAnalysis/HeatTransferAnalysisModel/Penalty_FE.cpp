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
// Note: This class was adapted from PenaltySP_FE 
#include <Penalty_FE.h>
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
#include <TemperatureBC.h>
//#include <typeinfo>


Matrix Penalty_FE::tang(1,1);
Vector Penalty_FE::resid(1);

Penalty_FE::Penalty_FE(int tag, HeatTransferDomain& theDomain, 
			           TemperatureBC& theTempBC, double pnum)
:HT_FE_Element(tag, 1,1), penalty_num(pnum), 
 the_temp_bc(&theTempBC), the_node(0)
{
    // get a pointer to the Node
    the_node = theDomain.getNode(the_temp_bc->getNodeTag());
    if (the_node == 0) {
		opserr << "FATAL Penalty_FE::Penalty_FE() - no HeatTransferNode: ";
		opserr << the_temp_bc->getNodeTag() << "in HeatTransferDomain\n";
		exit(-1);
		}

    // set the DOF_Group tags
    HT_DOF_Group* dofGrpPtr = the_node->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
		myDOF_Groups(0) = dofGrpPtr->getTag();	    
}


Penalty_FE::~Penalty_FE()
{

}    

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.

int
Penalty_FE::setID(void)
{
    HT_DOF_Group* theNodesDOFs = the_node->getDOF_GroupPtr();
    if (theNodesDOFs == 0) {
		opserr << "WARNING Penalty_FE::setID(void) - no HT_DOF_Group with HeatTransferNode\n";
		return -2;
		}    
    myDOF_Groups(0) = theNodesDOFs->getTag();
    
    int restrainedDOF = the_temp_bc->getDOF_Number();
    if (restrainedDOF < 0 || restrainedDOF >= the_node->getNumberDOF()) {
		opserr << "WARNING PenaltySP_FE::setID(void) - unknown DOF ";
		opserr << restrainedDOF << " at HeatTransferNode\n";
		return -3;
		}    	

    const ID& theNodesID = theNodesDOFs->getID();
    if (restrainedDOF >= theNodesID.Size()) {
		opserr << "WARNING Penalty_FE::setID(void) - ";
		opserr << " ID in HT_DOF_Group too small\n";
		return -4;
		}    		
    
    myID(0) = theNodesID(restrainedDOF);

    return 0;
}


const Matrix&
Penalty_FE::getTangent(HeatTransferIntegrator* theNewIntegrator)
{
    tang(0,0) = penalty_num;
    return tang;
}


const Vector&
Penalty_FE::getResidual(HeatTransferIntegrator* theNewIntegrator)
{
    double constraint = the_temp_bc->getValue();
    int constrainedDOF = the_temp_bc->getDOF_Number();
    const Vector& nodeTemp = the_node->getTrialTemperature(); // getTrialTdot()
	
    if (constrainedDOF < 0 || constrainedDOF >= nodeTemp.Size()) {
		opserr << "WARNING Penalty_FE::getResidual() - ";	
		opserr << " constrained DOF " << constrainedDOF << " outside range\n";
		resid(0) = 0;
		}

    // (*resid)(0) = penalty_num * (constraint - nodeDisp(constrainedDOF));    
    // is replace with the following to remove possible problems with
    // subtracting very small numbers

	//bool vflag = theNewIntegrator->getUpdatingFlag();
 //   
	//// only transient integrator uses v-form updating
	//if (vflag == true) {
	//	HT_TransientIntegrator* theIntegrator = (HT_TransientIntegrator*)theNewIntegrator;
	//	double alpha = theIntegrator->getAlphaDeltat();
	//	resid(0) = penalty_num * (constraint - nodeTemp(constrainedDOF)) / alpha;
	//	} else {
	//		resid(0) = penalty_num * (constraint - nodeTemp(constrainedDOF));   
	//	}

    HT_TransientIntegrator* tranIntegrator = 0;
	//const char* name1 = "HT_TransientIntegrator";
	//const char* name2 = typeid(*theNewIntegrator).name();
	if (tranIntegrator = dynamic_cast<HT_TransientIntegrator*>(theNewIntegrator)) {
		//tranIntegrator = static_cast<HT_TransientIntegrator*>(theNewIntegrator);
		bool vflag = tranIntegrator->getUpdatingFlag();
		if (vflag == true) {
			double alpha = tranIntegrator->getAlphaDeltat();
			resid(0) = penalty_num * (constraint - nodeTemp(constrainedDOF)) / alpha;
			} else {
				resid(0) = penalty_num * (constraint - nodeTemp(constrainedDOF));   
			}
		} else {
			resid(0) = penalty_num * (constraint - nodeTemp(constrainedDOF));   
		}

    return resid;
}