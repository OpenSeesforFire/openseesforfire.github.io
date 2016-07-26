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
// Note: This class was adapted from IncrementalIntegrator  
#include <HeatTransferIntegrator.h>
#include <HT_FE_Element.h>
#include <LinearSOE.h>
#include <HT_AnalysisModel.h>
#include <Vector.h>
#include <HT_DOF_Group.h>
#include <HT_FE_EleIter.h>
#include <HT_DOF_GrpIter.h>

HeatTransferIntegrator::HeatTransferIntegrator()
:the_soe(0), the_model(0), the_test(0)
{

}


HeatTransferIntegrator::~HeatTransferIntegrator()
{

}


void
HeatTransferIntegrator::setLinks(HT_AnalysisModel& theModel, 
								 LinearSOE& theLinSOE, 
								 HT_ConvergenceTest* theConvergenceTest)
{
    the_model = &theModel;
    the_soe = &theLinSOE;
    the_test = theConvergenceTest;
}


int
HeatTransferIntegrator::domainChanged()
{
    return 0;
}


int 
HeatTransferIntegrator::formTangent()
{
    int result = 0;

    if (the_model == 0 || the_soe == 0) {
		opserr << "WARNING HeatTransferIntegrator::formTangent() -";
		opserr << " no HT_AnalysisModel or LinearSOE have been set\n";
		return -1;
		}

    // zero the A matrix of the linearSOE
    the_soe->zeroA();
    
    // the loops to form and add the tangents are broken into two for 
    // efficiency when performing parallel computations - CHANGE

    // loop through the FE_Elements adding their contributions to the tangent
    HT_FE_Element* elePtr;
    HT_FE_EleIter& theEles = the_model->getFEs();    
    while((elePtr = theEles()) != 0)     
		if (the_soe->addA(elePtr->getTangent(this),elePtr->getID()) < 0) {
			opserr << "WARNING HeatTransferIntegrator::formTangent -";
			opserr << " failed in addA for ID " << elePtr->getID();	    
			result = -3;
			}

    return result;
}


int 
HeatTransferIntegrator::formUnbalance(void)
{
    if (the_model == 0 || the_soe == 0) {
		opserr << "WARNING IncrementalIntegrator::formUnbalance -";
		opserr << " no HT_AnalysisModel or LinearSOE has been set\n";
		return -1;
		}
    
    the_soe->zeroB();
    
    if (this->formElementResidual() < 0) {
		opserr << "WARNING IncrementalIntegrator::formUnbalance ";
		opserr << " - this->formElementResidual failed\n";
		return -1;
		}  

    return 0;
}
    

int
HeatTransferIntegrator::newStep(double deltaT)
{
    return 0;
}


int
HeatTransferIntegrator::initialize(void)
{
    return 0;
}

int
HeatTransferIntegrator::commit(void) 
{
    if (the_model == 0) {
		opserr << "WARNING HeatTransferIntegrator::commit() -";
		opserr << "no HT_AnalysisModel object associated with this object\n";	
		return -1;
		}    

    return the_model->commitDomain();
}   


int
HeatTransferIntegrator::revertToLastStep(void) 
{
    return 0;
}   


LinearSOE*
HeatTransferIntegrator::getLinearSOE(void) const
{
    return the_soe;
}   

HT_ConvergenceTest*
HeatTransferIntegrator::getConvergenceTest(void) const
{
    return the_test;
}   

HT_AnalysisModel*
HeatTransferIntegrator::getModel(void) const
{
    return the_model;
}


int 
HeatTransferIntegrator::formElementResidual(void)
{
    // loop through the FE_Elements and add the residual
    HT_FE_Element* elePtr;

    int res = 0;    

    HT_FE_EleIter& theEles = the_model->getFEs();    
    while((elePtr = theEles()) != 0) {
		if (the_soe->addB(elePtr->getResidual(this),elePtr->getID()) < 0) {
			opserr << "WARNING HeatTransferIntegrator::formElementResidual -";
			opserr << " failed in addB for ID " << elePtr->getID();
			res = -2;
			}
		}
	//opserr << "deltaQ: " << the_soe->getB();
    return res;	    
}

