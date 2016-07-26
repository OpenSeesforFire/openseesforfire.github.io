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

#include <BackwardDifference.h>
#include <HT_FE_Element.h>
#include <LinearSOE.h>
#include <HT_AnalysisModel.h>
#include <Vector.h>
#include <HT_DOF_GrpIter.h>
#include <HT_DOF_Group.h>



BackwardDifference::BackwardDifference(bool vflag)
:Tt(0), Ttdot(0), T(0), Tdot(0), alpha(0.0),
 v_form(vflag)
{
    
}


BackwardDifference::~BackwardDifference()
{
    // clean up the memory created 
    if (Tt != 0)
		delete Tt;
    if (Ttdot != 0)
        delete Ttdot;
    if (T != 0)
        delete T;
    if (Tdot != 0)
        delete Tdot;
}

int
BackwardDifference::formEleTangent(HT_FE_Element* theEle)
{
    theEle->zeroTangent();
	if (v_form == true) { 
		theEle->addMcToTang();
		theEle->addMkAndMqToTang(alpha);
		} else {
			theEle->addMcToTang(1/alpha);
			theEle->addMkAndMqToTang();
		}

    return 0;
}


int BackwardDifference::newStep(double deltaT)
{
    // get a pointer to the AnalysisModel
    HT_AnalysisModel* theModel = this->getModel();
	alpha = deltaT;
    
    if (T == 0)  {
		opserr << "BackwardDifference::newStep() - domainChange() failed or hasn't been called\n";
		return -3;	
		}
    
    // set response at t to be that at t+deltaT of previous step
    (*Tt) = *T;        
    (*Ttdot) = *Tdot;  

	// set the value of Tdot 0 at beginning of each timestep
	(*Tdot).Zero();  

	theModel->setTemp(*T);
	theModel->setTdot(*Tdot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
		opserr << "BackwardDifference::newStep() - failed to update the domain\n";
		return -4;
		}
    
    return 0;
}


int BackwardDifference::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next newStep
    if (T != 0)  {
		(*T) = *Tt;        
		(*Tdot) = *Ttdot;  
		}

    return 0;
}


int BackwardDifference::domainChanged()
{
    HT_AnalysisModel* myModel = this->getModel();
    LinearSOE* theLinSOE = this->getLinearSOE();
    const Vector& x = theLinSOE->getX();
    int size = x.Size();
     
    // create the new Vector objects
    if (Tt == 0 || Tt->Size() != size)  {
        // if sizes unmatch, delete the old
		if (Tt != 0)
			delete Tt;
		if (Ttdot != 0)
			delete Ttdot;
		if (T != 0)
			delete T;
		if (Tdot != 0)
			delete Tdot;
        
        // create the new
        Tt = new Vector(size);
        Ttdot = new Vector(size);
        T = new Vector(size);
        Tdot = new Vector(size);
        
        // check we obtained the new
        if (Tt == 0 || Tt->Size() != size ||
			Ttdot == 0 || Ttdot->Size() != size ||
			T == 0 || T->Size() != size ||
			Tdot == 0 || Tdot->Size() != size)  {
				opserr << "BackwardDifference::domainChanged - ran out of memory\n";
				// if sizes unmatch, delete the old
				if (Tt != 0)
					delete Tt;
				if (Ttdot != 0)
					delete Ttdot;
				if (T != 0)
					delete T;
				if (Tdot != 0)
					delete Tdot;

				Tt = 0; 
				Ttdot = 0;
				T = 0; 
				Tdot = 0;

				return -1;
			}
    }        
    
    // now go through and populate T, Tdot by iterating through
    // the HT_DOF_Groups and getting the last committed T&Tdot
    HT_DOF_GrpIter& theDOFs = myModel->getDOFs();
    HT_DOF_Group* dofPtr;
    while ((dofPtr = theDOFs()) != 0)  {
		const ID& id = dofPtr->getID();
		int idSize = id.Size();

		int i;
		const Vector& nodalTemp = dofPtr->getCommittedTemp();	
		for (i = 0; i < idSize; i++)  {
			int loc = id(i);
			if (loc >= 0)  {
				(*T)(loc) = nodalTemp(i);		
				}
			}

		const Vector& nodalTdot = dofPtr->getCommittedTdot();
		for (i = 0; i < idSize; i++)  {
			int loc = id(i);
			if (loc >= 0)  {
				(*Tdot)(loc) = nodalTdot(i);
				}
			}
		}    
    
    return 0;
}


int BackwardDifference::update(const Vector& deltaT)
{
    HT_AnalysisModel* theModel = this->getModel();
    if (theModel == 0)  {
		opserr << "WARNING BackwardDifference::update() - no HT_AnalysisModel set\n";
		return -1;
		}	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Tt == 0)  {
		opserr << "WARNING BackwardDifference::update() - domainChange() failed or not called\n";
		return -2;
		}	
    
    // check deltaU is of correct size
    if (deltaT.Size() != T->Size())  {
		opserr << "WARNING BackwardDifference::update() - Vectors of incompatible size ";
		opserr << " expecting " << T->Size() << " obtained " << deltaT.Size() << endln;
		return -3;
		}

	if (v_form == true) {
		// determine the response at t+deltaT
		*Tdot += deltaT;
		*T = *Tt + alpha * (*Tdot);  // may be changed for better efficiency
		} else {
			*T += deltaT;
			(*Tdot) = (*T - *Tt) / alpha;
		}
        
    // update the response at the DOFs
    theModel->setResponse(*T,*Tdot);
    if (theModel->updateDomain() < 0)  {
		opserr << "BackwardDifference::update() - failed to update the domain\n";
		return -4;
		}
    
    return 0;
}    
