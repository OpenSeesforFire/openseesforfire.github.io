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
// Note: This class was adapted from FE_Element
#include <HT_FE_Element.h>
#include <stdlib.h>

#include <HeatTransferElement.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HT_DOF_Group.h>
#include <HeatTransferIntegrator.h>
#include <HT_AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>

#define MAX_NUM_DOF 64

// static variables initialisation
Matrix HT_FE_Element::errMatrix(1,1);
Vector HT_FE_Element::errVector(1);
Matrix** HT_FE_Element::theMatrices; // pointers to class wide matrices
Vector** HT_FE_Element::theVectors;  // pointers to class widde vectors
int HT_FE_Element::numFEs(0);        // number of objects

//  FE_Element(Element *, Integrator *theIntegrator);
//	construictor that take the corresponding model element.
HT_FE_Element::HT_FE_Element(int tag, HeatTransferElement* ele)
:TaggedObject(tag),
myDOF_Groups((ele->getExternalNodes()).Size()), myID(ele->getNumDOF()), 
numDOF(ele->getNumDOF()), theModel(0), myEle(ele), 
theResidual(0), theTangent(0), theIntegrator(0)
{
    if (numDOF <= 0) {
		opserr << "HT_FE_Element::HT_FE_Element(HeatTransferElement* ) ";
		opserr << " element must have 1 dof " << *ele;  // need to find out how the << operator is overloaded
		exit(-1);                                       // to output the info for Node or Element
		}   

    // get elements domain & check it is valid
    HeatTransferDomain* theDomain = ele->getDomain();
    if (theDomain == 0) {
		opserr << "FATAL HT_FE_Element::HT_FE_Element() - element has no domain "<< *ele;
		exit(-1);
		}

    // keep a pointer to all DOF_Groups
    int numGroups = ele->getNumExternalNodes();
    const ID& nodes = ele->getExternalNodes();

    for (int i = 0; i < numGroups; i++) {
		HeatTransferNode* nodePtr = theDomain->getNode(nodes(i));
		if (nodePtr == 0) {
			opserr << "FATAL HT_FE_Element::HT_FE_Element() - HeatTransferNode: ";
			opserr <<  nodes(i) << "does not exist in the HeatTransferDomain\n";
			opserr << *ele;
			exit(-1);
			}

		HT_DOF_Group* dofGrpPtr = nodePtr->getDOF_GroupPtr();
		if (dofGrpPtr != 0) 
			myDOF_Groups(i) = dofGrpPtr->getTag();	
		else {
			opserr << "FATAL HT_FE_Element::HT_FE_Element() - HeatTransferNode: ";
			opserr <<  *nodePtr <<  " has no HT_DOF_Group associated with it\n";
			exit(-1);
			}
		}

    // if this is the first FE_Element we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numFEs == 0) {
		theMatrices = new Matrix *[MAX_NUM_DOF+1];
		theVectors  = new Vector *[MAX_NUM_DOF+1];

		if (theMatrices == 0 || theVectors == 0) {
			opserr << "HT_FE_Element::HT_FE_Element(HeatTransferElement* ) ";
			opserr << " ran out of memory";	    
			}
		for (int i = 0; i < MAX_NUM_DOF; i++) {
			theMatrices[i] = 0;
			theVectors[i] = 0;
			}
		}
	
	// if Elements are not subdomains, set up pointers to
	// objects to return tangent Matrix and residual Vector.

    if (numDOF <= MAX_NUM_DOF) {
        // use class wide objects
        if (theVectors[numDOF] == 0) {
            theVectors[numDOF] = new Vector(numDOF);
            theMatrices[numDOF] = new Matrix(numDOF,numDOF);
            theResidual = theVectors[numDOF];
            theTangent = theMatrices[numDOF];
            if (theResidual == 0 || theResidual->Size() != numDOF ||	
                theTangent == 0 || theTangent->noCols() != numDOF)	{  
                    opserr << "HT_FE_Element::HT_FE_Element(HeatTransferElement* ) ";
                    opserr << " ran out of memory for vector/Matrix of size :";
                    opserr << numDOF << endln;
                    exit(-1);
                }
            } else {
                theResidual = theVectors[numDOF];
                theTangent = theMatrices[numDOF];
            }
        } else {
            // create matrices and vectors for each object instance
            theResidual = new Vector(numDOF);
            theTangent = new Matrix(numDOF, numDOF);
            if (theResidual == 0 || theTangent ==0 ||
                theTangent ==0 || theTangent->noRows() ==0) {

                    opserr << "HT_FE_Element::HT_FE_Element(HeatTransferElement* ) ";
                    opserr << " ran out of memory for vector/Matrix of size :";
                    opserr << numDOF << endln;
                    exit(-1);
                }
        } 
 
    // increment number of FE_Elements by 1
    numFEs++;
}


HT_FE_Element::HT_FE_Element(int tag, int numDOF_Group, int ndof)
  :TaggedObject(tag),
   myDOF_Groups(numDOF_Group), myID(ndof), numDOF(ndof), theModel(0),
   myEle(0), theResidual(0), theTangent(0), theIntegrator(0)
{
    // this is for a subtype, the subtype must set the myDOF_Groups ID array
    numFEs++;

    // if this is the first FE_Element we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numFEs == 0) {
		theMatrices = new Matrix *[MAX_NUM_DOF+1];
		theVectors  = new Vector *[MAX_NUM_DOF+1];

		if (theMatrices == 0 || theVectors == 0) {
			opserr << "HT_FE_Element::HT_FE_Element(HeatTransferElement* ) ";
			opserr << " ran out of memory";	    
			}
		for (int i=0; i<MAX_NUM_DOF; i++) {
			theMatrices[i] = 0;
			theVectors[i] = 0;
			}
		}
}



// ~FE_Element();    
//	destructor.
HT_FE_Element::~HT_FE_Element()
{
    // decrement number of FE_Elements
    numFEs--;

    // delete tangent and residual if created specially
	if (numDOF > MAX_NUM_DOF) {
		if (theTangent != 0) delete theTangent;
		if (theResidual != 0) delete theResidual;
		}

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numFEs == 0) {
		for (int i=0; i<MAX_NUM_DOF; i++) {
			if (theVectors[i] != 0)
				delete theVectors[i];
			if (theMatrices[i] != 0)
				delete theMatrices[i];
			}	
		delete[] theMatrices;
		delete[] theVectors;
		}
}    


const ID&
HT_FE_Element::getDOFtags(void) const 
{
    return myDOF_Groups;
}


const ID&
HT_FE_Element::getID(void) const
{
    return myID;
}


void 
HT_FE_Element::setModel(HT_AnalysisModel& the_model)
{
    theModel = &the_model;
}


int
HT_FE_Element::setID(void)
{
    int current = 0;
    
    if (theModel == 0) {
		opserr << "WARNING HT_FE_Element::setID() - no HT_AnalysisModel set\n";
		return -1;
		}
    
    int numGrps = myDOF_Groups.Size();
    for (int i = 0; i < numGrps; i++) {
		int tag = myDOF_Groups(i);

		HT_DOF_Group* dofPtr = theModel->getDOF_GroupPtr(tag);
		if (dofPtr == 0) {
			opserr << "WARNING HT_FE_Element::setID: 0 HT_DOF_Group Pointer\n";
			return -2;
			}

		const ID& theDOFid = dofPtr->getID();

		for (int j = 0; j < theDOFid.Size(); j++)  
			if (current < numDOF)
				myID(current++) = theDOFid(j);
			else {
				opserr << "WARNING HT_FE_Element::setID() - numDOF and";
				opserr << " number of dof at the DOF_Groups\n";
				return -3;
				}		
		}

    return 0;
}


const Matrix&
HT_FE_Element::getTangent(HeatTransferIntegrator* theNewIntegrator)
{
    theIntegrator = theNewIntegrator;
    
    if (myEle == 0) {
		opserr << "FATAL HT_FE_Element::getTangent() - no HeatTransferElement* given ";
		opserr << "- subclasses must provide implementation - ";
		exit(-1);
		}

    if (theNewIntegrator != 0)
		theNewIntegrator->formEleTangent(this);	    	    

		return *theTangent;
}

const Vector&
HT_FE_Element::getResidual(HeatTransferIntegrator* theNewIntegrator)
{
    theIntegrator = theNewIntegrator;

    if (theIntegrator == 0)
		return *theResidual;

    if (myEle == 0) {
		opserr << "FATAL HT_FE_Element::getTangent() - no HeatTransferElement* given ";
		opserr << "- subclasses must provide implementation - ";
		exit(-1);
		}    

      theNewIntegrator->formEleResidual(this);
      return *theResidual;
}


void  
HT_FE_Element::zeroTangent(void)
{
    if (myEle != 0)
		theTangent->Zero();
}

void  
HT_FE_Element::addMkAndMqToTang(double fact)
{
    if (myEle != 0) {
		// check for a quick return	
		if (fact == 0.0) 
			return;
		else   
			theTangent->addMatrix(1.0, myEle->getConductionTangent(), fact);

		if (myEle->hasConvection() == true)
			theTangent->addMatrix(1.0, myEle->getConvectionTangent(), fact);
		if (myEle->hasRadiation() == true)
			theTangent->addMatrix(1.0, myEle->getRadiationTangent(), fact);
		}
}

void  
HT_FE_Element::addMcToTang(double fact)
{
    if (myEle != 0) {
		// check for a quick return	
		if (fact == 0.0) 
			return;
		else  	    
			theTangent->addMatrix(1.0, myEle->getCapacityTangent(),fact);   	    	    	
		}
}
  

void  
HT_FE_Element::zeroResidual(void)
{
    if (myEle != 0) {
		theResidual->Zero();  	    	    
		}
	else {
		opserr << "WARNING HT_FE_Element::zeroResidual() - no HeatTransferElement* given ";
		opserr << "- subclasses must provide implementation - ";
		}    
}


void  
HT_FE_Element::addQkAndQqToResidual(double fact)
{
    if (myEle != 0) {
		// check for a quick return	
		if (fact == 0.0) {
			return;}
		else {
			const Vector& Qk = myEle->get_Q_Conduction();
			theResidual->addVector(1.0, Qk, fact);

			if ((myEle->hasConvection()) == true){
				const Vector& Q_qc = myEle->get_Q_Convection();
				theResidual->addVector(1.0, Q_qc, fact);
				}

			if(myEle->hasRadiation() == true){
				const Vector& Q_qr = myEle->get_Q_Radiation();
				theResidual->addVector(1.0, Q_qr, fact);
				}
			}
		}
    else {
		opserr << "WARNING HT_FE_Element::addQkAndQqToResidual() - no HeatTransferElement* given ";
		opserr << "- subclasses must provide implementation\n";
		}    	        
}


void  
HT_FE_Element::addQcToResidual(double fact)
{
    if (myEle != 0) {
		// check for a quick return	
		if (fact == 0.0) 
			return;
		else{
			const Vector& Qc = myEle->get_Q_Transient();
			theResidual->addVector(1.0, Qc, fact);
			} 	    	    	
		}
    else {
		opserr << "WARNING HT_FE_Element::addQcToResidual() - no HeatTransferElement* given ";
		opserr << "- subclasses must provide implementation\n";
		}    	        
}


void
HT_FE_Element::addQtotalToResidual(double fact)
{
    this->addQkAndQqToResidual(fact);
	this->addQcToResidual(fact);
}
