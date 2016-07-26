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
// Note: This class was adapted from AnalysisModel  
#include <stdlib.h>

#include <ArrayOfTaggedObjects.h>
#include <HT_AnalysisModel.h>
#include <HeatTransferDomain.h>
#include <HT_FE_Element.h>
#include <HT_DOF_Group.h>
#include <HT_DOF_GrpIter.h>
#include <HT_FE_EleIter.h>
#include <Graph.h>
#include <Vertex.h>
#include <HeatTransferNode.h>
#include <HT_NodeIter.h>
#include <TemperatureBCHandler.h>


#include <MapOfTaggedObjects.h>

#define START_EQN_NUM 0
#define START_VERTEX_NUM 0

//  AnalysisModel();
//	constructor

HT_AnalysisModel::HT_AnalysisModel()
:myDomain(0), myHandler(0),
 myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0)
{
    theFEs     = new ArrayOfTaggedObjects(256);
    theDOFs    = new ArrayOfTaggedObjects(256);
    theFEiter  = new HT_FE_EleIter(theFEs);
    theDOFiter = new HT_DOF_GrpIter(theDOFs);
} 


// ~AnalysisModel();    
HT_AnalysisModel::~HT_AnalysisModel()
{
    if (theFEs != 0) {
		theFEs->clearAll();
		delete theFEs;
		}

    if (theDOFs != 0) {
		theDOFs->clearAll();
		delete theDOFs;
		}

    if (theFEiter != 0)
		delete theFEiter;

    if (theDOFiter != 0)
		delete theDOFiter;

    if (myGroupGraph != 0) {
		delete myGroupGraph;    
		}	
  
    if (myDOFGraph != 0) {
		delete myDOFGraph;
		}
}    


void
HT_AnalysisModel::setLinks(HeatTransferDomain& theDomain, TemperatureBCHandler& theHandler)
{
    myDomain = &theDomain;
    myHandler = &theHandler;
}


bool
HT_AnalysisModel::addHT_FE_Element(HT_FE_Element* theFE_Element)
{
    // check we don't add a null pointer or this is a subclass
    // trying to use this method when it should'nt
    if (theFE_Element == 0 || theFEs == 0)
		return false;

    // check if an Element with a similar tag already exists in the Domain
    int tag = theFE_Element->getTag();
    TaggedObject* other = theFEs->getComponentPtr(tag);
    if (other != 0) {
		opserr << "HT_AnalysisModel::addHT_FE_Element - HT_FE_Element with tag " << tag << "already exists in model\n"; 
		return false;
		}

    // add the element to the container object for the elements
    bool result = theFEs->addComponent(theFE_Element);
    if (result == true) {
		theFE_Element->setModel(*this);
		numFE_Ele++;
		return true;  // o.k.
    } else
		return false;

    return result;
}


// void addDOF_Group(DOF_Group *);
//	Method to add an element to the model.

bool
HT_AnalysisModel::addHT_DOF_Group(HT_DOF_Group* theDOF_Group)
{

  // check we don't add a null pointer or this is a subclass trying
  // to use a method it should'nt be using
    if (theDOF_Group == 0 || theDOFs == 0)
		return false;
  

    // check if an Element with a similar tag already exists in the Domain
    int tag = theDOF_Group->getTag();
    TaggedObject* other = theDOFs->getComponentPtr(tag);
    if (other != 0) {
		opserr << "HT_AnalysisModel::addHT_DOF_Group - group with tag " << tag << "already exists in model\n"; 
		return false;
		}

    // add the element to the container object for the elements
    bool result = theDOFs->addComponent(theDOF_Group);
    if (result == true) {
		numDOF_Grp++;
		return true;  // o.k.
    } else
		return false;
}


void
HT_AnalysisModel::clearAll(void) 
{
    // if the graphs have been constructed delete them
    if (myDOFGraph != 0)
		delete myDOFGraph;

    if (myGroupGraph != 0)
		delete myGroupGraph;    

    theFEs->clearAll();
    theDOFs->clearAll();

    myDOFGraph = 0;
    myGroupGraph = 0;
    
    numFE_Ele =0;
    numDOF_Grp = 0;
    numEqn = 0;    
}


int
HT_AnalysisModel::getNumDOF_Groups(void) const
{
    return numDOF_Grp;
}


HT_DOF_Group*
HT_AnalysisModel::getDOF_GroupPtr(int tag)
{
    TaggedObject* other = theDOFs->getComponentPtr(tag);
    if (other == 0) {
		return 0;
		}
    HT_DOF_Group* result = (HT_DOF_Group *)other;
    return result;
}


HT_FE_EleIter&
HT_AnalysisModel::getFEs()
{
    theFEiter->reset();
    return *theFEiter;
}


HT_DOF_GrpIter&
HT_AnalysisModel::getDOFs()
{
    theDOFiter->reset();
    return *theDOFiter;
}


void 
HT_AnalysisModel::setNumEqn(int theNumEqn)
{
    numEqn = theNumEqn;
}


int 
HT_AnalysisModel::getNumEqn(void) const
{
    return numEqn;
}


Graph&
HT_AnalysisModel::getDOFGraph(void)
{
    if (myDOFGraph == 0) {
		int numVertex = this->getNumDOF_Groups();

		//  myDOFGraph = new Graph(numVertex);
		MapOfTaggedObjects* graphStorage = new MapOfTaggedObjects();
		myDOFGraph = new Graph(*graphStorage);

		//
		// create a vertex for each dof
		//

		HT_DOF_Group* dofPtr = 0;
		HT_DOF_GrpIter& theDOFs = this->getDOFs();
		while ((dofPtr = theDOFs()) != 0) {
			const ID& id = dofPtr->getID();
			int size = id.Size();
			for (int i = 0; i < size; i++) {
				int dofTag = id(i);
				if (dofTag >= START_EQN_NUM) {
					Vertex* vertexPtr = myDOFGraph->getVertexPtr(dofTag);
					if (vertexPtr == 0) {
						Vertex *vertexPtr = new Vertex(dofTag, dofTag);      
						if (vertexPtr == 0) {
							opserr << "WARNING HT_AnalysisModel::getDOFGraph";
							opserr << " - Not Enough Memory to create " << i+1 << "th Vertex\n";
							return *myDOFGraph;
							}
						if (myDOFGraph->addVertex(vertexPtr, false) == false) {
							opserr << "WARNING HT_AnalysisModel::getDOFGraph - error adding vertex\n";
							return *myDOFGraph;
							}
						}
					}
				}
			}

		// now add the edges, by looping over the FE_elements, getting their
		// IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM

		HT_FE_Element* elePtr =0;
		HT_FE_EleIter& eleIter = this->getFEs();
		int cnt = 0;

		while((elePtr = eleIter()) != 0) {
			const ID& id = elePtr->getID();
			cnt++;
			int size = id.Size();
			for (int i = 0; i < size; i++) {
				int eqn1 = id(i);

				// if eqnNum of DOF is a valid eqn number add an edge
				// to all other DOFs with valid eqn numbers.

				if (eqn1 >=START_EQN_NUM) {
					for (int j = i+1; j < size; j++) {
						int eqn2 = id(j);
						if (eqn2 >= START_EQN_NUM)
							myDOFGraph->addEdge(eqn1-START_EQN_NUM+START_VERTEX_NUM,
							eqn2-START_EQN_NUM+START_VERTEX_NUM);
						}
					}
				}
			}
		}    

  return *myDOFGraph;
}


Graph&
HT_AnalysisModel::getDOFGroupGraph(void)
{
  if (myGroupGraph == 0) {
    int numVertex = this->getNumDOF_Groups();

    if (numVertex == 0) {
		opserr << "WARNING HT_AnalysisModel::getGroupGraph";
		opserr << "  - 0 vertices, has the HeatTransferDomain been populated?\n";
		exit(-1);
		}	

    //    myGroupGraph = new Graph(numVertex);
    MapOfTaggedObjects* graphStorage = new MapOfTaggedObjects();
    myGroupGraph = new Graph(*graphStorage);

    if (numVertex == 0) {
		opserr << "WARNING HT_AnalysisModel::getGroupGraph";
		opserr << "  - out of memory\n";
		exit(-1);
		}	
	
    HT_DOF_Group* dofPtr;

    // now create the vertices with a reference equal to the DOF_Group number.
    // and a tag which ranges from 0 through numVertex-1

    HT_DOF_GrpIter& dofIter2 = this->getDOFs();
    int count = START_VERTEX_NUM;
    while ((dofPtr = dofIter2()) != 0) {
		int DOF_GroupTag = dofPtr->getTag();
		int DOF_GroupNodeTag = dofPtr->getNodeTag();
		int numDOF = dofPtr->getNumFreeDOF();
		Vertex* vertexPtr = new Vertex(DOF_GroupTag, DOF_GroupNodeTag, 0, numDOF);

		if (vertexPtr == 0) {
			opserr << "WARNING HT_AnalysisModel::getDOFGroupGraph";
			opserr << " - Not Enough Memory to create ";
			opserr << count << "th Vertex\n";
			return *myGroupGraph;
			}
		myGroupGraph->addVertex(vertexPtr);
		}

    // now add the edges, by looping over the Elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    HT_FE_Element* elePtr;
    HT_FE_EleIter& eleIter = this->getFEs();

    while((elePtr = eleIter()) != 0) {
		const ID& id = elePtr->getDOFtags();
		int size = id.Size();
		for (int i = 0; i < size; i++) {
			int dof1 = id(i);
			for (int j = 0; j < size; j++) 
				if (i != j) {
					int dof2 = id(j);
					myGroupGraph->addEdge(dof1,dof2);
					}
			}
		}
	  }

  return *myGroupGraph;
}


void 
HT_AnalysisModel::setResponse(const Vector& T,
			                  const Vector& Tdot)
{
    HT_DOF_GrpIter& theDOFGrps = this->getDOFs();
    HT_DOF_Group* dofPtr;

    while ((dofPtr = theDOFGrps()) != 0) {
		dofPtr->setNodeTemp(T);
		dofPtr->setNodeTot(Tdot);	
		}
}	
	
void 
HT_AnalysisModel::setTemp(const Vector& T)
{
    HT_DOF_GrpIter& theDOFGrps = this->getDOFs();
    HT_DOF_Group* dofPtr;

    while ((dofPtr = theDOFGrps()) != 0) 
		dofPtr->setNodeTemp(T);
}	
	
void 
HT_AnalysisModel::setTdot(const Vector& Tdot)
{
    HT_DOF_GrpIter& theDOFGrps = this->getDOFs();
    HT_DOF_Group* dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
		dofPtr->setNodeTot(Tdot);
}	
	

void 
HT_AnalysisModel::incrTemp(const Vector& T)
{
    HT_DOF_GrpIter& theDOFGrps = this->getDOFs();
    HT_DOF_Group* dofPtr;

    while ((dofPtr = theDOFGrps()) != 0) 
		dofPtr->incrNodeTemp(T);
}	


void 
HT_AnalysisModel::incrTdot(const Vector& Tdot)
{
    HT_DOF_GrpIter& theDOFGrps = this->getDOFs();
    HT_DOF_Group* dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
		dofPtr->incrNodeTdot(Tdot);
}		


int
HT_AnalysisModel::updateDomain(void)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::updateDomain. No HeatTransferDomain linked.\n";
		return -1;
		}

    // invoke the method
    int res = myDomain->update();

    return 0;
}


int
HT_AnalysisModel::updateDomain(double newTime, double dT)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::updateDomain. No HeatTransferDomain linked.\n";
		return -1;
		}

    // invoke the method
    myDomain->applyBCs(newTime);

    return 0;
}


int
HT_AnalysisModel::analysisStep(double dT)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::newStep. No HeatTransferDomain linked.\n";
		return -1;
		}

    // invoke the method
    return myDomain->analysisStep(dT);
}


int
HT_AnalysisModel::commitDomain(void)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::commitDomain. No HeatTransferDomain linked.\n";
		return -1;
		}

    // invoke the method
    if (myDomain->commit() < 0) {
		opserr << "WARNING: HT_AnalysisModel::commitDomain - HeatTransferDomain::commit() failed\n";
		return -2;
		}

    return 0;
}


int
HT_AnalysisModel::revertDomainToLastCommit(void)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::revertDomainToLastCommit.";
		opserr << " No HeatTransferDomain linked.\n";
		return -1;
		}

    // invoke the method
    if (myDomain->revertToLastCommit() < 0) {
		opserr << "WARNING: HT_AnalysisModel::revertDomainToLastCommit.";
		opserr << " HeatTransferDomain::revertToLastCommit() failed.\n";
		return -2;
		}	
    return 0;
}


double
HT_AnalysisModel::getCurrentDomainTime(void)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::getCurrentDomainTime.";
		opserr << " No HeatTransferDomain linked.\n";
		return 0.0;
		}

    // invoke the method
    return myDomain->getCurrentTime();
}


void
HT_AnalysisModel::setCurrentDomainTime(double newTime)
{
    if (myDomain == 0) {
		opserr << "WARNING: HT_AnalysisModel::getCurrentDomainTime.";
		opserr << " No HeatTransferDomain linked.\n";
		}

    // invoke the method
    myDomain->setCurrentTime(newTime);
}


HeatTransferDomain*
HT_AnalysisModel::getDomainPtr(void) const
{
    return myDomain;
}
