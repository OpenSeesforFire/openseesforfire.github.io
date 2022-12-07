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
// Editted by Liming Jiang(liming.jiang@ed.ac.uk)
//
// Note: This class was adapted from Domain   
#include <HeatTransferDomain.h>
#include <stdlib.h>
#include <OPS_Globals.h>
#include <HT_ElementIter.h>
#include <HT_NodeIter.h>
#include <HeatFluxBCIter.h>
#include <HeatTransferElement.h>
#include <HeatTransferNode.h>
#include <MP_TemperatureBC.h>
#include <TemperatureBC.h>
#include <HeatFluxBC.h>
#include <BoundaryPattern.h>
#include <TempBCIter.h>
#include <MP_TempBCIter.h>
#include <BoundaryPatternIter.h>
#include <AllTempBCIter.h>
#include <HeatTransferAnalysis.h>
#include <Vertex.h>
#include <Graph.h>
#include <HTRecorder.h>
//#include <MeshRegion.h>
//#include <FE_Datastore.h>
//#include <FEM_ObjectBroker.h>
#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

//Domain       *ops_TheActiveDomain = 0;

HeatTransferDomain::HeatTransferDomain()
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), 
 currentGeoTag(0), hasDomainChangedFlag(false),
 theBounds(6)
{
    // init the arrays for storing the domain components
    theElements = new MapOfTaggedObjects();
    theNodes    = new MapOfTaggedObjects();
    theTempBCs  = new MapOfTaggedObjects();
	theMPTempBCs = new MapOfTaggedObjects();
    thePatterns = new MapOfTaggedObjects();

    // init the iters    
    theEleIter = new HT_ElementIter(theElements);    
    theNodIter = new HT_NodeIter(theNodes);
    theTempBC_Iter = new TempBCIter(theTempBCs);
	theMP_TempBCIter = new MP_TempBCIter(theMPTempBCs);
    thePatternIter = new BoundaryPatternIter(thePatterns);
    allTempBC_Iter = new AllTempBCIter(*this);
    
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || theTempBCs == 0 || theMPTempBCs == 0
		|| thePatterns == 0 || theEleIter == 0 || theNodIter == 0
		|| theTempBC_Iter == 0 || thePatternIter == 0 || theMP_TempBCIter == 0
		|| allTempBC_Iter == 0 ) {
			opserr << "HeatTransferDomain::HeatTransferDomain() - out of memory\n";
			exit(-1);
		}      
	
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0; 
}


HeatTransferDomain::~HeatTransferDomain()
{
    // delete the objects in the domain
    this->HeatTransferDomain::clearAll();

    // delete all the storage objects
    // SEGMENT FAULT WILL OCCUR IF THESE OBJECTS WERE NOT CONSTRUCTED
    // USING NEW
  
    if (theElements != 0)
		delete theElements;    
  
    if (theNodes != 0)
		delete theNodes;
  
    if (theTempBCs != 0)
		delete theTempBCs;

	if (theMPTempBCs != 0)
		delete theMPTempBCs;

    if (thePatterns != 0)
		delete thePatterns;
  
    if (theEleIter != 0)
		delete theEleIter;
  
    if (theNodIter != 0)
		delete theNodIter;
  
    if (theTempBC_Iter != 0)
		delete theTempBC_Iter;

	if (theMP_TempBCIter != 0)
		delete theMP_TempBCIter;
  
    if (allTempBC_Iter != 0)
		delete allTempBC_Iter;
  
    if (thePatternIter != 0)  // be careful, its counterpart is not deleted in OpenSees' Domain
		delete thePatternIter;

	for (int i=0; i<numRecorders; i++) 
		if (theRecorders[i] != 0)
			delete theRecorders[i];

	theRecorders = 0;
	numRecorders = 0;
}


// void addElement(Element *);
//	Method to add an element to the model.


int
HeatTransferDomain::addElement(HeatTransferElement* element)
{
    int eleTag = element->getTag();
    // check all the elements nodes exist in the domain
    const ID& nodes = element->getExternalNodes();
    int numDOF = 0;
    for (int i=0; i < nodes.Size(); i++) {
		int nodeTag = nodes(i);
		HeatTransferNode* nodePtr = this->getNode(nodeTag);
		if (nodePtr == 0) {
			// NOTE: NEED TO DO SOMETHING TO THE OVERLOADED << OPERATOR
			opserr << "WARNING HeatTransferDomain::addElement - In element " << *element;
			opserr << "\n no Node " << nodeTag << " exists in the HeatTransferDomain\n";
			return -1;
			}
		numDOF += nodePtr->getNumberDOF();
		}   

    // check if an Element with a similar tag already exists in the Domain
    TaggedObject* other = theElements->getComponentPtr(eleTag);
    if (other != 0) {
		opserr << "HeatTransferDomain::addElement - element with tag " << eleTag << "already exists in model\n"; 
		return -1;
		}

    // add the element to the container object for the elements
    bool result = theElements->addComponent(element);
    if (result == true) {
		element->setDomain(this);
		element->update();

    // finally check the ele has correct number of dof
#ifdef _G3DEBUG
    if (numDOF != element->getNumDOF()) { 
      
      opserr << "Domain::addElement - element " << eleTag << " - #DOF does not match with number at nodes\n";
      theElements->removeComponent(eleTag);
      return false;
    }
#endif      
	// mark the Domain as having been changed
	this->domainChange();
	} 
	else {
		opserr << "HeatTransferDomain::addElement - element " << eleTag << "could not be added to container\n";
	}
			     
	 
		return 0;
}


// void addNode(Node *);
//	Method to add a Node to the model.

int
HeatTransferDomain::addNode(HeatTransferNode* node)
{
    int nodTag = node->getTag();

    TaggedObject* other = theNodes->getComponentPtr(nodTag);
    if (other != 0) {
		opserr << "HeatTransferDomain::addNode - node with tag " << nodTag << "already exists in model\n"; 
		return false;
		}
  
    bool result = theNodes->addComponent(node);
    if (result == true) {
		node->setDomain(this);
		this->domainChange();

		// see if the physical bounds are changed
		// note this assumes 0,0,0,0,0,0 as startup min,max values
		const Vector &crds = node->getCrds();
		int dim = crds.Size();
		if (dim >= 1) {
			double x = crds(0);
			if (x < theBounds(0)) theBounds(0) = x;
			if (x > theBounds(3)) theBounds(3) = x;
			} 
		if (dim >= 2) {
			double y = crds(1);
			if (y < theBounds(1)) theBounds(1) = y;
			if (y > theBounds(4)) theBounds(4) = y;	  
			} 
		if (dim == 3) {
			double z = crds(2);
			if (z < theBounds(2)) theBounds(2) = z;
			if (z > theBounds(5)) theBounds(5) = z;	  
			}

		} 
	else {
		opserr << "HeatTransferDomain::addNode - node with tag " << nodTag << "could not be added to container\n";
		return -1;
	}

    return 0;
}


// void addSP_Constraint(SP_Constraint *);
//	Method to add a constraint to the model.
//

bool
HeatTransferDomain::addTemperatureBC(TemperatureBC* temp_bc)
{
#ifdef _G3DEBUG    
    // check the Node exists in the Domain
    int nodeTag = spConstraint->getNodeTag();
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as node node with tag" <<
	nodeTag << "does not exist in model\n";       	
      return false;
    }

    // check that the DOF specified exists at the Node
    int numDOF = nodePtr->getNumberDOF();
    if (numDOF < spConstraint->getDOF_Number()) {
	opserr << "Domain::addSP_Constraint - cannot add as node node with tag" << 
	  nodeTag << "does not have associated constrained DOF\n"; 
	return false;
    }      
#endif

    // check that no other object with similar tag exists in model
    int tag = temp_bc->getTag();
    TaggedObject* other = theTempBCs->getComponentPtr(tag);
    if (other != 0) {
		opserr << "HeatTransferDomain::addTemperatureBC - cannot add as TemperatureBC with tag" << 
		tag << "already exists in model\n";             
		return false;
		}
  
    bool result = theTempBCs->addComponent(temp_bc);
    if (result == false) {
		opserr << "HeatTransferDomain::addTemperatureBC - cannot add TemperatureBC with tag" << 
		tag << "to the container\n";             
		return false;
		} 

	temp_bc->setDomain(this);
	temp_bc->setNodalValue();
	this->domainChange(); 

    return true;
}


// void addMP_TemperatureBC(MP_Constraint *);
//	Method to add a constraint to the model.
//

bool
HeatTransferDomain::addMP_TemperatureBC(MP_TemperatureBC *mpTempBC)
{

  // check that no other object with similar tag exists in model
  int tag = mpTempBC->getTag();
  TaggedObject *other = theMPTempBCs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "Domain::addMP_TemperatureBC - cannot add as constraint with tag" <<
      tag << "already exists in model";             
			      
    return false;
  }
  
  bool result = theMPTempBCs->addComponent(mpTempBC);
  if (result == true) {
      mpTempBC->setDomain(this);
      this->domainChange();
  } else
    opserr << "Domain::addMP_Constraint - cannot add constraint with tag" << 
      tag << "to the container\n";                   
			      
  return result;
}


//
bool 
HeatTransferDomain::addBoundaryPattern(BoundaryPattern* the_pattern)
{
    // first check if a pattern with a similar tag exists in model
    int tag = the_pattern->getTag();
    TaggedObject* other = thePatterns->getComponentPtr(tag);
    if (other != 0) {
		opserr << "HeatTransferDomain::addBoundaryPattern - cannot add as BoundaryPattern with tag" <<
		tag << "already exists in model\n";             		
		return false;
		}    

    // now we add the pattern to the container
    bool result = thePatterns->addComponent(the_pattern);
    if (result == true) {
		the_pattern->setDomain(this);
		this->domainChange();
		}
    else 
		opserr << "HeatTransferDomain::addBoundaryPattern - cannot add BoundaryPattern with tag" <<
		tag << "to the container\n";                   	
			      
    return result;
}    


bool
HeatTransferDomain::addTemperatureBC(TemperatureBC* temp_bc, int pattern)
{
#ifdef _G3DEBUG
    // check the Node exists in the Domain
    int nodeTag = spConstraint->getNodeTag();
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as node with tag" <<
	nodeTag << "does not exist in model\n";
	return false;
    }

    // check that the DOF specified exists at the Node
    int numDOF = nodePtr->getNumberDOF();
    if (numDOF < spConstraint->getDOF_Number()) {
      opserr << "Domain::addSP_Constraint - cannot add as node with tag" <<
	nodeTag << "does not have associated constrained DOF\n"; 

	return false;
    }      
#endif

    // now add it to the pattern
    TaggedObject* thePattern = thePatterns->getComponentPtr(pattern);
    if (thePattern == 0) {
		opserr << "HeatTransferDomain::addTemperatureBC - cannot add as pattern with tag" <<
		pattern << "does not exist in domain\n"; 
		return false;
		}
    BoundaryPattern* the_pattern = (BoundaryPattern *)thePattern;
    bool result = the_pattern->addTemperatureBC(temp_bc);
    if (result == false) {
		opserr << "HeatTransferDomain::addTemperatureBC - " << pattern 
			   << "pattern could not add the addTemperatureBC\n";  
		return false;
		}

    temp_bc->setDomain(this);
    this->domainChange();  

    return true;
}


bool 
HeatTransferDomain::addHeatFluxBC(HeatFluxBC* flux, int pattern)
{
    // now add it to the pattern
    TaggedObject* thePattern = thePatterns->getComponentPtr(pattern);
    if (thePattern == 0) {
		opserr << "HeatTransferDomain::addHeatFluxBC() - no pattern with tag " << pattern 
			   << "exits in  the model, not adding the flux " << *flux << endln;
		return false;
		}

    BoundaryPattern* the_pattern = (BoundaryPattern *)thePattern;
    bool result = the_pattern->addHeatFluxBC(flux);
    if (result == false) {
		opserr << "HeatTransferDomain::addHeatFluxBC() - no pattern with tag" << 
		pattern << "in  the model, not adding the ele flux" << *flux << endln;
		return false;
		}

    this->domainChange();
    return result;
}


void
HeatTransferDomain::clearAll(void) {
    // clear the loads and constraints from any load pattern
    BoundaryPatternIter& the_patterns = this->getBoundaryPatterns();
    BoundaryPattern* thePattern;
    while ((thePattern = the_patterns()) != 0)
		thePattern->clearAll();

    // clean out the containers
    theElements->clearAll();
    theNodes->clearAll();
    theTempBCs->clearAll();
    thePatterns->clearAll();

	for (int i=0; i<numRecorders; i++)
		if (theRecorders[i] != 0)
			delete theRecorders[i];

	numRecorders = 0; 

	//if (theRecorders != 0) {
		//delete [] theRecorders;
		//theRecorders = 0;
		//}
  
    // set the time back to 0.0
    currentTime = 0.0;
    committedTime = 0.0;
    dT = 0.0;

    // set the bounds around the origin
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;        
  
    currentGeoTag = 0;
    // rest the flag to be as initial
    hasDomainChangedFlag = false;
}


HeatTransferElement*
HeatTransferDomain::removeElement(int tag)
{
    // remove the object from the container    
    TaggedObject* myEle = theElements->removeComponent(tag);
  
    // if not there return 0
    if (myEle == 0) 
		return 0;

    // otherwise mark the domain as having changed
    this->domainChange();
  
    // perform a downward cast to an Element (safe as only Element added to
    // this container, 0 the Elements DomainPtr and return the result of the cast  
    HeatTransferElement* result = (HeatTransferElement *)myEle;
    //  result->setDomain(0);
    return result;
}

HeatTransferNode*
HeatTransferDomain::removeNode(int tag)
{
    // remove the object from the container
    TaggedObject* myNode = theNodes->removeComponent(tag);
  
    // if not there return 0
    if (myNode == 0) 
		return 0;  

    // mark the domain has having changed 
    this->domainChange();
  
    // perform a downward cast to a Node (safe as only Node added to
    // this container and return the result of the cast
    HeatTransferNode* result = (HeatTransferNode *)myNode;
    // result->setDomain(0);
    return result;
}


TemperatureBC*
HeatTransferDomain::removeTemperatureBC(int tag)
{
    // remove the object from the container    
    TaggedObject* myTemp = theTempBCs->removeComponent(tag);
    
    // if not there return 0    
    if (myTemp == 0) 
	return 0;

    // mark the domain as having changed    
    this->domainChange();
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast    
    TemperatureBC* result = (TemperatureBC *)myTemp;
    // result->setDomain(0);

    // should check that theLoad and result are the same    
    return result;
}


BoundaryPattern*
HeatTransferDomain::removeBoundaryPattern(int tag)
{
    // remove the object from the container            
    TaggedObject* myPattern = thePatterns->removeComponent(tag);
    
    // if not there return 0    
    if (myPattern == 0)
	return 0;
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast            
    BoundaryPattern* result = (BoundaryPattern *)myPattern;
    // result->setDomain(0);

    //
    // now set the Domain pointer for all loads and SP constraints 
    // in the loadPattern to be 0
    //

    HeatFluxBC* myFluxes;
    HeatFluxBCIter& theFluxes = result->getHeatFluxBCs();
    while ((myFluxes = theFluxes()) != 0) {
      // theElementalLoad->setDomain(0);  // does nothing?
    }

    int num_tem_bcs = 0;
    TemperatureBC* myTemp;
    TemperatureBCIter& theTemps = result->getTemperatureBCs();
    while ((myTemp = theTemps()) != 0) {
	num_tem_bcs++;
	// theSP_Constraint->setDomain(0);
    }

    // mark the domain has having changed if numSPs > 0
    // as the constraint handlers have to be redone
    if (num_tem_bcs > 0)
      this->domainChange();

    // finally return the load pattern
    return result;    
}    


HeatFluxBC*
HeatTransferDomain::removeHeatFluxBC(int tag, int Pattern)
{
    // remove the object from the container            
    BoundaryPattern* thePattern = this->getBoundaryPattern(Pattern);
    
    // if not there return 0    
    if (thePattern == 0)
		return 0;
    
    return thePattern->removeHeatFluxBC(tag);
}    


TemperatureBC*
HeatTransferDomain::removeTemperatureBC(int tag, int Pattern)
{
    // remove the object from the container            
    BoundaryPattern* thePattern = this->getBoundaryPattern(Pattern);
    
    // if not there return 0    
    if (thePattern == 0)
		return 0;
    
    TemperatureBC* theTempBC = thePattern->removeTemperatureBC(tag);
    if (theTempBC != 0)
		this->domainChange();

    return theTempBC;
} 

//Providing method to remove mpTemperatureBC
MP_TemperatureBC*
HeatTransferDomain::removeMP_TemperatureBC(int tag)
{
    // remove the object from the container        
    TaggedObject *mc = theMPTempBCs->removeComponent(tag);
    
    // if not there return 0    
    if (mc == 0) 
	return 0;

    // mark the domain as having changed    
    this->domainChange();
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast        
    MP_TemperatureBC *result = (MP_TemperatureBC *)mc;
    // result->setDomain(0);
    return result;
}  


HT_ElementIter&
HeatTransferDomain::getElements()
{
    theEleIter->reset();    
    return *theEleIter;
}


HT_NodeIter&
HeatTransferDomain::getNodes()
{
    theNodIter->reset();    
    return *theNodIter;
}


TemperatureBCIter&
HeatTransferDomain::getTemperatureBCs()
{
    theTempBC_Iter->reset();
    return *theTempBC_Iter;
}


TemperatureBCIter&
HeatTransferDomain::getAllTempBCs() // Note: I changed the func name, referenced by PenaltyBC_Handler::handle()
{
    allTempBC_Iter->reset();
    return *allTempBC_Iter;
}

MP_TemperatureBCIter&
HeatTransferDomain::getMPTemperatureBCs()
{
	theMP_TempBCIter->reset();
	return *theMP_TempBCIter;
}


BoundaryPatternIter&
HeatTransferDomain::getBoundaryPatterns()
{
    thePatternIter->reset();
    return *thePatternIter;
}


HeatTransferElement*
HeatTransferDomain::getElement(int tag) 
{
    TaggedObject* myEle = theElements->getComponentPtr(tag);
    // if not there return 0 otherwise perform a cast and return that
    if (myEle == 0) 
		return 0;
    HeatTransferElement* result = (HeatTransferElement *)myEle;
    return result;
}


HeatTransferNode*
HeatTransferDomain::getNode(int tag) 
{
    TaggedObject* myNode = theNodes->getComponentPtr(tag);

    // if not there return 0 otherwise perform a cast and return that  
    if (myNode == 0) 
		return 0;  
    HeatTransferNode* result = (HeatTransferNode *)myNode;
    return result;
}


TemperatureBC*
HeatTransferDomain::getTemperatureBC(int tag) 
{
    TaggedObject* myTempBC = theTempBCs->getComponentPtr(tag);

    // if not there return 0 otherwise perform a cast and return that  
    if (myTempBC == 0) 
		return 0;
    TemperatureBC* result = (TemperatureBC *)myTempBC;
    return result;
}


MP_TemperatureBC*
HeatTransferDomain::getMPTemperatureBC(int tag)
{
   TaggedObject* myMPTempBC = theMPTempBCs->getComponentPtr(tag);
   // if not there return 0 otherwise perform a cast and return that
   if(myMPTempBC==0)
	   return 0;
   MP_TemperatureBC* result = (MP_TemperatureBC *)myMPTempBC;
   return result;

}


BoundaryPattern*
HeatTransferDomain::getBoundaryPattern(int tag) 
{
    TaggedObject* myPattern = thePatterns->getComponentPtr(tag);
    // if not there return 0 otherwise perform a cast and return that  
    if (myPattern == 0) 
		return 0;
    BoundaryPattern* result = (BoundaryPattern *)myPattern;
    return result;
}


double
HeatTransferDomain::getCurrentTime(void) const
{
    return currentTime;
}

int
HeatTransferDomain::getCommitTag(void) const
{
    return commitTag;
}


int 
HeatTransferDomain::getNumElements(void) const
{
    return theElements->getNumComponents();
}
int 
HeatTransferDomain::getNumNodes(void) const
{
    return theNodes->getNumComponents();
}

int 
HeatTransferDomain::getNumTempBCs(void) const
{
    return theTempBCs->getNumComponents();
}

int 
HeatTransferDomain::getNumMPTempBCs(void) const
{
	return theMPTempBCs->getNumComponents();
}

int 
HeatTransferDomain::getNumBoundaryPatterns(void) const
{
    return thePatterns->getNumComponents();
}


const Vector&
HeatTransferDomain::getPhysicalBounds(void)
{
    return theBounds;
}



void
HeatTransferDomain::setCommitTag(int newTag)
{
    commitTag = newTag;
}

void
HeatTransferDomain::setCurrentTime(double newTime)
{
    currentTime = newTime;
    dT = currentTime - committedTime;
}

void
HeatTransferDomain::setCommittedTime(double newTime)
{
    committedTime = newTime;
    dT = currentTime - committedTime;
}


void
HeatTransferDomain::applyBCs(double timeStep)
{
    currentTime = timeStep;
    dT = currentTime - committedTime; 
    //
    // first loop over nodes and elements getting them to first zero their loads
    //
    HeatTransferElement* elePtr;
    HT_ElementIter& theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0)
		elePtr->zeroFlux();    

    // now loop over load patterns, invoking applyLoad on them
    BoundaryPattern* thePattern;
    BoundaryPatternIter& thePatterns = this->getBoundaryPatterns();
    while((thePattern = thePatterns()) != 0)
		thePattern->applyBCs(timeStep);
    //
    // finally loop over TemperatureBCs of the domain
    //
    TemperatureBCIter& theTempBCs = this->getTemperatureBCs();
    TemperatureBC* tempBC;
    while ((tempBC = theTempBCs()) != 0) {
		tempBC->applyTemperatureBC(timeStep);
		}

	MP_TemperatureBCIter& theMPTempBCs = this->getMPTemperatureBCs();
    MP_TemperatureBC* MPtempBC;
    while ((MPtempBC = theMPTempBCs()) != 0)
	MPtempBC->applyConstraint(timeStep);
}


void
HeatTransferDomain::setFluxConstant(void)
{
    // loop over all the load patterns that are currently added to the domain
    // getting them to set their loads as now constant
    BoundaryPattern* thePattern;
    BoundaryPatternIter& thePatterns = this->getBoundaryPatterns();
    while((thePattern = thePatterns()) != 0)
		thePattern->setFluxConstant();
}


void
HeatTransferDomain::unsetFluxConstant(void)
{
    // loop over all the load patterns that are currently added to the domain
    // getting them to set their loads as now constant
    BoundaryPattern* thePattern;
    BoundaryPatternIter& thePatterns = this->getBoundaryPatterns();
    while((thePattern = thePatterns()) != 0)
		thePattern->unsetFluxConstant();
}


int
HeatTransferDomain::setInitial(double value, int dofTag)
{
    HeatTransferNode* nodePtr;
    HT_NodeIter& theNodeIter = this->getNodes();
	
	int ok = 0;

    while ((nodePtr = theNodeIter()) != 0) {
		ok += nodePtr->setResponse(value, dofTag);
		}

	if (ok != 0){
		opserr << "HeatTransferDomain::setInitial - error in setting initial conditions\n";
		return ok;
		}
	
//	this->update();  // ask elements to update their material model objects
  
	return ok;
}


int
HeatTransferDomain::setInitial(int nodeTag, double value, int dofTag)
{
    opserr << "HeatTransferDomain::setInitial(int nodeTag, int dofTag, double value)" 
           << " not implemented!" << endln;

    return -1;
}


int
HeatTransferDomain::commit(void)
{
    // 
    // first invoke commit on all nodes and elements in the domain
    //
    HeatTransferNode* nodePtr;
    HT_NodeIter& theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0) {
		nodePtr->commitState();
		}

    HeatTransferElement* elePtr;
    HT_ElementIter& theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0) {
		elePtr->commitState();
		}

    // set the new committed time in the domain
    committedTime = currentTime;
    dT = 0.0;

	// invoke record on all recorders
	for (int i=0; i<numRecorders; i++){
		if (theRecorders[i] != 0)
			theRecorders[i]->record(currentTime);
		}

    // update the commitTag
    commitTag++;
    return 0;
}


int
HeatTransferDomain::revertToLastCommit(void)
{
    // 
    // first invoke revertToLastCommit  on all nodes and elements in the domain
    //
    HeatTransferNode* nodePtr;
    HT_NodeIter& theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
		nodePtr->revertToLastCommit();

    HeatTransferElement* elePtr;
    HT_ElementIter& theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0) {
		elePtr->revertToLastCommit();
		}

    // set the current time and load factor in the domain to last committed
    currentTime = committedTime;
    dT = 0.0;

    // apply load for the last committed time
    this->applyBCs(currentTime);

    return this->update();
}


int
HeatTransferDomain::update(void)
{
    int ok = 0;

    // invoke update on all the ele's
    HT_ElementIter& theEles = this->getElements();
    HeatTransferElement* theEle;

    while ((theEle = theEles()) != 0) {
		ok += theEle->update();
		}

    if (ok != 0)
		opserr << "HeatTransferDomain::update - HeatTransferDomain failed in update\n";

    return ok;
}


int
HeatTransferDomain::update(double newTime, double dT)
{
    this->applyBCs(newTime);
    this->update();

    return 0;
}


int
HeatTransferDomain::analysisStep(double dT)
{
    return 0;
}


void
HeatTransferDomain::setDomainChangeStamp(int newStamp)
{
    currentGeoTag = newStamp;
}


void
HeatTransferDomain::domainChange(void)
{
    hasDomainChangedFlag = true;
}


bool 
HeatTransferDomain::getDomainChangeFlag(void)
{
    return hasDomainChangedFlag;
}


int
HeatTransferDomain::hasDomainChanged(void)
{	
    // if the flag indicating the domain has changed since the
    // last call to this method has changed, increment the integer
    // and reset the flag
    bool result = hasDomainChangedFlag;
    hasDomainChangedFlag = false;
    if (result == true) {
		currentGeoTag++;
		}

    // return the integer so user can determine if domain has changed 
    // since their last call to this method
    return currentGeoTag;
}


int
HeatTransferDomain::addRecorder(HTRecorder& theRecorder)
{
    if (theRecorder.setDomain(*this) != 0) {
		opserr << "HeatTransferDomain::addRecorder() - recorder could not be added\n";
		return -1;
		}

	for (int i=0; i<numRecorders; i++) {
		if (theRecorders[i] == 0) {
			theRecorders[i] = &theRecorder;
			return 0;
			}
		}

	HTRecorder** newRecorders = new HTRecorder*[numRecorders + 1]; 
	if (newRecorders == 0) {
		opserr << "HeatTransferDomain::addRecorder() - ran out of memory\n";
		return -1;
		}

	for (int i=0; i<numRecorders; i++)
		newRecorders[i] = theRecorders[i];
	newRecorders[numRecorders] = &theRecorder;

	if (theRecorders != 0)
		delete [] theRecorders;

	theRecorders = newRecorders;
	numRecorders++;
	return 0;
}


int
HeatTransferDomain::removeRecorders(void)
{
    for (int i=0; i<numRecorders; i++)  
      if (theRecorders[i] != 0)
	delete theRecorders[i];
    
    if (theRecorders != 0) {
      delete [] theRecorders;
    }
  
    theRecorders = 0;
    numRecorders = 0;
    return 0;
}

int
HeatTransferDomain::removeRecorder(int tag)
{
  for (int i=0; i<numRecorders; i++) {
    if (theRecorders[i] != 0) {
      if (theRecorders[i]->getTag() == tag) {
	delete theRecorders[i];
	theRecorders[i] = 0;
	return 0;
      }
    }    
  }
  
  return -1;
}


int
HeatTransferDomain::SelectingNodes(ID& NodesRange, int crdTag, double MinValue, double MaxValue, double Tolerance)
{
    //const ID& TempNodesRange = NodesRange; 
    int iniNumOfNodes = 0;
    bool iniSelecting = true;
    //check if it is selecting nodes in the whole domain
    if (NodesRange.Size() != 0)
    {
        iniNumOfNodes = NodesRange.Size();
        iniSelecting = false;
    }
    else
    {
        iniNumOfNodes = this->getNumNodes();
        iniSelecting = true;
    }


    vector<int> SelectedNodes;
    int NodeTag = 0;
    for (int i = 0; i < iniNumOfNodes; i++) {
        
            if (iniSelecting)
                NodeTag = i+1;
            else
                NodeTag = NodesRange(i);

            double NodalCrd = (this->getNode(NodeTag)->getCrds())(crdTag);
            if ((NodalCrd <= MaxValue + Tolerance) && (NodalCrd >= MinValue - Tolerance)) {
                SelectedNodes.push_back(NodeTag);
            }
       

    }
    int NewIDsize = SelectedNodes.size(); 
    //Mhd Anwar Orabi 2021
    if (NewIDsize == 0)
        NodesRange = 0;
    else 
        NodesRange.resize(NewIDsize);
    
    
    
    for (int i = 0; i < NewIDsize; i++) {
        NodesRange(i) = SelectedNodes[i];
    }
#ifdef _DEBUG
    opserr << "Debug mode node range currently selected is:" << NodesRange << endln;
#endif
    return 0;
}