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
// Modified by Liming Jiang (Liming.jiang@ed.ac.uk)
//

#include <HTNodeRecorder.h>
#include <OPS_Globals.h>


#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HT_NodeIter.h>

#include <fstream>
#include <string>
#include <sstream>
#include <Vector.h>
using std::ios;


HTNodeRecorder::HTNodeRecorder(int tag)
:HTRecorder(tag),response(0),
 theNodalTags(0), theNodes(0),theOutputHandler(0),
 theDomain(0), initializationDone(false), initialRecording(false),numValidNodes(0)
{
    //lastRecorderTag++;
}

HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes, 
			           HeatTransferDomain &theDom)
:HTRecorder(tag),response(0),
 theNodalTags(0), theNodes(0), theDomain(&theDom), theOutputHandler(0),ThePathSeries(0),
 initializationDone(false),initialRecording(false), numValidNodes(0)
{

  // 
  // create memory to hold nodal ID's
  //
    if (nodes != 0) {
		int numNode = nodes->Size();
		if (numNode != 0) {
			theNodalTags = new ID(*nodes);
			if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
				opserr << "HTNodeRecorder::HTNodeRecorder - out of memory\n";
				}
			}
		} 

	//lastRecorderTag++;
}


//if using output handler to printout the recorded values
HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes, 
			           HeatTransferDomain &theDom,OPS_Stream &theOutputHandler)
:HTRecorder(tag),response(0), ThePathSeries(0),
 theNodalTags(0), theNodes(0), theDomain(&theDom), theOutputHandler(&theOutputHandler),
 initializationDone(false),initialRecording(false), numValidNodes(0)
{

  // 
  // create memory to hold nodal ID's
  //
    if (nodes != 0) {
		int numNode = nodes->Size();
		if (numNode != 0) {
			theNodalTags = new ID(*nodes);
			if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
				opserr << "HTNodeRecorder::HTNodeRecorder - out of memory\n";
				}
			}
		} 

	//lastRecorderTag++;
}

//if using internal data transition via pathtimeseries
HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes,
                               HeatTransferDomain &theDom,PathTimeSeriesThermal* thePathSeries)
:HTRecorder(tag),response(0),
theNodalTags(0), theNodes(0), theDomain(&theDom), theOutputHandler(0),ThePathSeries(thePathSeries),
initializationDone(false),initialRecording(false), numValidNodes(0)
{
  
  //
  // create memory to hold nodal ID's
  //
  if (nodes != 0) {
    int numNode = nodes->Size();
    if (numNode != 0) {
      theNodalTags = new ID(*nodes);
      if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
        opserr << "HTNodeRecorder::HTNodeRecorder - out of memory\n";
      }
    }
		} 
  
  //lastRecorderTag++;
}

//if using internal data transition via pathtimeseries and a file
HTNodeRecorder::HTNodeRecorder(int tag, const ID* nodes,
                               HeatTransferDomain &theDom,OPS_Stream &theOutputHandler,PathTimeSeriesThermal* thePathSeries)
:HTRecorder(tag),response(0), ThePathSeries(thePathSeries),
theNodalTags(0), theNodes(0), theDomain(&theDom), theOutputHandler(&theOutputHandler),
initializationDone(false),initialRecording(false), numValidNodes(0)
{
  
  //
  // create memory to hold nodal ID's
  //
  if (nodes != 0) {
    int numNode = nodes->Size();
    if (numNode != 0) {
      theNodalTags = new ID(*nodes);
      if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
        opserr << "HTNodeRecorder::HTNodeRecorder - out of memory\n";
      }
    }
		} 
  
  //lastRecorderTag++;
}

HTNodeRecorder::~HTNodeRecorder()
{
  if (theOutputHandler != 0) {
    theOutputHandler->endTag(); // Data
    delete theOutputHandler;
  }

  if (theNodalTags != 0)
    delete theNodalTags;

  if (theNodes != 0)
    delete [] theNodes;
}


int 
HTNodeRecorder::record(double timeStamp)
{
  int HandlerType =0;
	if (theOutputHandler != 0)
  {
    if (ThePathSeries!=0)
      HandlerType =3;
    else
		HandlerType = 1;
  }
	else
  {
    if (ThePathSeries!=0)
      HandlerType=2;
		else
      HandlerType = 0;
  }
	
  
  if (theDomain == 0) {
    opserr << "HTNodeRecorder::record() - failed to access the heat transfer domain\n";
    return -1;
  }
	
    
  if (initializationDone != true) {
    if (this->initialize() != 0) {
			opserr << "HTNodeRecorder::record() - failed in initialize()\n";
			return -2;
		}
  }
  
  if(HandlerType==1){
        	
      response.resize(numValidNodes+1);
			double t = theDomain->getCurrentTime();
			response(0)= t;

			for (int i=0; i<numValidNodes; i++) {
				HeatTransferNode* theNode = theNodes[i];
				//Vector crds = theNode->getCrds();
				//opserr<<crds;
				double T = (theNode->getTemperature())(0);
				response(i+1)=T-273.15;	
			}
	        theOutputHandler->write(response);	
  } //end of using outputHandler;
  else if(HandlerType == 2){
    
    Vector resData(numValidNodes);
   
    double t = theDomain->getCurrentTime();
    //response(0)= t;
    
    for (int i=0; i<numValidNodes; i++) {
      HeatTransferNode* theNode = theNodes[i];
      //Vector crds = theNode->getCrds();
      double T = (theNode->getTemperature())(0);
      resData(i)=T-293.15;
    }
    
    //PathTimeSeriesThermal::WriteResults(double currentTime, const Vector& newData)
    ThePathSeries->WriteResults(t,resData);
  }
  else if(HandlerType == 3){
    response.resize(numValidNodes+1);
    Vector resData(numValidNodes);
    
    double t = theDomain->getCurrentTime();
    response(0)= t;
    
    for (int i=0; i<numValidNodes; i++) {
      HeatTransferNode* theNode = theNodes[i];
      //Vector crds = theNode->getCrds();
      double T = (theNode->getTemperature())(0);
      resData(i)=T-293.15;
      response(i+1)=T-273.15;
    }
    
    //PathTimeSeriesThermal::WriteResults(double currentTime, const Vector& newData)
    ThePathSeries->WriteResults(t,resData);
    theOutputHandler->write(response);
  }
  else
  {
        	
			std::string recorder = "HTNodeRecorder";
		std::string recorderSuffix = ".dat";
		int Rtag = this->getTag();
		std::stringstream f;
		f << recorder << Rtag << recorderSuffix;
		std::string filename;
		filename = f.str();
		//static const char* fname = filename.c_str();
		const char* fname = filename.c_str();
		
		double t = theDomain->getCurrentTime();
		std::ofstream outpt;
		if(initialRecording){
				outpt.open(fname, ios::out);
				initialRecording =false;
		} 
		else {
			outpt.open(fname, ios::out | ios::app);
		}
		
		//outpt.open(fname);
		if (!outpt.is_open()) {
			opserr << "HTNodeRecorder::record(double timestamp) "
				<< "- could not open the file for output" << endln;
			exit(-1);
			}

		for (int i=0; i<numValidNodes; i++) {
			HeatTransferNode* theNode = theNodes[i];
			//Vector crds = theNode->getCrds();

			//opserr<<"NodeRecorder: "<<theNode->getTag()<<endln;
			double T = (theNode->getTemperature())(0);
			/*
			int dim = crds.Size();
			double x = crds(0);
			double y = crds(1);		
			if  (dim == 2) {
				outpt << t << "  " << x << "  " << y << "  " << (T-273.15) << endln;
			} else if (dim == 3) {
				//double z = crds(2);
				//outpt << t << "  " << x << "  " << y <<  "  " << z << "  " 
				outpt<< (T-273.15) ;
				if(i==numValidNodes-1)
					oupt<<endln;
			*/
				//Modified by Liming,2013	
				if(i==0){
					outpt << t <<"   "<< (T-273.15)<< " ";
				}
				else if (i == numValidNodes-1) {
					outpt << (T-273.15) << endln;
				}
				else {
					outpt<< (T-273.15)<< " ";
			
			
			}
		}//end of for loop;
		}//end of not using outputHandler;
  
		return 0;
}


int 
HTNodeRecorder::setDomain(HeatTransferDomain& theDom)
{
    theDomain = &theDom;
    return 0;
}


int
HTNodeRecorder::initialize(void)
{
    if (theDomain == 0) {
		opserr << "HTNodeRecorder::initialize() - either nodes or domain has not been set\n";
		return -1;
		}

	//
	// create & set nodal array pointer
	//

	if (theNodes != 0) 
		delete [] theNodes;

	numValidNodes = 0;

    // record values for nodes given in theNodalTags
	if (theNodalTags != 0) {
		int numNode = theNodalTags->Size();
		theNodes = new HeatTransferNode* [numNode];
		if (theNodes == 0) {
			opserr << "HTNodeRecorder::initialize() - out of memory\n";
			return -1;
			}

		for (int i=0; i<numNode; i++) {
			int nodeTag = (*theNodalTags)(i);
			HeatTransferNode* theNode = theDomain->getNode(nodeTag);
			if (theNode != 0) {
				theNodes[numValidNodes] = theNode;
				numValidNodes++;
				}
			}

		} else {          // record values for all the nodes in the domain
			int numNodes = theDomain->getNumNodes();
			theNodes = new HeatTransferNode* [numNodes];
			if (theNodes == 0) {
				opserr << "HTNodeRecorder::initialize() - out of memory\n";
				return -1;
				}
			HT_NodeIter& theDomainNodes = theDomain->getNodes();
			HeatTransferNode* theNode;
			numValidNodes = 0;
			while (((theNode = theDomainNodes()) != 0) && (numValidNodes < numNodes)) {
				theNodes[numValidNodes] = theNode;
				numValidNodes++;
				}
		}

	initializationDone = true;
	initialRecording =true;
	return 0;
}
