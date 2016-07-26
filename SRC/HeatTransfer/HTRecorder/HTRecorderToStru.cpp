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
// Created by Liming Jiang (liming.jiang@ed.ac.uk)
//

#include <HTRecorderToStru.h>
#include <OPS_Globals.h>


#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HT_NodeIter.h>
#include <Vector.h>
#include <vector>
#include <stdlib.h>

#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
using std::ios;
using std::vector;


HTRecorderToStru::HTRecorderToStru(int tag)
:HTRecorder(tag),theNodes(0),theNodes1(0),
initialRecording(false), response(0), theOutputHandler(0),
 theDomain(0), initializationDone(false), numValidNodes(0), Tolerance(0)
{
    //lastRecorderTag++;
}


HTRecorderToStru::HTRecorderToStru(int tag, const Vector& theCrds, 
			           HeatTransferDomain &theDom, double tolerance)
:HTRecorder(tag), theNodes(0),theNodes1(0), theDomain(&theDom),theOutputHandler(0),
					 initialRecording(false),response(0),
 					initializationDone(false), numValidNodes(0), Tolerance(tolerance)
{
// create memory to hold nodal ID's
  // FInding the nodes for required Crds
  int nRows = theCrds.Size()-1;
  xCrd = theCrds(0);
  yCrds.resize(nRows);
  for(int i=0; i<nRows; i++){
    yCrds(i)= theCrds(i+1);
  }
#ifdef _DEBUG
  opserr<<"HTRecorderToStru:XcRD "<<xCrd<<endln;
#endif

    
	vector<int> xCrdNodes;

	ID RefNodeTagsL =ID(nRows);
	ID RefNodeTagsR = ID(nRows);
	ID RefNodeFound = ID(nRows);

	int SelectedNum =0;
	double xCrdReference =0;
		//Finding nodes with xCrd==xCrdReference
		
	xCrdReference = xCrd;
	xCrdNodes.clear();
	
	for(int i=0;i<theDomain->getNumNodes();i++){
		HeatTransferNode* theHTNode= theDomain->getNode(i+1);
		if((theHTNode->getCrds()(0)<= xCrdReference+Tolerance)&&(theHTNode->getCrds()(0)>= xCrdReference-Tolerance)){
			xCrdNodes.push_back(i);
		}
              
	}
	//Finding two reference nodes from xCrdNodes for each ycrd
	if (xCrdNodes.size()==0){
			opserr<< "HTRecorderToStru:No node is selected for xCrd: "<<xCrdReference<<endln;
		}

	for (int m=0;m<nRows;m++){	
		for(int j=1; j<xCrdNodes.size();j++) {
		//Need to check order if using other mesh tools
			HeatTransferNode* theNode = theDomain->getNode(xCrdNodes[j-1]);
			HeatTransferNode* theNode1 = theDomain->getNode(xCrdNodes[j]);
			double yLocTolerance = (theNode1->getCrds()(1))-(theNode->getCrds()(1));
			double yDif = yCrds(m)-(theNode->getCrds()(1));
			double yDif1 = (theNode1->getCrds()(1))-yCrds(m);
			if((yLocTolerance>=yDif)&&(yLocTolerance>=yDif1)){
				//opserr<<"ycrd: "<<yCrds(m)<<"NodeCrd: "<<theNode->getCrds()(1)<<"NodeCrd1: "<<theNode1->getCrds()(1)<<endln;
				RefNodeTagsL(m) = xCrdNodes[j-1];
				RefNodeTagsR(m)= xCrdNodes[j];
				RefNodeFound(m)=1;
				//if reference yCrd is larger than ycrd of node j-1 and less than node j;put j-1,j into reference NodeTagID
			}
		}
	}
			

	for(int i =0; i<nRows; i++){
		if (RefNodeFound(i)==0)
			opserr<<"HTRecorderToStru:fail to find a Node for y: "<<yCrds(i)<<endln;
	}

	theNodeTags = RefNodeTagsL;
	theNodeTags1 = RefNodeTagsR;
#ifdef _DEBUG
	opserr<<"theNodeTags "<<theNodeTags<<endln;
	opserr<<"theNodeTags1 "<<theNodeTags1<<endln;
#endif

}

//

HTRecorderToStru::HTRecorderToStru(int tag, const Vector& theCrds, 
			           HeatTransferDomain &theDom, OPS_Stream &theOutputHandler, double tolerance)
:HTRecorder(tag), theNodes(0),theNodes1(0), theDomain(&theDom),theOutputHandler(&theOutputHandler),
					 initialRecording(false),response(0),
 					initializationDone(false), numValidNodes(0), Tolerance(tolerance)
{

  // create memory to hold nodal ID's
  // FInding the nodes for required Crds
  int nRows = theCrds.Size()-1;
  xCrd = theCrds(0);
  yCrds.resize(nRows);
  for(int i=0; i<nRows; i++){
    yCrds(i)= theCrds(i+1);
  }
#ifdef _DEBUG
  opserr<<"HTRecorderToStru:XcRD "<<xCrd<<endln;
#endif

    
	vector<int> xCrdNodes;

	ID RefNodeTagsL =ID(nRows);
	ID RefNodeTagsR = ID(nRows);
	ID RefNodeFound = ID(nRows);

	int SelectedNum =0;
	double xCrdReference =0;
		//Finding nodes with xCrd==xCrdReference
		
	xCrdReference = xCrd;
	xCrdNodes.clear();
	
	for(int i=1;i<=theDomain->getNumNodes();i++){
		HeatTransferNode* theHTNode= theDomain->getNode(i);
		if((theHTNode->getCrds()(0)<= xCrdReference+Tolerance)&&(theHTNode->getCrds()(0)>= xCrdReference-Tolerance)){
			xCrdNodes.push_back(i);
		}
              
	}
	//Finding two reference nodes from xCrdNodes for each ycrd
	if (xCrdNodes.size()==0){
			opserr<< "HTRecorderToStru:No node is selected for xCrd: "<<xCrdReference<<endln;
		}

	for (int m=0;m<nRows;m++){	
		for(int j=1; j<xCrdNodes.size();j++) {
		//Need to check order if using other mesh tools
			HeatTransferNode* theNode = theDomain->getNode(xCrdNodes[j-1]);
			HeatTransferNode* theNode1 = theDomain->getNode(xCrdNodes[j]);
			double yLocTolerance = (theNode1->getCrds()(1))-(theNode->getCrds()(1));
			double yDif = yCrds(m)-(theNode->getCrds()(1));
			double yDif1 = (theNode1->getCrds()(1))-yCrds(m);
			if((yLocTolerance>=yDif)&&(yLocTolerance>=yDif1)){
				//opserr<<"ycrd: "<<yCrds(m)<<"NodeCrd: "<<theNode->getCrds()(1)<<"NodeCrd1: "<<theNode1->getCrds()(1)<<endln;
				RefNodeTagsL(m) = xCrdNodes[j-1];
				RefNodeTagsR(m)= xCrdNodes[j];
				RefNodeFound(m)=1;
				//if reference yCrd is larger than ycrd of node j-1 and less than node j;put j-1,j into reference NodeTagID
			}
		}
	}
			

	for(int i =0; i<nRows; i++){
		if (RefNodeFound(i)==0)
			opserr<<"HTRecorderToStru:fail to find a Node for y: "<<yCrds(i)<<endln;
	}

	theNodeTags = RefNodeTagsL;
	theNodeTags1 = RefNodeTagsR;
#ifdef _DEBUG
	opserr<<"theNodeTags "<<theNodeTags<<endln;
	opserr<<"theNodeTags1 "<<theNodeTags1<<endln;
#endif

}



HTRecorderToStru::~HTRecorderToStru()
{
	if (theOutputHandler != 0) {

	    theOutputHandler->endTag(); // Data

	    delete theOutputHandler;

	  } 
	if (theNodes != 0) 
		delete [] theNodes;
	if (theNodes1 != 0) 
		delete [] theNodes1;
 
}


int 
HTRecorderToStru::record(double timeStamp)
{
	bool usingHandler=true;

		if (theOutputHandler == 0)
			usingHandler = false;
		else 
			usingHandler = true;
		
	if (theDomain == 0) {
		return 0;
		}

	if (initializationDone == false) {
		if (this->initialize() != 0) {
			opserr << "HTRecorderToStru::record() - failed in initialize()\n";
			return -1;
		}
	}

	if(usingHandler){

        	

	      response.resize(numValidNodes+1);
				double t = theDomain->getCurrentTime();
				response(0)= t;

				for (int i=0; i<numValidNodes; i++) {
					HeatTransferNode* theNodeL = theNodes[i];
					HeatTransferNode* theNodeR = theNodes1[i];
					double yCrdL = theNodeL->getCrds()(1);
					double yCrdR = theNodeR->getCrds()(1);

					double TLeft = (theNodeL->getTemperature())(0);
					double TRight = (theNodeR->getTemperature())(0);
			        double T=0;
					if(yCrdR-yCrdL>Tolerance){
						T =TLeft*(yCrdR-yCrds(i))/(yCrdR-yCrdL)+TRight*(yCrds(i)-yCrdL)/(yCrdR-yCrdL);
					} 
					else {
						opserr << "HTRecorderToStru::record(double timestamp), "
							<<"NodeL and NodeR overlap at y cordinate "<<yCrdR<<endln;
					}
					response(i+1)=T-273.15;	
				}
		        theOutputHandler->write(response);	

	} //end of using outputHandler;
	else
	{
	//
	// now we go get the temperatures from nodes and save them in a file
	//
	std::string recorder = "HTRecorderToStru";
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
			opserr << "HTRecorderToStru::record(double timestamp) "
				<< "- could not open the file for output" << endln;
			exit(-1);
			}

	for (int i=0; i<numValidNodes; i++) {
		HeatTransferNode* theNodeL = theNodes[i];
		HeatTransferNode* theNodeR = theNodes1[i];
		double yCrdL = theNodeL->getCrds()(1);
		double yCrdR = theNodeR->getCrds()(1);

		double TLeft = (theNodeL->getTemperature())(0);
		double TRight = (theNodeR->getTemperature())(0);
        double T=0;
		if(yCrdR-yCrdL>Tolerance){
			T =TLeft*(yCrdR-yCrds(i))/(yCrdR-yCrdL)+TRight*(yCrds(i)-yCrdL)/(yCrdR-yCrdL);
		} 
		else {
			opserr << "HTRecorderToStru::record(double timestamp), "
				<<"NodeL and NodeR overlap at y cordinate "<<yCrdR<<endln;
		}

		if(i==0){
			outpt << t <<"   "<< (T-273.15)<< " ";
		}
		else if (i == numValidNodes-1) {
			outpt << (T-273.15) << endln;
		}
		else {
			outpt<< (T-273.15)<< " ";
		}
	} //end of for loop of nodes
} //end of using no outputHandler
	return 0;
}


int 
HTRecorderToStru::setDomain(HeatTransferDomain& theDom)
{
    theDomain = &theDom;
    return 0;
}


int
HTRecorderToStru::initialize(void)
{
    if (theDomain == 0) {
		opserr << "HTRecorderToStru::initialize() - either nodes or domain has not been set\n";
		return -1;
		}

	// create & set nodal array pointer

	if (theNodes != 0) 
		delete [] theNodes;
	if (theNodes1 != 0) 
		delete [] theNodes1;

	numValidNodes = 0;
    // record values for nodes given in theNodalTags
	if (theNodeTags != 0&&theNodeTags1 != 0) {
		int numNode = theNodeTags.Size();
		theNodes = new HeatTransferNode* [numNode];
		theNodes1 = new HeatTransferNode* [numNode];

		if (theNodes == 0||theNodes1 == 0) {
			opserr << "HTRecorderToStru::initialize() - out of memory\n";
			return -1;
		}

		for (int i=0; i<numNode; i++) {
			int nodeTagL = theNodeTags(i);
			int nodeTagR = theNodeTags1(i);
			HeatTransferNode* theNode = theDomain->getNode(nodeTagL);
			HeatTransferNode* theNode1 = theDomain->getNode(nodeTagR);
			if (theNode != 0&&theNode != 0) {
				theNodes[numValidNodes] = theNode;
				theNodes1[numValidNodes] = theNode1;
				numValidNodes++;
			}
		}
	} 
	if(numValidNodes!=yCrds.Size()|| numValidNodes!=yCrds.Size()){         
		opserr<<"Fail to initialize the HTRecorderToStru\n";
	}

	initializationDone = true;
	initialRecording =true;

	return 0;
}
