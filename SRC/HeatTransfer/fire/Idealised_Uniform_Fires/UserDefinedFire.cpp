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
// Modiefied by Liming Jiang (liming.jiang@polyu.edu.hk)

#include <UserDefinedFire.h>
#include <Vector.h>
//#include <Channel.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <HeatFluxBC.h>
#include <Convection.h>
#include <Radiation.h>
#include <PrescribedSurfFlux.h>

using std::ios;
using std::ifstream;

Vector UserDefinedFire:: thefireData(2);

UserDefinedFire::UserDefinedFire(int tag)
:FireModel(tag, 8),theData(0), time(0), currentTimeLoc(-1), type_tag(-1)
{
	
	
}
		   
UserDefinedFire::UserDefinedFire(int tag,const Vector& theFireData, const Vector& theTime, int typeTag)
:FireModel(tag ,8),theData(0), time(0), thePar(0), currentTimeLoc(0), type_tag(typeTag)
{
    // check vectors are of same size
    if (theFireData.Size() != theTime.Size()) {
		opserr << "WARNING UserDefinedFire::UserDefinedFire() - vector containing data ";
		opserr << "points for temperature/flux and time are not of the same size\n";
		} else {
			// create copies of the vectors
			theData = new Vector(theFireData);
			time = new Vector(theTime);

			// ensure did not run out of memory creating copies
			if (theData == 0 || theData->Size() == 0 ||
				time == 0 || time->Size() == 0) {
					opserr << "WARNING UserDefinedFire::UserDefinedFire() - out of memory\n ";
					if (theData != 0)
						delete theData;
					if (time != 0)
						delete time;
					theData = 0;
					time = 0;
				}
		}
}


UserDefinedFire::UserDefinedFire(int tag, const char* theFireData, int typeTag)		       
:FireModel(tag ,8),theData(0), time(0), currentTimeLoc(0), type_tag(typeTag)
{
	//Type1: time-gasT
	//Type2: time-HT
	//Type3: time-gasT-convec h
	//Type: time-ht-covnec h
    // determine the number of data points
    int numDataPoints =0;
	int numRows = 0;
	double dataPoint;
	ifstream theFile;

	// first open and go through file containg path
	theFile.open(theFireData, ios::in);
	if (theFile.bad() || !theFile.is_open()) {
		opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
		opserr << " - could not open file " << theFireData << endln;
		} else {
			while (theFile >> dataPoint)
				numDataPoints++;
		}

	theFile.close();

	if(type_tag==3|| type_tag == 4)
		numRows = numDataPoints / 3;
	else
		numRows = numDataPoints / 2;

	// check number of data entries in both are the same

	if (numRows != 0) {
		// now create the two vector
		theData = new Vector(numRows);
		time = new Vector(numRows);
		if (type_tag == 3 || type_tag == 4)
			thePar = new Vector(numRows);

		// ensure did not run out of memory creating copies
		if (theData == 0 || theData->Size() == 0 ||time == 0 || time->Size() == 0) {

			opserr << "WARNING UserDefinedFire::UserDefinedFire() - out of memory\n ";
			if (theData != 0)
				delete theData;
			if (time != 0)
				delete time;
			theData = 0;
			time = 0;
		}

		// first open the file for temperature/flux and read in the data
		
		theFile.open(theFireData, ios::in);
	  // read in the path data and then do the time
			int count = 0;
			while (theFile >> dataPoint) {
				(*time)(count) = dataPoint;
				theFile >> dataPoint;
				(*theData)(count) = dataPoint;
				

				if (type_tag == 3 || type_tag == 4)
				{
					theFile >> dataPoint;
					(*thePar)(count) = dataPoint;
					count++;
				}
				else
					count++;
			}
#ifdef _FDEBUG
			opserr << (*time)(10) << "," << (*theData)(10) << "," << (*thePar)(10)<<endln;
			opserr << (*time)(180) << "," << (*theData)(180) << "," << (*thePar)(180)<<endln;
#endif // DEBUG

			
			// finally close the file
			theFile.close();

	}
}



UserDefinedFire::~UserDefinedFire()
{
  if (theData != 0)
    delete theData;
  if (time != 0)
    delete time;
  //Mhd Anwar Orabi 2021: added the type_tag check to protect against deleting thePar when it was no initialised.
  if (type_tag == 3 || type_tag == 4) {
	  if (thePar != 0)
		  delete thePar;
  };
}


const Vector&
UserDefinedFire::getData(double theTime)
{

    // check for a quick return
    if (theData == 0)
		return 0.0;

	// determine indexes into the data array whose boundary holds the time
	double time1 = (*time)(currentTimeLoc);

	// check for another quick return
	if (theTime == time1)
		return (*theData)[currentTimeLoc];

	int size = time->Size();
	int sizem1 = size - 1;
	int sizem2 = size - 2;

	// check we are not at the end
	if (theTime > time1 && currentTimeLoc == sizem1)
		return 0.0;

	if (theTime < time1 && currentTimeLoc == 0)
		return 0.0;

	// otherwise go find the current interval
	double time2 = (*time)(currentTimeLoc+1);
	if (theTime > time2) {
		while ((theTime > time2) && (currentTimeLoc < sizem2)) {
			currentTimeLoc++;
			time1 = time2;
			time2 = (*time)(currentTimeLoc+1);
			}
		// if pseudo time greater than ending time reurn 0
		if (theTime > time2)
			return 0.0;

		} else if (theTime < time1) {
			while ((theTime < time1) && (currentTimeLoc > 0)) {
				currentTimeLoc--;
				time2 = time1;	
				time1 = (*time)(currentTimeLoc);
				}
			// if starting time less than initial starting time return 0
			if (theTime < time1)
				return 0.0;
		}

	double value1 = (*theData)[currentTimeLoc];
	double value2 = (*theData)[currentTimeLoc+1];
	thefireData(0)= (value1 + (value2-value1)*(theTime-time1)/(time2 - time1));
	if (type_tag == 3 || type_tag == 4) {
		double par1 = (*thePar)[currentTimeLoc];
		double par2 = (*thePar)[currentTimeLoc+1];
		thefireData(1) = (par1 + (par2 - par1) * (theTime - time1) / (time2 - time1));
	}
	return thefireData;
}


double
UserDefinedFire::getDuration()
{
    if (theData == 0)
		{
		opserr << "WARNING -- UserDefinedFire::getDuration() on empty Vector" << endln;
		return 0.0;
		}

	int lastIndex = time->Size(); // index to last entry in time vector
	return ((*time)[lastIndex-1]);
}


double
UserDefinedFire::getPeakData()
{
    if (theData == 0)
		{
		opserr << "WARNING -- UserDefinedFire::getPeakData() on empty Vector" << endln;
		return 0.0;
		}

	double peak = fabs((*theData)[0]);
	int num = theData->Size();
	double temp;

	for (int i = 1; i < num; i++)
		{
		temp = fabs((*theData)[i]);
		if (temp > peak)
			peak = temp;
		}

	return peak;
}



void
UserDefinedFire::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    int flux_type = theFlux->getTypeTag();
    // need to identify the type of data information
    if (type_tag == 1|| type_tag == 3)
		{
		// need to determine convection type or radiation
		if (flux_type == 1) {
			Convection* convec = (Convection*) theFlux;
			if (type_tag == 3 || type_tag == 4) {
				convec->setParameter(this->getData(time)(1));
			}
			convec->setSurroundingTemp(this->getData(time)(0)+273.15);
			convec->applyFluxBC(time);
		} else if (flux_type == 2) {
				Radiation* rad = (Radiation*) theFlux;
				double bzm = rad->getBLZMConstant();
				double alpha = rad->getAbsorptivity();
				double temp = this->getData(time)(0)+273.15;
				double qir = bzm * pow(temp, 4.0);
				rad->setIrradiation(qir);
				rad->applyFluxBC(time);
		} else {
				opserr << "UserDefinedFire::applyFluxBC() - incorrect flux type provided\n";
		}
	} 
	//Heat flux based fire model
	else if (type_tag == 2||type_tag == 4){
		if (flux_type == 2) {
			Radiation* rad = (Radiation*) theFlux;
			//double qir = this->getData(time)(0);
			//rad->setIrradiation(qir);
			rad->applyFluxBC(time);
			} 
		else if (flux_type == 1) {
			Convection* convec = (Convection*)theFlux;
			if (type_tag == 3 || type_tag == 4) {
				convec->setParameter(this->getData(time)(1));
			}
			convec->applyFluxBC(time);
		}
		else if (flux_type == 3) {
			PrescribedSurfFlux* pflux = (PrescribedSurfFlux*)theFlux;

			//int flux_type = pflux->getTypeTag();
			int eleTag = pflux->getElementTag();
			int fTag = pflux->getFaceTag();
			HeatTransferDomain* theDomain = pflux->getDomain();
			if (theDomain == 0) {
				opserr << "Idealised_Local_Fire::applyFluxBC() - HeatFluxBC has not been associated with a domain";
				exit(-1);
			}

			HeatTransferElement* theEle = theDomain->getElement(eleTag);
			if (theEle == 0) {
				opserr << "Idealised_Local_Fire::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
				exit(-1);
			}

			const ID& faceNodes = theEle->getNodesOnFace(fTag);
			int size = faceNodes.Size();
			Vector nodalFlux(size);

			for (int i = 0; i < size; i++) {
				int nodTag = faceNodes(i);
				HeatTransferNode* theNode = theDomain->getNode(nodTag);
				if (theNode == 0) {
					opserr << "Idealised_Local_Fire::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
					exit(-1);
				}
				nodalFlux(i) = this->getData(time)(0);
#ifdef _FDEBUG
				opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
#endif

				
			}

			pflux->setData(nodalFlux);
			pflux->applyFluxBC();
		}
		else {
				opserr << "UserDefinedFire::applyFluxBC() - flux_type should be 2\n";
				exit(-1); 
		}
	} else {
		opserr << "UserDefinedFire::applyFluxBC() - incorrect input type provided\n";
		exit(-1); 
	}
}

