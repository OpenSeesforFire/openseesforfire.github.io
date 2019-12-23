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

#include <UserDefinedFire.h>
#include <Vector.h>
//#include <Channel.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <HeatFluxBC.h>
#include <Convection.h>
#include <Radiation.h>

using std::ios;
using std::ifstream;

UserDefinedFire::UserDefinedFire(int tag)
:FireModel(tag, 8),theData(0), time(0), currentTimeLoc(-1), type_tag(-1)
{
	
	
}
		   
UserDefinedFire::UserDefinedFire(int tag,const Vector& theFireData, const Vector& theTime, int typeTag)
:FireModel(tag ,8),theData(0), time(0), currentTimeLoc(0), type_tag(typeTag)
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


UserDefinedFire::UserDefinedFire(int tag, const char* theFireData, const char* theTime, int typeTag)		       
:FireModel(tag ,8),theData(0), time(0), currentTimeLoc(0), type_tag(typeTag)
{
    // determine the number of data points
    int numDataPoints1 =0;
	int numDataPoints2 =0;
	double dataPoint;
	ifstream theFile;

	// first open and go through file containg path
	theFile.open(theFireData, ios::in);
	if (theFile.bad() || !theFile.is_open()) {
		opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
		opserr << " - could not open file " << theFireData << endln;
		} else {
			while (theFile >> dataPoint)
				numDataPoints1++;
		}   
	theFile.close();

	// now open and go through file containg time
	ifstream theFile1;
	theFile1.open(theTime, ios::in);
	if (theFile1.bad() || !theFile1.is_open()) {
		opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
		opserr << " - could not open file " << theTime << endln;
		} else {
			while (theFile1 >> dataPoint)
				numDataPoints2++;
		}   
	theFile1.close();

	// check number of data entries in both are the same
	if (numDataPoints1 != numDataPoints2) {
		opserr << "WARNING UserDefinedFire::UserDefinedFire() - files containing data ";
		opserr << "points for path and time do not contain same number of points\n";
		} else {
			// create a vector and read in the data
			if (numDataPoints1 != 0) {
				// now create the two vector
				theData = new Vector(numDataPoints1);
				time = new Vector(numDataPoints1);

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

				// first open the file for temperature/flux and read in the data
				ifstream theFile2;
				theFile2.open(theFireData, ios::in);
				if (theFile2.bad() || !theFile2.is_open()) {
					opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
					opserr << " - could not open file " << theFireData << endln;
					delete theData;
					delete time;
					theData = 0;
					time =0;
					} else { // read in the path data and then do the time
						int count = 0;
						while (theFile2 >> dataPoint) {
							(*theData)(count) = dataPoint;
							count++;
							}

						// finally close the file
						theFile2.close();

						// now open the time file and read in the data
						ifstream theFile3;
						theFile3.open(theTime, ios::in);
						if (theFile3.bad() || !theFile3.is_open()) {
							opserr << "WARNING - UserDefinedFire::UserDefinedFire()";
							opserr << " - could not open file " << theTime << endln;
							delete theData;
							delete time;
							theData = 0;
							time =0;
							} else { // read in the data
								int count = 0;
								while (theFile3 >> dataPoint) {
									(*time)(count) = dataPoint;
									count++;
									}
								theFile3.close();
							} // read in the data 
					}   // read in the path data and then do the time
				}
		}
}


UserDefinedFire::~UserDefinedFire()
{
  if (theData != 0)
    delete theData;
  if (time != 0)
    delete time;
}


double
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
	return (value1 + (value2-value1)*(theTime-time1)/(time2 - time1));
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
    if (type_tag == 1)
		{
		// need to determine convection type or radiation
		if (flux_type == 1) {
			Convection* convec = (Convection*) theFlux;
			convec->setSurroundingTemp(this->getData(time)+273.15);
			} else if (flux_type == 2) {
				Radiation* rad = (Radiation*) theFlux;
				double bzm = rad->getBLZMConstant();
				double alpha = rad->getAbsorptivity();
				double temp = this->getData(time)+273.15;
				double qir = alpha * bzm * pow(temp, 4.0);
				rad->setIrradiation(qir);
				rad->applyFluxBC(time);
			} else {
				opserr << "UserDefinedFire::applyFluxBC() - incorrect flux type provided\n";
				}
		} else if (type_tag == 2){
			if (flux_type == 2) {
				Radiation* rad = (Radiation*) theFlux;
				double qir = this->getData(time);
				rad->setIrradiation(qir);
				rad->applyFluxBC(time);
				} else {
					opserr << "UserDefinedFire::applyFluxBC() - flux_type should be 2\n";
					exit(-1); 
				}
		} else {
			opserr << "UserDefinedFire::applyFluxBC() - incorrect input type provided\n";
			exit(-1); 
			}
}

