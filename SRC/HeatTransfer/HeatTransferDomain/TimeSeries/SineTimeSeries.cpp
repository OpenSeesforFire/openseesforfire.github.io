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

//Adapted by Yaqiang Jiang (y.jiang@ed.ac.uk)

#include <SineTimeSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

//#include <elementAPI.h>
//#define OPS_Export 


SineTimeSeries::SineTimeSeries(int tag, double startTime, double finishTime,
			                   double T, double peakAmplt)
:TimeSeries(tag, 1),
 tStart(startTime), tFinish(finishTime),
 period(T), peakAmplitude(peakAmplt)
{
	if (period == 0.0)  {
		opserr << "SineTimeSeries::SineTimeSeries -- input period is zero, setting period to 1\n";
		period = 1;
	}
}

SineTimeSeries::SineTimeSeries()
:TimeSeries(1),
 tStart(0.0), tFinish(0.0),
 period(1.0), peakAmplitude(1.0)
{
	// does nothing
}


SineTimeSeries::~SineTimeSeries()
{
	// does nothing
}

TimeSeries*
SineTimeSeries::getCopy(void) 
{
    return new SineTimeSeries(this->getTag(), tStart, tFinish, period, peakAmplitude);
}


double 
SineTimeSeries::getFactor(double time)
{	
	if (tStart <= time && time <= tFinish)  {
		double k = peakAmplitude * sin(2 * 3.1415926 * (time - tStart) / period);
		//if (k >= 0)
			return k;
		//else{
		//	opserr << "SineTimeSeries::getFactor -- got negative values, set to 0\n";
		//	return 0;
		//	}
		}
	else
		return 0;
}


int 
SineTimeSeries::sendSelf(int commitTag, Channel &theChannel)
{
    return -1;
}


int 
SineTimeSeries::recvSelf(int commitTag, Channel &theChannel, 
						  FEM_ObjectBroker &theBroker)
{
    return -1;
}


void 
SineTimeSeries::Print(OPS_Stream& s, int flag)
{
	s << "SineTimeSeries" << endln;
	s << "\tpeakAmplitude: " << peakAmplitude << endln;
	s << "\ttStart: " << tStart << endln;
	s << "\ttFinish: " << tFinish << endln;
	s << "\tPeriod: " << period << endln;
}
