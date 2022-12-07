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

//Created: Liming Jiang@ PolyU (openseesforfire.git.io)
//Date: Jan 2020

//Note: This class was adapted from PathTimeSeries class in OpenSees

#ifndef UserDefinedFire_h
#define UserDefinedFire_h

#include <FireModel.h>
#include <PathTimeSeriesThermal.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
class Vector;



class UserDefinedFire : public FireModel
{
    public:
		UserDefinedFire(int tag);
		UserDefinedFire(int tag, const Vector& theData, const Vector& theTime,
			            int dataType);

		UserDefinedFire(int tag, const char* fileNameData, int dataType);    


		~UserDefinedFire();

		void applyFluxBC(HeatFluxBC* theFlux, double time);
		double getDuration();
		double getPeakData();

    protected:

    private:
		const Vector& getData(double time);
		//double getIncidentFlux(double time);

		Vector* theData;      // vector containg the data points
		Vector* time;		  // vector containg the time values of data points
		Vector* thePar;      // Vector containing the heat transfer parameter
		int currentTimeLoc;   // current location in time
        int type_tag;         // 1 for gas temperature, 2 for incident radiative flux
		static Vector thefireData;
};

#endif
