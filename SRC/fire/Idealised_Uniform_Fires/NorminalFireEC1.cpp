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
// Modified by Liming Jiang (L.Jiang@ed.ac.uk)

#include <NorminalFireEC1.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <Convection.h>
#include <Radiation.h>


NorminalFireEC1::NorminalFireEC1(int tag, int typeTag,double startTime)
:FireModel(tag, 1),type_tag(typeTag),StartTime(startTime)
{

}


NorminalFireEC1::~NorminalFireEC1()
{
	
}

// default type 1:  standard temperature-time curve,
// type 2 : external fire curve, 
// type 3 : hydrocarbon curve
// type 4: ASTM E119
// type 5: exponent fire

double 
NorminalFireEC1::getGasTemperature(double time)
{
    static double T;
	static double Ta = 293.15;
	
	if (time < StartTime)
		return Ta;
	time = (time - StartTime) / 60.0;

	if (type_tag == 1) {
		static double a1 = 345.0;
		static double a2 = 8.0;
		static double a3 = 1.0;
		//T=600;
		T = Ta + a1 * log10(a2 * time + a3);
		return T;
	} else if (type_tag == 2) {
		static double b1 = 660.0;
		static double b2 = 1.0;
		static double b3 = -0.687;
		static double b4 = -0.32;
		static double b5 = - 0.313;
		static double b6 = -3.8;
		T = Ta + b1 * ( b2 + b3 * exp(b4 * time) + b5 * exp(b6 * time));
		return T;
	} else if (type_tag == 3) {
		static double d1 = 1080.0;
		static double d2 = 1.0;
		static double d3 = -0.325;
		static double d4 = -0.167;
		static double d5 = - 0.675;
		static double d6 = -2.5;
		T = Ta + d1 * ( d2 + d3 * exp(d4 * time) + d5 * exp(d6 * time));
		  return T;
	} else if (type_tag == 4) {
		//ASTM E119
		//T=600;
		T = Ta + 750*(1-exp(-3.79553*sqrt(time/60)))+170.41*sqrt(time/60);
		return T;
	}
	else if (type_tag == 5) {
		//Exponent fire
		//T=600;
		time = time * 60;
		double alpha = 0.005;
		double duration = 3600;  //60min
		double decay = 1400;
		if (time < duration)
			T = Ta + (800 - Ta)*(1 - exp(-alpha*(time)));
		else if (time < duration + decay)
			T = 800 - (800 - Ta)*(time - duration) / decay;
		else
			T = Ta;
		return T;
	} else {
			opserr << "NorminalFireEC1::NorminalFireEN1991(int ) - incorrect fire type provided\n";
		    exit(-1);
	}
}


void
NorminalFireEC1::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    // need to determine convection type or radiation
    int flux_type = theFlux->getTypeTag();
	if (flux_type == 1) {
		Convection* convec = (Convection*) theFlux;
		convec->setSurroundingTemp(this->getGasTemperature(time));
		convec->applyFluxBC(time);
	} else if (flux_type == 2) {
		Radiation* rad = (Radiation*) theFlux;
		 double bzm = 5.67 * 1e-008;
		//double alpha = rad->getAbsorptivity();
		double temp = this->getGasTemperature(time);
		double qir = bzm * pow(temp, 4.0);
		rad->setIrradiation(qir);
		rad->applyFluxBC(time);
	} else {
		opserr << "NorminalFireEC1::applyFluxBC() - incorrect flux type provided\n";
		exit(-1);
	}
}

