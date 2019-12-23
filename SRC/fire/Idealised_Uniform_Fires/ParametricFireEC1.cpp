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

#include <ParametricFireEC1.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <Convection.h>
#include <Radiation.h>


ParametricFireEC1::ParametricFireEC1(int tag, double I, double Av, double H, double At,
										   double Af, double Qf, double tlim, double startTime)
:FireModel(tag,2), ThI(I), Avent(Av), Hvent(H), Atotal(At), Afire(Af), Qfire(Qf), tlimit(tlim/3600), StartTime(startTime)
{
	//Qf in MJ
}

ParametricFireEC1::ParametricFireEC1(int tag)
:FireModel(tag, 2),ThI(-1.0), Avent(-1.0), Hvent(-1.0), Atotal(-1.0), Afire(-1.0), Qfire(-1.0), tlimit(-1.0), StartTime(0)
{
	
	
}


ParametricFireEC1::~ParametricFireEC1()
{

}


double 
ParametricFireEC1::getGasTemperature(double time)
{
    if ((Afire > 500.0) || (Hvent > 4.0))
		{ 
		opserr << "ParametricFireEC1::getGasTemperature() - floor area or compartment height is invalid.\n"
			<< "Note: Afire < 500m2, H < 4m2";
		exit(-1);
		}
	if ((ThI < 100.0) || (ThI > 2200.0)) 
		{
		opserr << "ParametricFireEC1::getGasTemperature() - thermal inertia value is invalid.\n"
			<< "Note: 100.0 <= ThermalI <=2200.0";
		exit(-1);
		}
	static double T, Ta, Ofactor, xx, x1, Gamma, tasterisk, Qtotal, tmax, tmax_ast;

	Ta = 20.0;

	if (time < StartTime)
		return Ta;

    double theTime = (time-StartTime) / 3600.0;
	

	
	Ofactor = Avent * sqrt(Hvent) / Atotal;
	xx = 1160/0.04;
	x1 = Ofactor / ThI;
	Gamma = x1 * x1 * xx * xx;
	tasterisk = theTime * Gamma;
	Qtotal = Qfire * Afire / Atotal;
	tmax = 0.2 * 1e-3 * Qtotal / Ofactor;

	if ((Ofactor < 0.02) || (Ofactor > 0.2))
		{ 
		opserr << "ParametricFireEC1::getGasTemperature() - floor area or compartment height is invalid.\n"
			<< "Note: Afire < 500m2, H < 4m2";
		exit(-1);
		}

	tmax_ast = tmax * Gamma;
	static double TempMax;

	if (tmax > tlimit) {
		if (theTime <= tmax) {
			T = Ta + 1325 * ( 1 - 0.324 * exp(-0.2 * tasterisk) - 0.204 * exp(-1.7 * tasterisk) 
				- 0.472 * exp(-19.0 * tasterisk));
			} else {
				TempMax = Ta + 1325 * ( 1 - 0.324 * exp(-0.2 * tmax_ast) - 0.204 * exp(-1.7 * tmax_ast) 
					- 0.472 * exp(-19.0 * tmax_ast));
				if (tmax_ast <= 0.5) {
					T = TempMax - 625 * (tasterisk - tmax_ast);
					} else if ((tmax_ast > 0.5) && (tmax_ast < 2.0)) {
						T = TempMax - 250 * (3.0 - tmax_ast) * (tasterisk - tmax_ast);
					} else if (tmax_ast >= 2.0) {
						T = TempMax - 250 * (tasterisk - tmax_ast);
						} 
					if (T < 20.0) T = 20.0;
			}

		return T;

		} else {
			double Olmt = 0.1 * 1e-3 * Qtotal / tlimit;
			double pt1 = Olmt / ThI;	
			double GammaLmt = pt1 * pt1 * xx * xx;
			if ((Ofactor > 0.04) && (Qtotal < 75.0) && (ThI < 1160.0)) {
				double k = 1 + (Ofactor / 0.04 - 1.0) * (Qtotal / 75.0 - 1.0) * (1 - ThI / 1160.0);
				GammaLmt *= k;
				}
			if (theTime <= tlimit) {
				tasterisk = theTime * GammaLmt;
				T = Ta + 1325 * ( 1 - 0.324 * exp(-0.2 * tasterisk) - 0.204 * exp(-1.7 * tasterisk) 
					- 0.472 * exp(-19.0 * tasterisk));
				return T;
				} else {
					double tlmt_ast = GammaLmt * tlimit;
					double tlmt_ast2 = Gamma * tlimit;
					TempMax = Ta + 1325 * ( 1 - 0.324 * exp(-0.2 * tlmt_ast) - 0.204 * exp(-1.7 * tlmt_ast) 
						- 0.472 * exp(-19.0 * tlmt_ast));
					if (tmax_ast <= 0.5) {
						T = TempMax - 625 * (tasterisk - tlmt_ast2);
						} else if ((tmax_ast > 0.5) && (tmax_ast < 2.0)) {
							T = TempMax - 250 * (3.0 - tmax_ast) * (tasterisk - tlmt_ast2);
						} else if (tmax_ast >= 2.0) {
							T = TempMax - 250 * (tasterisk - tlmt_ast2);
							} 
						if (T < 20.0) T = 20.0;
						return T;
				}
		}
}


void
ParametricFireEC1::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    // need to determine convection type or radiation
    int flux_type = theFlux->getTypeTag();
	if (flux_type == 1) 
		{
		Convection* convec = (Convection*) theFlux;
		convec->setSurroundingTemp(this->getGasTemperature(time)+273.15);
		double Temp;
		Temp= this->getGasTemperature(time);
		convec->applyFluxBC(time);
		} else if (flux_type == 2) {
			Radiation* rad = (Radiation*) theFlux;
			static const double bzm = 5.67 * 1e-008;
			//double alpha = rad->getAbsorptivity();
			double temp = this->getGasTemperature(time)+273.15;
			double qir = bzm * pow(temp, 4.0);
			rad->setIrradiation(qir);
			rad->applyFluxBC(time);
		} else {
			opserr << "ParametricFireEC1::applyFluxBC() - incorrect flux type provided.\n";
			exit(-1);
			}
}

