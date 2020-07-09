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

#include <NWConcreteEC2.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>


NWConcreteEC2::NWConcreteEC2(int tag, double moisture)
:HeatTransferMaterial(tag), trial_temp(0.0), ini_temp(0.0), 
 rho_a(2300.0), rho(2300.0),cp(900.0), enthalpy(0.0), moist(moisture)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "CarbonSteelEC3::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


NWConcreteEC2::~NWConcreteEC2()
{
    if (k != 0)
		delete k;
}

int 
NWConcreteEC2::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
NWConcreteEC2::getConductivity(void)
{
    if (trial_temp <= 20.0) {
		(*k)(0,0) = 1.333028;
		(*k)(1,1) = 1.333028;
		(*k)(2,2) = 1.333028;
		} else if ((20.0 < trial_temp) && (trial_temp < 1200.0)) {
			double c = trial_temp / 100.0;
			double kc = 1.36 - 0.136 * c + 0.0057 * c * c;
			(*k)(0,0) = kc;
			(*k)(1,1) = kc;
			(*k)(2,2) = kc;
		} else {
			(*k)(0,0) = 0.5488;
			(*k)(1,1) = 0.5488;
			(*k)(2,2) = 0.5488;
		 }

    return *k;
}


double  
NWConcreteEC2::getRho(void)
{
	if (trial_temp <= 115)
		rho = rho_a;
	else if (trial_temp <= 200)
		rho = rho_a * (1 - 0.02 * (trial_temp - 115) / 85);
	else if (trial_temp <= 400)
		rho = rho_a * (0.98 - 0.03 * (trial_temp - 200) / 200);
	else if (trial_temp <= 1200)
		rho = rho_a * (0.95 - 0.07 * (trial_temp - 400) / 800);
	else
		rho = 0.88 * rho_a;
	
	return rho;
}


double 
NWConcreteEC2::getSpecificHeat(void)
{
    if (trial_temp <= 100.0) {
		cp = 900.0;
		} else if ((100.0 < trial_temp) && (trial_temp <= 200.0)) {
			if (moist == 0.0) {
				cp = trial_temp + 800.0;
				} else if (moist == 0.015) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						cp = 1470.0;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							cp = -5.529411765 * trial_temp + 2105.882353;
						}
				} else if (moist == 0.03) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						cp = 2020.0;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							cp = -12.0 * trial_temp + 3400.0;
						} 
					} else {
						opserr << "NWConcreteEC2::getSpecificHeat() - specific heat values "
							<< "are not available for moisture level " << moist << " .\n";
						exit(-1);
					}
		} else if ((200.0 < trial_temp) && (trial_temp <= 400.0)) {
			cp = 0.5 * trial_temp + 900.0;
			} else {
				cp = 1100.0;
			}

			return cp;
}


double
NWConcreteEC2::getEnthalpy()
{
    // The enthalpy data here is obtained by integrating the spesific heat over
    // the temperature range [20,1200] given in EC 2
    if (trial_temp <= 100.0) {
		double c1 = 2070000.0;
		double c2 = 41400000.0;
		enthalpy = c1 * trial_temp - c2;
		} else if ((100.0 < trial_temp) && (trial_temp <= 200.0)) {
			if (moist == 0.0) {
				double c01 = 1150.0;
				double c02 = 1840000.0;
				double c03 = -29900000.0;
				enthalpy = c01 * trial_temp * trial_temp + c02 * trial_temp + c03;
				} else if (moist == 0.015) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						double c11 = 3381000.0;
						double c12 = -172500000.0;
						enthalpy = c11 * trial_temp + c12;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							double c13 = -6358.823529;
							double c14 = 4843529.412;
							double c15 = -256595441.2;
							enthalpy = c13 * trial_temp * trial_temp + c14 * trial_temp
								+ c15;
						}
				} else if (moist == 0.03) {
					if ((100.0 < trial_temp) && (trial_temp <= 115.0)) {
						double c31 = 4646000.0;
						double c32 = -299000000.0;
						enthalpy = c31 * trial_temp + c32;
						} else if ((115.0 < trial_temp) && (trial_temp <= 200.0)) {
							double c33 = -13800.0;
							double c34 = 7820000.0;
							double c35 = -481505000.0;
							enthalpy = c33 * trial_temp * trial_temp + c34 * trial_temp
								+ c35;
						} else {
							opserr << "NWConcreteEC2::getEnthalpy() - enthalpy values "
								<< "are not available for moisture level " << moist << " .\n";
							exit(-1);
							}
					}
		} else if ((200.0 < trial_temp) && (trial_temp <= 400.0)) {
			double cc0 = 575.0;
			double cc1 = 2070000.0;
			if (moist == 0.0) {
				double c4 = -52900000.0;
				enthalpy = cc0 * trial_temp * trial_temp + cc1 * trial_temp + c4;
				} else if (moist == 0.015) {
					double c5 = 20757500.0;
					enthalpy = cc0 * trial_temp * trial_temp + cc1 * trial_temp + c5;
				} else if (moist == 0.03) {
					double c6 = 93495000.0;
					enthalpy = cc0 * trial_temp * trial_temp + cc1 * trial_temp + c6;
					} else {
						opserr << "NWConcreteEC2::getEnthalpy() - enthalpy values "
							<< "are not available for moisture level " << moist << " .\n";
						exit(-1);
					}
			} else if (400.0 < trial_temp) {
				double cc2 = 2530000.0;
				if (moist == 0.0) {
					double c7 = -144900000.0;
					enthalpy = cc2 * trial_temp + c7;
					} else if (moist == 0.015) {
						double c8 = -71242500.0;
						enthalpy = cc2 * trial_temp + c8;
					} else if (moist == 0.03) {
						double c9 = 1495000.0;
						enthalpy = cc2 * trial_temp + c9;
						} else {
							opserr << "NWConcreteEC2::getEnthalpy() - enthalpy values "
								<< "are not available for moisture level " << moist << " .\n";
							exit(-1);
						}
			} 

			return enthalpy;	
}


double
NWConcreteEC2::getEnthalpy(double temp)
{
    double enthp;
	double nod_temp = temp - 273.15;

	// The temperature range is expanded other than the original one [20,1200] in Eurocode.
	// The reason is, for an analysis with initial temperature at 20, the solution could be lower than
	// 20 after initial iterations. Eventhough, the slope of H-T within the expanded temperature range is 
	// kept constant, the same as the heat capacity at T = 20;
	//if ((0.0 <= nod_temp) && (nod_temp <= 100.0)) {
	if (nod_temp <= 100.0) {
		double c1 = 2070000.0;
		double c2 = 41400000.0;
		enthp = c1 * nod_temp - c2;
		} else if ((100.0 < nod_temp) && (nod_temp <= 200.0)) {
			if (moist == 0.0) {
				double c01 = 1150.0;
				double c02 = 1840000.0;
				double c03 = -29900000.0;
				enthp = c01 * nod_temp * nod_temp + c02 * nod_temp + c03;
				} else if (moist == 0.015) {
					if ((100.0 < nod_temp) && (nod_temp <= 115.0)) {
						double c11 = 3381000.0;
						double c12 = -172500000.0;
						enthp = c11 * nod_temp + c12;
						} else if ((115.0 < nod_temp) && (nod_temp <= 200.0)) {
							double c13 = -6358.823529;
							double c14 = 4843529.412;
							double c15 = -256595441.2;
							enthp = c13 * nod_temp * nod_temp + c14 * nod_temp
								+ c15;
						}
				} else if (moist == 0.03) {
					if ((100.0 < nod_temp) && (nod_temp <= 115.0)) {
						double c31 = 4646000.0;
						double c32 = -299000000.0;
						enthp = c31 * nod_temp + c32;
						} else if ((115.0 < nod_temp) && (nod_temp <= 200.0)) {
							double c33 = -13800.0;
							double c34 = 7820000.0;
							double c35 = -481505000.0;
							enthp = c33 * nod_temp * nod_temp + c34 * nod_temp
								+ c35;
						} 
					} else {
						opserr << "NWConcreteEC2::getEnthalpy(double ) - enthalpy values "
							<< "are not available for moisture level " << moist << " .\n";
					}
		} else if ((200.0 < nod_temp) && (nod_temp <= 400.0)) {
			double cc0 = 575.0;
			double cc1 = 2070000.0;
			if (moist == 0.0) {
				double c4 = -52900000.0;
				enthp = cc0 * nod_temp * nod_temp + cc1 * nod_temp + c4;
				} else if (moist == 0.015) {
					double c5 = 20757500.0;
					enthp = cc0 * nod_temp * nod_temp + cc1 * nod_temp + c5;
				} else if (moist == 0.03) {
					double c6 = 93495000.0;
					enthp = cc0 * nod_temp * nod_temp + cc1 * nod_temp + c6;
					} else {
						opserr << "NWConcreteEC2::getEnthalpy(double ) - enthalpy values "
							<< "are not available for moisture level " << moist << " .\n";
					}
				//} else if ((400.0 < nod_temp) && (nod_temp <= 1200.0)) {
			} else if (400.0 < nod_temp) {
				double cc2 = 2530000.0;
				if (moist == 0.0) {
					double c7 = -144900000.0;
					enthp = cc2 * nod_temp + c7;
					} else if (moist == 0.015) {
						double c8 = -71242500.0;
						enthp = cc2 * nod_temp + c8;
					} else if (moist == 0.03) {
						double c9 = 1495000.0;
						enthp = cc2 * nod_temp + c9;
						} else {
							opserr << "NWConcreteEC2::getEnthalpy(double ) - enthalpy values "
								<< "are not available for moisture level " << moist << " .\n";
						}
			}

			return enthp;	
}


HeatTransferMaterial*
NWConcreteEC2::getCopy(void)
{
    NWConcreteEC2* theCopy = new NWConcreteEC2(this->getTag(), moist);
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
NWConcreteEC2::update()
{
    return; 
}


int 
NWConcreteEC2::commitState(void)
{
    return 0;
}


int 
NWConcreteEC2::revertToLastCommit(void)
{
    return 0;
}


int 
NWConcreteEC2::revertToStart(void)
{
    return 0;
}


