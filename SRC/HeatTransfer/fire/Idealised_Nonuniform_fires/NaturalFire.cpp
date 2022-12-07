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
// Written by Liming Jiang (liming.jiang@poly.edu.hk)
//

#include <NaturalFire.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <PathTimeSeriesThermal.h>


NaturalFire::NaturalFire(int tag, double D,
	double Q, double H, int lineTag, double smokeTemp, PathTimeSeriesThermal* FireParPath)
	:FireModel(tag, 7), FireParPath(FireParPath), fireLocs(3), d(D), hc(0), absorp(0),
	q(Q), h(H), smokeT(smokeTemp), addq(1e5), centerLine(lineTag), fmode(0)
{
	// check the direction of central line of a Hasemi fire
	// 1 indicates it is parrallel to x1 axis, 2 indicates
	// parallelt to x2 axis, 3 indicates parallel to x3 axis.
	if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "NaturalFire::NaturalFire - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
	}
	if ((d > 10.0) || (q > 5e7)) {
		opserr << "NaturalFire::NaturalFire - error in specifying the fire, diameter "
			<< " shoudn't be greater than 10m, fire size shoudn't be greater than 50MW.\n";
	}
}

NaturalFire::NaturalFire(int tag, int lineTag, PathTimeSeriesThermal* FireParPath)
	:FireModel(tag, 7), FireParPath(FireParPath), fireLocs(3), d(0.0), hc(0), absorp(0),
	q(0.0), h(0.0), smokeT(0.0), addq(0.0), centerLine(lineTag), fmode(0)
{
	// check the direction of central line of a Hasemi fire
	// 1 indicates it is parrallel to x1 axis, 2 indicates
	// parallelt to x2 axis, 3 indicates parallel to x3 axis.
	if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "NaturalFire::NaturalFire - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
	}

}

NaturalFire::~NaturalFire()
{

}



int
NaturalFire::setFirePars(double time, const Vector& firePars) {

	Vector FirePars = 0;
	if (FireParPath != 0) {
		FirePars = FireParPath->getFactors(time);
	}
	else if (FireParPath == 0 && firePars != 0)
		FirePars = firePars;


	if (FirePars.Size() == 3) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
	}
	else if (FirePars.Size() == 4) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
	}
	else if (FirePars.Size() == 5) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
		d = FirePars(4);
	}
	else if (FirePars.Size() == 6) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
		d = FirePars(4);
		smokeT = FirePars(5) + 273.15;
	}
	else if (FirePars.Size() == 7) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
		d = FirePars(4);
		smokeT = FirePars(5) + 273.15;
		addq = FirePars(6);

	}
	else if (FirePars.Size() == 8) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
		d = FirePars(4);
		h = FirePars(5);
		smokeT = FirePars(6) + 273.15;
		addq = FirePars(7);

	}
	else {
		opserr << "WARNING! NaturalFire::getFlux failed to get the location of fire origin" << endln;
		return -1;
	}
#ifdef _DEBUG
	//opserr << FirePars << endln;
#endif // DEBUG


	return 0;
}

double
NaturalFire::getFirePars(int ParTag) {
	if (ParTag == 1)
		return fireLocs(0);
	else if (ParTag == 2)
		return fireLocs(1);
	else if (ParTag == 3)
		return fireLocs(2);
	else if (ParTag == 4)
		return q;
	else if (ParTag == 5)
		return d;
	else if (ParTag == 6)
		return smokeT;
	else if (ParTag == 7)
		return addq;
	else {
		opserr << "WARNING! invalid tag for NaturalFire::getFirePars " << ParTag << endln;
		return -1;
	}


}

double
NaturalFire::getFireOut(double time, const Vector& coords)
{
	if (FireParPath != 0) {
		if (this->setFirePars(time) < 0)
			exit(-1);
	}

	double gas_t = 0;
	double q_dot = 0;
	double Addqs = addq;


	// now calculate r	

	int size = coords.Size();

	double deltaX1, deltaX2, sum, r;

	if (centerLine == 1) {
		deltaX1 = fireLocs(1) - coords(1);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = fireLocs(2) - coords(2);
		}
		else {
			// then heat transfer coordinate is 2D
			deltaX2 = fireLocs(2);
		}
	}
	else if (centerLine == 2) {
		deltaX1 = fireLocs(0) - coords(0);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = fireLocs(2) - coords(2);
		}
		else {
			// then heat transfer coordinate is 2D
			deltaX2 = fireLocs(2);
		}
	}
	else if (centerLine == 3) {
		deltaX1 = fireLocs(0) - coords(0);
		deltaX2 = fireLocs(1) - coords(1);
	}

	sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
	r = sqrt(sum);
	//--------------completion of calculating r-------------------------------

	if (abs(h) < 1e-6)
		opserr << "Travelling fire: h got a zero " << endln;

	if (q < 0) {
		if (-q < 600) {
			fmode = 4;  //deca mode if gas temperature is lower than 600oC
			hc = 25;
		}
		else if (-q > 600) {
			fmode = 3;  //localised flashover:Mode 3
			hc = 15;
		}
		
		gas_t = -q + 273.15;    //absolute temperature
		
		if (r < d / 2.0)
			q_dot = absorp * 5.67e-8 * (pow(gas_t, 4) - pow(293.15, 4)) + hc * (gas_t - 293.15);
		else if (r < d / 2.0 + 0.5) {
			double q_dot_n = absorp * 5.67e-8 * (pow(gas_t, 4) - pow(293.15, 4)) + hc * (gas_t - 293.15);
			double q_dot_s = q_dot + absorp * 5.67e-8 * (pow(smokeT, 4) - pow(293.15, 4)) + hc * (smokeT - 293.15);
			q_dot = q_dot_n + (r - d / 2.0) / 0.5 * (q_dot_s - q_dot_n);
		}

		else
			q_dot = absorp * 5.67e-8 * (pow(smokeT, 4) - pow(293.15, 4)) + 25 * (smokeT - 293.15);

	}
	else {
		//localised fire
	    // first calculate flame length
		double Lf = 0.0148 * pow(q, 0.4) - 1.02 * d;

		double constant = 1.11 * 1e6 * pow(d, 2.5);
		double Qd_ast = q / constant;
		double z_acute;

		if (Qd_ast < 1.0) {
			double a = 2.0 / 3.0;
			double term = pow(Qd_ast, 0.4) - pow(Qd_ast, a);
			z_acute = 2.4 * d * term;
		}
		else {
			double term = 1.0 - pow(Qd_ast, 0.4);
			z_acute = 2.4 * d * term;
		}
		double term = 1.11 * 1e6 * pow(h, 2.5);
		double Qh_ast = q / term;

		// now calculate H plus Lh
		double Lt = 2.9 * h * pow(Qh_ast, 0.33);


		double q_smoke = absorp * 5.67e-8 * (pow(smokeT, 4) - pow(293.15, 4)) + hc * (smokeT - 293.15);

		if (Lf < h) {
			//------------------not impinge ceiling---Alpert-------------------------------------------
			fmode = 1;
			hc = 15;
			if (r / h > 0.18) {
				gas_t = 5.38 * pow(q / 1000.0 / r, 2.0 / 3.0) / h;
			}
			else {
				gas_t = 16.9 * pow(q / 1000.0, 2.0 / 3.0) / pow(h, 5.0 / 3.0);

			}

			//opserr << " GAS: " << gas_t << " Q: "<< q<<" h "<< pow(q / 1000.0 / r, 2.0 / 3.0);
			gas_t = gas_t + 293.15;


			//if (Addqs > q_smoke)
				//Addqs = q_smoke;

			//if (r < d)
			q_dot = absorp * 5.67e-8 * (pow(gas_t, 4) - pow(293.15, 4)) + hc * (gas_t - 293.15) + q_smoke;
			//else
				//q_dot = 0.85 * 5.67e-8 * (pow(gas_t, 4) - pow(293.15, 4)) + 35 * (gas_t - 293.15);


			if (q_dot < q_smoke) {
#ifdef _DEBUG
				opserr << "Natural fire: q_dot " << q_dot << "q_smoke: " << q_smoke << endln;
#endif
				q_dot = q_smoke;

			}

		}
		else {
			//---------------------impinge ceiling--------------------------------------
			fmode = 2;
			// now calculate y
			hc = 15;
			double y = (r + h + z_acute) / (Lt + z_acute);

			// now determine the flux

			if (y <= 0.3) {
				q_dot = 100000;
			}
			else if ((y > 0.3) && (y < 1.0)) {
				q_dot = 136300 - 121000 * y;
			}
			else if (y >= 1.0) {

				if (r / h > 0.18) {
					gas_t = 5.38 * pow(q / 1000.0 / r, 2.0 / 3.0) / h;
				}
				else {
					gas_t = 16.9 * pow(q / 1000.0, 2.0 / 3.0) / pow(h, 5.0 / 3.0);

				}

				//opserr << " GAS: " << gas_t << " Q: "<< q<<" h "<< pow(q / 1000.0 / r, 2.0 / 3.0);
				gas_t = gas_t + 293.15;

				double q_dot_ec = 15000 * pow(y, -3.7);
				double q_dot_gt = absorp * 5.67e-8 * (pow(gas_t, 4) - pow(293.15, 4)) + hc * (gas_t - 293.15);

				if (q_dot_ec < q_dot_gt) {
					q_dot = q_dot_gt;
				}
				else {
					q_dot = q_dot_ec;
				}
			}

			//if (Addqs > q_smoke)
				//Addqs = q_smoke;
			q_dot = q_dot + q_smoke; //modify the maximum q only by smoke induced q

			if (q_dot > 120000)
				q_dot = 120000;

			//adibadic temperature principle: eps*qr -eps*sigma*T^4 +h (smokeT-T)=0
			// Gauge heat flux = eps*qr-eps*sigma*Tg^4+h(smokeT-Tg)
			if (q_dot < q_smoke) {
#ifdef _DEBUG
				opserr << "Travelling fire: q_dot " << q_dot << "q_smoke: " << q_smoke << endln;
#endif
				q_dot = q_smoke;

			}
		}


#ifdef _DEBUG
		//int tag = node->getTag();
		//opserr<<" Travelling fire: "<<fireLocs(1)<<","<<x3<<", q"<<q<<" r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif


	}

	return q_dot;

}


double
NaturalFire::getFlux(HeatTransferNode* node, double time)
{
	const Vector& coords = node->getCrds();
	double qdot = this->getFireOut(time, coords);
	return qdot;
}


void
NaturalFire::applyFluxBC(HeatFluxBC* theFlux, double time)
{
	int flux_type = theFlux->getTypeTag();

	if (flux_type == 1) {
		Convection* convec = (Convection*)theFlux;
		if (fmode == 1 || fmode == 2)
			hc = 15;
		else if (fmode == 3)
			hc = 25;
		else if (fmode == 4)
			hc = 15;
		else
			opserr << "NaturalFire:: incorrect fire mode" << endln;
		convec->setParameter(hc); //set hc for convection;
		convec->applyFluxBC(time);
	}
	else if (flux_type == 2) {
		Radiation* rad = (Radiation*)theFlux;
		//double bzm = 5.67e-8;
		//double temp = rad->;
		//double temp = this->getGasTemperature(time);
		//double qir = bzm * pow(ambTemp, 4.0);
		//rad->setIrradiation(qir);
		if (abs(absorp) < 1e-5)
			absorp = rad->getAbsorptivity();
		rad->applyFluxBC(time);
	}
	else if (flux_type == 3)
	{
		PrescribedSurfFlux* pflux = (PrescribedSurfFlux*)theFlux;

		//int flux_type = pflux->getTypeTag();
		int eleTag = pflux->getElementTag();
		int fTag = pflux->getFaceTag();
		HeatTransferDomain* theDomain = pflux->getDomain();
		if (theDomain == 0) {
			opserr << "NaturalFire::applyFluxBC() - HeatFluxBC has not been associated with a domain";
			exit(-1);
		}

		HeatTransferElement* theEle = theDomain->getElement(eleTag);
		if (theEle == 0) {
			opserr << "NaturalFire::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
			exit(-1);
		}

		const ID& faceNodes = theEle->getNodesOnFace(fTag);
		int size = faceNodes.Size();
		Vector nodalFlux(size);
		if (abs(absorp) < 1e-5)
			absorp = 0.85;     //if not assigned

		if (abs(hc) < 1e-5)
			hc = 25;      //if not assigned.

		for (int i = 0; i < size; i++) {
			int nodTag = faceNodes(i);
			HeatTransferNode* theNode = theDomain->getNode(nodTag);
			if (theNode == 0) {
				opserr << "NaturalFire::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
				exit(-1);
			}
			nodalFlux(i) = this->getFlux(theNode, time);
			//start the prescribed heat flux
			//opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
		}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
	}
	else {
		opserr << "NaturalFire::applyFluxBC() - incorrect flux type "
			<< flux_type << " provided\n";
		exit(-1);
	}
}

