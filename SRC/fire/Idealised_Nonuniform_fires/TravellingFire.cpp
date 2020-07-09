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

#include <TravellingFire.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <PathTimeSeriesThermal.h>


TravellingFire::TravellingFire(int tag, double D,
	double Q, double H, int lineTag, double smokeTemp, PathTimeSeriesThermal* fireLocPath)
	:FireModel(tag, 7), FireLocPath(fireLocPath), fireLocs(3), d(D), 
	q(Q), h(H), smokeT(smokeTemp),maxq(1e5),centerLine(lineTag)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "TravellingFire::TravellingFire - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}
	if ((d > 10.0) || (q > 5e7)) {
		opserr << "TravellingFire::TravellingFire - error in specifying the fire, diameter "
			<< " shoudn't be greater than 10m, fire size shoudn't be greater than 50MW.\n";
		}
}


TravellingFire::~TravellingFire()
{

}



int
TravellingFire::setFirePars(double time, const Vector& firePars) {
	
	Vector FirePars = 0;
	if (FireLocPath != 0) {
		FirePars = FireLocPath->getFactors(time);
	}
	else if (FireLocPath==0&&firePars != 0)
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
		smokeT = FirePars(5);
	}
	else if (FirePars.Size() == 7) {
		fireLocs(0) = FirePars(0);
		fireLocs(1) = FirePars(1);
		fireLocs(2) = FirePars(2);
		q = FirePars(3);
		d = FirePars(4);
		smokeT = FirePars(5);
		maxq = FirePars(6);

	}
	else {
			opserr << "WARNING! TravellingFire::getFlux failed to get the location of fire origin" << endln;
			return -1;
		}
#ifdef _DEBUG
	opserr << FirePars << endln;
#endif // DEBUG

	
		return 0;
}

double
TravellingFire::getFirePars(int ParTag) {
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
		return maxq;
	else {
		opserr << "WARNING! invalid tag for TravellingFire::getFirePars " << ParTag << endln;
		return -1;
	}
		

}

double 
TravellingFire::getFireOut( double time, const Vector& coords)
{
	if (FireLocPath != 0) {
		if (this->setFirePars(time) < 0)
			exit(-1);
	}


	double q_dot = 0;

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

	// now calculate y
	double y = (r + h + z_acute) / (Lt + z_acute);

	// now determine the flux

	if (y <= 0.3) {
		q_dot = 100000;
	}
	else if ((y > 0.3) && (y < 1.0)) {
		q_dot = 136300 - 121000 * y;
	}
	else if (y >= 1.0) {
		q_dot = 15000 * pow(y, -3.7);
	}

	q_dot = q_dot * maxq / 1e5; //modify the maximum q

	if (Lf < h) {
		opserr << "Lf: "<<Lf<<" h: "<<h<<endln;
		q_dot = 0.0;
	}


	double q_smoke = 0.8 * 5.67e-8 * (pow(smokeT, 4) - pow(293.15, 4)) + 35 * (smokeT - 293.15);
	if (q_dot < q_smoke) {
#ifdef _DEBUG
		opserr << "Travelling fire: q_dot " << q_dot << "q_smoke: " << q_smoke << endln;
#endif
		q_dot = q_smoke;

	}



#ifdef _DEBUG
	//int tag = node->getTag();
	//opserr<<" Travelling fire: "<<fireLocs(1)<<","<<x3<<", q"<<q<<" r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif
	return q_dot;

}


double 
TravellingFire::getFlux(HeatTransferNode* node, double time)
{
	const Vector& coords = node->getCrds();
	double qdot = this->getFireOut(time,coords);
	return qdot;
}


void
TravellingFire::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    int flux_type = theFlux->getTypeTag();
    if (flux_type == 3) 
		{
		PrescribedSurfFlux* pflux = (PrescribedSurfFlux*) theFlux;

		//int flux_type = pflux->getTypeTag();
		int eleTag = pflux->getElementTag();
		int fTag = pflux->getFaceTag();
		HeatTransferDomain* theDomain = pflux->getDomain();
		if (theDomain == 0) {
			opserr << "TravellingFire::applyFluxBC() - HeatFluxBC has not been associated with a domain";
			exit(-1);
			}

		HeatTransferElement* theEle = theDomain->getElement(eleTag);
		if (theEle == 0) {
			opserr << "TravellingFire::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
			exit(-1);
			}

		const ID& faceNodes = theEle->getNodesOnFace(fTag);
		int size = faceNodes.Size();
		Vector nodalFlux(size);

		for (int i = 0; i < size; i++) {
			int nodTag = faceNodes(i);
			HeatTransferNode* theNode = theDomain->getNode(nodTag);
			if (theNode == 0) {
				opserr << "TravellingFire::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
				exit(-1);
				}
			nodalFlux(i) = this->getFlux(theNode,time); 
			//opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
			}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
		} else {
			opserr << "TravellingFire::applyFluxBC() - incorrect flux type "
				<< flux_type << " provided\n";
			exit(-1);
		}
}

