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

#include <LocalizedFireEC1.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>


LocalizedFireEC1::LocalizedFireEC1(int tag, double crd1, double crd2, double crd3, double D,
								   double Q, double H, int lineTag,bool forColumn,double startTime)
:FireModel(tag,3),x1(crd1), x2(crd2), x3(crd3), d(D), q(Q), h(H), centerLine(lineTag),ForColumn(forColumn)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "LocalizedFireEC1::LocalizedFireEC1 - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}
	if ((d > 10.0) || (q > 5e7)) {
		opserr << "LocalizedFireEC1::LocalizedFireEC1 - error in specifying the fire, diameter "
			<< " shoudn't be greater than 10m, fire size shoudn't be greater than 50MW.\n";
		}
}


LocalizedFireEC1::~LocalizedFireEC1()
{

}


double 
LocalizedFireEC1::getFlux(HeatTransferNode* node, double time)
{
    // first calculate flame length
    double Lf = 0.0148 * pow(q,0.4) - 1.02 * d;
	if (Lf < h) {
		opserr << "LocalizedFireEC1::getFlux() - flame is not impinging ceiling, method has not implemented.\n";
		return -1;
		}
	double constant = 1.11 * 1e6 * pow(d,2.5);
	double Qd_ast = q / constant;
	double z_acute;

	if (Qd_ast < 1.0) {
		double a = 0.66666666666666666666666666666667;
		double term = pow(Qd_ast,0.4) - pow(Qd_ast,a);
		z_acute = 2.4 * d * term;
		} else {
			double term = 1.0 - pow(Qd_ast, 0.4);
			z_acute = 2.4 * d * term;
		}
	double term = 1.11 * 1e6 * pow(h,2.5);
	double Qh_ast = q / term;

	// now calculate H plus Lh
	double Lt = 2.9 * h * pow(Qh_ast,0.33);

	// now calculate r	
	const Vector& coords = node->getCrds();
	int size = coords.Size();

	double deltaX1, deltaX2, sum, r;

	if (centerLine == 1) {
		deltaX1 = x2 - coords(1);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = x3 - coords(2);
			} else {
				// then heat transfer coordinate is 2D
				deltaX2 = x3;
			}
	} else if (centerLine == 2) {
		deltaX1 = x1 - coords(0);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = x3 - coords(2);
		} else {
			// then heat transfer coordinate is 2D
			deltaX2 = x3;
		}
	} else if (centerLine == 3) {
			deltaX1 = x1 - coords(0);
			deltaX2 = x2 - coords(1);
			}

    sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
    r = sqrt(sum);

	// now calculate y
	double y = (r + h + z_acute) / (Lt + z_acute);

	// now determine the flux
	double q_dot;
	if (y <= 0.3) {
		q_dot = 100000;
	} else if ((y > 0.3) && (y < 1.0)) {
		q_dot = 136300 - 121000 * y;
	} else if (y >= 1.0) {
		q_dot = 15000 * pow(y,-3.7);
	}
	int tag = node->getTag();
#ifdef _DEBUG
	//opserr<<" NodeTag: "<<node->getTag()<< " r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif
	return q_dot;
}


void
LocalizedFireEC1::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    int flux_type = theFlux->getTypeTag();
   
	if (flux_type == 1) {
		Convection* convec = (Convection*)theFlux;
		//convec->setSurroundingTemp(this->getGasTemperature(time));
		int eleTag = theFlux->getElementTag();
		int fTag = theFlux->getFaceTag();
		HeatTransferDomain* theDomain = theFlux->getDomain();
		if (theDomain == 0) {
			opserr << "LocalizedFireSFPE::applyFluxBC() - HeatFluxBC has not been associated with a domain";
			exit(-1);
		}

		HeatTransferElement* theEle = theDomain->getElement(eleTag);
		if (theEle == 0) {
			opserr << "LocalizedFireSFPE::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
			exit(-1);
		}

		const ID& faceNodes = theEle->getNodesOnFace(fTag);
		HeatTransferNode* theNode = theDomain->getNode(faceNodes(0));
		convec->applyFluxBC(time);
	}
	else if (flux_type == 2) {
		Radiation* rad = (Radiation*)theFlux;
		rad->applyFluxBC(time);
	}
	else if (flux_type == 3) 
	{
		PrescribedSurfFlux* pflux = (PrescribedSurfFlux*) theFlux;

		//int flux_type = pflux->getTypeTag();
		int eleTag = pflux->getElementTag();
		int fTag = pflux->getFaceTag();
		HeatTransferDomain* theDomain = pflux->getDomain();
		if (theDomain == 0) {
			opserr << "LocalizedFireEC1::applyFluxBC() - HeatFluxBC has not been associated with a domain";
			exit(-1);
		}

		HeatTransferElement* theEle = theDomain->getElement(eleTag);
		if (theEle == 0) {
			opserr << "LocalizedFireEC1::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
			exit(-1);
		}

		const ID& faceNodes = theEle->getNodesOnFace(fTag);
		int size = faceNodes.Size();
		Vector nodalFlux(size);

		for (int i = 0; i < size; i++) {
			int nodTag = faceNodes(i);
			HeatTransferNode* theNode = theDomain->getNode(nodTag);
			if (theNode == 0) {
				opserr << "LocalizedFireEC1::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
				exit(-1);
			}
			nodalFlux(i) = this->getFlux(theNode,time); 
			//opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
		}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
	}
	else {
			opserr << "LocalizedFireEC1::applyFluxBC() - incorrect flux type "
				<< flux_type << " provided\n";
			exit(-1);
	}
}

void 
LocalizedFireEC1::Print(OPS_Stream& s, int i)
{
	opserr << "Fire model: "<<this->getTag()<< " , Type:EC1 localised" << endln
		<< "Origin: " << x1 << ", " << x2 << ", " << x3 << endln
		<< "D: " << d << " ,Q: " << q << " ,H: " << h <<" centerline: "<< centerLine << endln;

}