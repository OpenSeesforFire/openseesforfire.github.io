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
**   Liming Jiang(liming.jiang@ed.ac.uk)                              **
** ****************************************************************** */

//
// Added by Liming Jiang (y.jiang@ed.ac.uk)
//

#include <Idealised_Local_Fire.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>


Idealised_Local_Fire::Idealised_Local_Fire(int tag,double crd1, double crd2, double crd3,
								   double Q, double D1, double D2, int lineTag)
:FireModel(tag,6),x1(crd1), x2(crd2), x3(crd3), d1(D1),d2(D2), k1(0),k2(0), q(Q), centerLine(lineTag), CurveType(1), FireLocPath(0)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "Idealised_Local_Fire::Idealised_Local_Fire - invalid line tag\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}

	if (D2 - D1 > 1e-3)
		CurveType = 1;
	else
		CurveType = 0;
	
}

Idealised_Local_Fire::Idealised_Local_Fire(int tag,double crd1, double crd2, double crd3, double Q, double D1, double D2, double K1, double K2, int lineTag)
:FireModel(tag, 6),x1(crd1), x2(crd2), x3(crd3), d1(D1),d2(D2), q(Q), k1(K1), k2(K2), centerLine(lineTag), FireLocPath(0)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "Idealised_Local_Fire::Idealised_Local_Fire - invalid line tag\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}
	if(D2-D1>1e-3)
		CurveType =2;
	else
		CurveType =1;
	
}

Idealised_Local_Fire::Idealised_Local_Fire(int tag,double crd1, double crd2, double crd3, double Q, double D1, double D2, double factor, int lineTag)
:FireModel(tag,6),x1(crd1), x2(crd2), x3(crd3), d1(D1),d2(D2), Factor(factor), q(Q), centerLine(lineTag), CurveType(3), FireLocPath(0)
{
  // check the direction of central line of a Hasemi fire
  // 1 indicates it is parrallel to x1 axis, 2 indicates
  // parallelt to x2 axis, 3 indicates parallel to x3 axis.
  if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
    opserr << "Idealised_Local_Fire::Idealised_Local_Fire - invalid line tag provided for Hasemi fire.\n"
    << " Only 1, or 2, or 3 is correct.\n";
		}
  
}

//adding moving fire source
Idealised_Local_Fire::Idealised_Local_Fire(int tag, PathTimeSeriesThermal* fireLocPath, double Q, double D1, double D2, double factor, int lineTag)
	:FireModel(tag, 6), FireLocPath(fireLocPath), d1(D1), d2(D2), Factor(factor), q(Q), centerLine(lineTag), CurveType(3)
{
	// check the direction of central line of a Hasemi fire
	// 1 indicates it is parrallel to x1 axis, 2 indicates
	// parallelt to x2 axis, 3 indicates parallel to x3 axis.
	if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "Idealised_Local_Fire::Idealised_Local_Fire - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
	}
	x1 = 0; x2 = 0; x3 = 0;

}


Idealised_Local_Fire::~Idealised_Local_Fire()
{

}


double 
Idealised_Local_Fire::getFlux(HeatTransferNode* node, double time)
{
	// now determine the flux
	double q_dot;
  
	if (CurveType == 0) {
		q_dot = q;
		//Uniform distribution of heat flux

		return q_dot;

	}

	// now calculate r	
	const Vector& coords = node->getCrds();
	int size = coords.Size();

	

	double deltaX1, deltaX2, sum, r;
	//------------moving ..................
	Vector fireLocs = 0;
	if (FireLocPath != 0)
	{
		fireLocs = FireLocPath->getFactors(time);

		if (fireLocs.Size() == 1) {
			x1 = fireLocs(0);
			x2 = 0;
		}
		else if (fireLocs.Size() == 2) {
			x1 = fireLocs(0);
			x2 = fireLocs(1);
		}
		else if (fireLocs.Size() == 3) {
			x1 = fireLocs(0);
			x2 = fireLocs(1);
			x3 = fireLocs(2);
		}
		else {
			opserr << "WARNING! IdealisedLoalFire::getFlux failed to get the location of fire origin" << endln;
		}
	}
	//-----------------
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

	
  if (r <= d1) {
		q_dot = q;
	}
  
  else if ((r > d1) && (r <= d2)) {
		
	  
	 if(CurveType==1){
		q_dot = q - q * (r-d1)/(d2-d1);
	 	
		}
     else if(CurveType==2){
		double L=d2-d1;
		double A = 2*(k2-2*k1+1)/L/L;
		double B = (-k2+4*k1-3)/L;
		double X=r-d1;
		q_dot = (A*X*X+B*X+1)*q;
    }
    else if(CurveType==3){
      double X= (r-d1)/(d2-d1);
      q_dot = q* (exp(-Factor*X));
    }
		
	} else if (r > d2) {
		q_dot = 0;
	}
	int tag = node->getTag();
#ifdef _DEBUG
	opserr<<" NodeTag: "<<node->getTag()<< " r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif
	return q_dot;
}


void
Idealised_Local_Fire::applyFluxBC(HeatFluxBC* theFlux, double time)
{
    int flux_type = theFlux->getTypeTag();
	if (flux_type == 1) {
		Convection* convec = (Convection*)theFlux;
		//convec->setSurroundingTemp(this->getGasTemperature(time));
		int eleTag = theFlux->getElementTag();
		int fTag = theFlux->getFaceTag();
		HeatTransferDomain* theDomain = theFlux->getDomain();
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
			nodalFlux(i) = this->getFlux(theNode,time); 
			//opserr << "Flux at node " << nodTag << " is " << nodalFlux(i) << endln;
			}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
	} else {
			opserr << "Idealised_Local_Fire::applyFluxBC() - incorrect flux type "
				<< flux_type << " provided\n";
			exit(-1);
	}
}

