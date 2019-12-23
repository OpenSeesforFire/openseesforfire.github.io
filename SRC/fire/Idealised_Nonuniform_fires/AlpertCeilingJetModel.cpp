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

#include <AlpertCeilingJetModel.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>


AlpertCeilingJetModel::AlpertCeilingJetModel(int tag,
                        double crd1, double crd2, double crd3,
											 double Q, double H, double Ta, 
											 double Tflame, int lineTag)
:FireModel(tag,5), x1(crd1), x2(crd2), x3(crd3), q(Q/1000.0), h(H), T0(Ta),
 Tf(Tflame), centerLine(lineTag)
{
    // check the direction of central line of the Alpert Ceiling Jet Model
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "AlpertCeilingJetModel::AlpertCeilingJetModel - "
			<< "invalid line tag provided for Alpert fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}

}


AlpertCeilingJetModel::~AlpertCeilingJetModel()
{

}


double 
AlpertCeilingJetModel::getGasTemperature(HeatFluxBC* the_flux, double time)
{
    //first identify nodes on each face
    int eleTag = the_flux->getElementTag();
	int fTag = the_flux->getFaceTag();
	HeatTransferDomain* theDomain = the_flux->getDomain();
	if (theDomain == 0) {
		opserr << "AlpertCeilingJetModel::applyFluxBC() - HeatFluxBC has not been associated with a domain";
		exit(-1);
		}
	
    HeatTransferElement* theEle = theDomain->getElement(eleTag);
	if (theEle == 0) {
		opserr << "AlpertCeilingJetModel::applyFluxBC() - no element with tag " << eleTag << " exists in the domain";
		exit(-1);
		}
    const ID& faceNodes = theEle->getNodesOnFace(fTag);
	int size = faceNodes.Size();
    
	double crd1 = 0.0;
	double crd2 = 0.0;
	double crd3 = 0.0;
	static int size2;

	for (int i = 0; i < size; i++) {
		int nodTag = faceNodes(i);
		HeatTransferNode* theNode = theDomain->getNode(nodTag);
		if (theNode == 0) {
			opserr << "AlpertCeilingJetModel::getGasTemperature() - no node with tag " << nodTag 
				<< " exists in the domain";
			exit(-1);
			}

		const Vector& coords = theNode->getCrds();
		size2 = coords.Size();
		crd1 = crd1 + coords(0);
		crd2 = crd2 + coords(1);
		if(size2 == 3)
			crd3 = crd3 + coords(2);
		}

	crd1 = crd1 / size;
	crd2 = crd2 / size;
	if(size2 == 3)
		crd3 = crd3 /size;

    // first calculate r	
	double deltaX1, deltaX2, sum, r;

	if (centerLine == 1) {
		//deltaX1 = x2 - crd1;
		deltaX1 = x2 - crd2;
		if (size2 == 3) {
			// then heat transfer coordinate is 3D
			//deltaX2 = x3 - crd2;
			deltaX2 = x3 - crd3;
			} else {
				// then heat transfer coordinate is 2D
				deltaX2 = x3;
			}
		} else if (centerLine == 2) {
			deltaX1 = x1 - crd1;
			if (size2 == 3) {
				// then heat transfer coordinate is 3D
				//deltaX2 = x3 - crd2;
				deltaX2 = x3 - crd3;
				} else {
					// then heat transfer coordinate is 2D
					deltaX2 = x3;
				}
		} else if (centerLine == 3) {
			deltaX1 = x1 - crd1;
			deltaX2 = x2 - crd2;
			}

		sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
		r = sqrt(sum);

		// now calculate ratio
		double ratio = r / h;

		double dT, Tmax;

		if(ratio <= 0.18){
			dT = 16.9 * pow(q,2.0/3.0) / pow(h,5.0/3.0);
			}else if(ratio > 0.18){
			dT = 5.38 * pow(q/r, 2.0/3.0)/h;
			}
		Tmax = dT + T0;
		if(Tmax > Tf)
			Tmax = Tf;
		return Tmax;
}


double 
AlpertCeilingJetModel::getGasTemperature(double xx1, double xx2, double xx3)
{
	
	// first calculate r	
	double deltaX1, deltaX2, sum, r;

	if (centerLine == 1) {
		//deltaX1 = x2 - crd1;
		deltaX1 = x2 - xx2;
		deltaX2 = x3 - xx3;
		} else if (centerLine == 2) {
			deltaX1 = x1 - xx1;
			deltaX2 = x3 - xx3;
		} else if (centerLine == 3) {
			deltaX1 = x1 - xx1;
			deltaX2 = x2 - xx2;
			}

		sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
		r = sqrt(sum);

		// now calculate ratio
		double ratio = r / h;

		double dT, Tmax;

		if(ratio <= 0.18){
			dT = 16.9 * pow(q,2.0/3.0) / pow(h,5.0/3.0);
			}else if(ratio > 0.18){
			dT = 5.38 * pow(q/r, 2.0/3.0)/h;
			}
		Tmax = dT + T0;
		if(Tmax > Tf)
			Tmax = Tf;
		return Tmax;
}

void
AlpertCeilingJetModel::applyFluxBC(HeatFluxBC* theFlux, double time)
{

    // need to determine convection type or radiation
    int flux_type = theFlux->getTypeTag();
	if (flux_type == 1) {
		Convection* convec = (Convection*) theFlux;
		convec->setSurroundingTemp(this->getGasTemperature(theFlux,time));
		convec->applyFluxBC(time);
		} else if (flux_type == 2) {
			Radiation* rad = (Radiation*) theFlux;
			static const double bzm = 5.67 * 1e-008;
			//double alpha = rad->getAbsorptivity();
			double temp = this->getGasTemperature(theFlux,time);
			double qir = bzm * pow(temp, 4.0);
			rad->setIrradiation(qir);
			rad->applyFluxBC(time);
		} else {
			opserr << "AlpertCeilingJetModel::applyFluxBC() - incorrect flux type provided.\n";
			}
	}


