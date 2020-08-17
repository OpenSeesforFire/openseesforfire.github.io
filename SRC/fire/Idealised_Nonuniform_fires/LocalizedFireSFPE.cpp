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
// Added by Liming Jiang (liming.jiang@ed.ac.uk)
// Modified based on LocalisedFireEC1.cpp
//

#include <LocalizedFireSFPE.h>
#include <cmath>
#include <HeatFluxBC.h>
#include <PrescribedSurfFlux.h>
#include <ID.h>
#include <Vector.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>


LocalizedFireSFPE::LocalizedFireSFPE(int tag, double crd1, double crd2, double crd3, double D,
								   double Q, double HC, double HB, int lineTag, double qsquare)
:FireModel(tag, 4),x1(crd1), x2(crd2), x3(crd3), d(D), ini_q(Q), hc(HC),hb(HB), centerLine(lineTag), Qsquare(qsquare)
{
    // check the direction of central line of a Hasemi fire
    // 1 indicates it is parrallel to x1 axis, 2 indicates
    // parallelt to x2 axis, 3 indicates parallel to x3 axis.
    if ((lineTag != 1) && (lineTag != 2) && (lineTag != 3)) {
		opserr << "LocalizedFireSFPE::LocalizedFireSFPE - invalid line tag provided for Hasemi fire.\n"
			<< " Only 1, or 2, or 3 is correct.\n";
		}
	if ((d > 10.0) || (ini_q > 5e7)) {
		opserr << "LocalizedFireSFPE::LocalizedFireSFPE - error in specifying the fire, diameter "
			<< " shoudn't be greater than 10m, fire size shoudn't be greater than 50MW.\n";
		}
	if(Qsquare<0){
		opserr << "LocalizedFireSFPE::LocalizedFireSFPE - error in specifying the square factor "
			<< " shoudn't be negative.\n";

	}
	//FireType = 0: local fire impingeing on the ceiling
	//FireType = 1;
}


LocalizedFireSFPE::~LocalizedFireSFPE()
{

}


double 
LocalizedFireSFPE::getFlux(HeatTransferNode* node, double time, int FireType)
{
    // first calculate flame length
	double q;
	//if Qsquare is activated
	if(Qsquare>1)
		q = (time/Qsquare)*(time/Qsquare)*1e6;
	else 
		q = ini_q;
	//applying Qsquare
	if(q>ini_q)
		q = ini_q;

	double constant = 1.11 * 1e6 * pow(d,2.5);
	double Qd_ast = q / constant;
	
	double z_acute, Lf;

	if (Qd_ast < 1.0) {
		double term = pow(Qd_ast,0.4) - pow(Qd_ast,0.66667) ;
		z_acute = 2.4 * d * term;
		Lf=3.5*d*pow(Qd_ast, 0.66667);
	} else {
		double term = 1.0 - pow(Qd_ast, 0.4);
		z_acute = 2.4 * d * term;
		Lf=3.5*d*pow(Qd_ast, 0.4);
	}
	
	
	//double Lf = 0.0148 * pow(q,0.4) - 1.02 * d;
	
	if (Lf < hc) {
#ifdef _DEBUG
		opserr << "LocalizedFireSFPE::getFlux() - flame is not impinging ceiling, method has not implemented.\n";
		//exit(-1);
#endif
		return 0;
		}
	
    //Dimensionless HRR
	double termb = 1.11 * 1e6 *pow(hb,2.5);
	double Qhb_ast = q / termb;

	double termc = 1.11 * 1e6 * pow(hc,2.5);
	double Qhc_ast = q / termc;
	


	// now calculate H plus Lh
	//double Lt = 2.9 * h * pow(Qh_ast,0.33);
   double LB = (2.3* pow(Qhb_ast, 0.3) -1)*hb;       // 2.3
	double  LC = (2.9* pow(Qhc_ast, 0.4) -1)*hc;      //2.9
	double LC0 = (2.89 * pow(Qhc_ast, 1.0/3.0) - 1) * hc;

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

	// now calculate w for beam
	double w,q_dot;
	if(FireType ==0){
		//for ceiling
		w = (r + hc + z_acute) / (LC0 +hc + z_acute);
		q_dot = 518.8* exp(-3.7*w)*1000;
	}else if(FireType==1){
		// for downward face of upper flange
		 w = (r+hc+z_acute)/(LC+ hc + z_acute);
		q_dot = 100.5* exp(-2.85*w)*1000;
	}else if(FireType==2){
		// for web of beam
		 w = (r+hc+z_acute)/(LC + hc + z_acute);
		q_dot = 148.1* exp(-2.75*w)*1000;
	}else if(FireType==3){
		// for upward face of lower flange
		w = (r+hc+z_acute)/(LC + hc + z_acute);
		q_dot = 148.1* exp(-2.75*w)*1000;
	}else if(FireType==4){
		// for downward face of lower flange
		w = (r+hb+z_acute)/(LB + hb + z_acute);
		q_dot = 518.8* exp(-3.7*w)*1000;
	}else{
		opserr<<"Localised fireSFPE::Incorrect fire type tag"<<endln;
	}


	// now determine the flux
	
	int tag = node->getTag();
#ifdef _FDEBUG
	opserr<<" NodeTag: "<<node->getTag()<< " r:  "<<r<<"____ "<< " q:  "<<q_dot<<"____ ";
#endif
	return q_dot;
}



double
LocalizedFireSFPE::determineFireConvec(HeatTransferNode* node, double time, int FireType) {
	double q, convech;
	//if Qsquare is activated
	if (Qsquare > 1)
		q = (time / Qsquare) * (time / Qsquare) * 1e6;
	else
		q = ini_q;
	//applying Qsquare
	if (q > ini_q)
		q = ini_q;

	double constant = 1.11 * 1e6 * pow(d, 2.5);
	double Qd_ast = q / constant;

	double z_acute, Lf;

	if (Qd_ast < 1.0) {
		double term = pow(Qd_ast, 0.4) - pow(Qd_ast, 0.66667);
		z_acute = 2.4 * d * term;
		Lf = 3.5 * d * pow(Qd_ast, 0.66667);
	}
	else {
		double term = 1.0 - pow(Qd_ast, 0.4);
		z_acute = 2.4 * d * term;
		Lf = 3.5 * d * pow(Qd_ast, 0.4);
	}


	//double Lf = 0.0148 * pow(q,0.4) - 1.02 * d;

	if (Lf < hc) {
#ifdef _DEBUG
		opserr << "LocalizedFireSFPE::getFlux() - flame is not impinging ceiling, method has not implemented.\n";
		//exit(-1);
#endif
		return 0;
	}

	//Dimensionless HRR
	double termb = 1.11 * 1e6 * d * pow(hb, 1.5);
	double Qhb_ast = q / termb;

	double termc = 1.11 * 1e6 * d * pow(hc, 1.5);
	double Qhc_ast = q / termc;



	// now calculate H plus Lh
	//double Lt = 2.9 * h * pow(Qh_ast,0.33);
	double LB = (2.3 * pow(Qhb_ast, 0.3) - 1) * hb;
	double  LC = (2.9 * pow(Qhc_ast, 0.4) - 1) * hc;

	// now calculate r	
	const Vector& coords = node->getCrds();
	int size = coords.Size();

	double deltaX1, deltaX2, sum, r;

	if (centerLine == 1) {
		deltaX1 = x2 - coords(1);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = x3 - coords(2);
		}
		else {
			// then heat transfer coordinate is 2D
			deltaX2 = x3;
		}
	}
	else if (centerLine == 2) {
		deltaX1 = x1 - coords(0);
		if (size == 3) {
			// then heat transfer coordinate is 3D
			deltaX2 = x3 - coords(2);
		}
		else {
			// then heat transfer coordinate is 2D
			deltaX2 = x3;
		}
	}
	else if (centerLine == 3) {
		deltaX1 = x1 - coords(0);
		deltaX2 = x2 - coords(1);
	}

	sum = deltaX1 * deltaX1 + deltaX2 * deltaX2;
	r = sqrt(sum);

	if (FireType == 0 || FireType == 1 || FireType == 2 || FireType == 3) {
		//for ceiling
		if (r <= LC)
			convech = 35;
		else
			convech = 35;
	}
	else if (FireType == 4) {
		// for downward face of lower flange
		if (r <= LB)
			convech = 35;
		else
			convech = 35;
	}
	else {
		opserr << "Localised fireSFPE::Incorrect fire type tag" << endln;
	}
}

void
LocalizedFireSFPE::applyFluxBC(HeatFluxBC* theFlux, double time)
{
	double ambTemp = 293.15;
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
		//double hconvec = this->determineFireConvec(theNode, time, 4);
		//convec->setParameter(hconvec);
		convec->applyFluxBC(time);
	}
	else if (flux_type == 2) {
		Radiation* rad = (Radiation*)theFlux;
		//double bzm = 5.67e-8;
		//double temp = rad->;
		//double temp = this->getGasTemperature(time);
		//double qir = bzm * pow(ambTemp, 4.0);
		//rad->setIrradiation(qir);
		rad->applyFluxBC(time);
	}
	
	else if (flux_type == 3) 
	{
		PrescribedSurfFlux* pflux = (PrescribedSurfFlux*) theFlux;
        int FireType = pflux->getFireType();
		//int flux_type = pflux->getTypeTag();
		int eleTag = pflux->getElementTag();
		int fTag = pflux->getFaceTag();
		HeatTransferDomain* theDomain = pflux->getDomain();
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
		int size = faceNodes.Size();
		Vector nodalFlux(size);

		for (int i = 0; i < size; i++) {
			int nodTag = faceNodes(i);
			HeatTransferNode* theNode = theDomain->getNode(nodTag);
			if (theNode == 0) {
				opserr << "LocalizedFireSFPE::applyFluxBC() - no node with tag " << nodTag << " exists in the domain";
				exit(-1);
				}
			nodalFlux(i) = this->getFlux(theNode,time, FireType);
			//double r = (theNode->getCrds()(0)) * (theNode->getCrds()(0)) + (theNode->getCrds()(2)-1.8) * (theNode->getCrds()(2)-1.8);
			//opserr << "Flux at node " << nodTag <<", r:"<<sqrt(r)<< " ,q=  " << nodalFlux(i) << endln;
			}

		pflux->setData(nodalFlux);
		pflux->applyFluxBC();
	} else 
	{
			opserr << "LocalizedFireSFPE::applyFluxBC() - incorrect flux type "
				<< flux_type << " provided\n";
			exit(-1);
	}
}

