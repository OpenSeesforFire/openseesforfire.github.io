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
// Written by Liming Jiang (liming.jiang@ed.ac.uk)
//

#include <LineTwo.h>
#include <HeatTransferNode.h>
#include <HeatTransferMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//#include <Renderer.h>
#include <HeatTransferDomain.h>
#include <cmath>
#include <Radiation.h>
#include <Convection.h>
#include <PrescribedSurfFlux.h>

double LineTwo::matrixData[4];
Matrix LineTwo::K(matrixData, 2, 2);
Vector LineTwo::Q(2);
double LineTwo::shp[2];
double LineTwo::pts[2];
double LineTwo::wts[2];
double LineTwo::shpd[2];




using std::list;

LineTwo::LineTwo(int tag, int nd1, int nd2,
				   HeatTransferMaterial& m,
				   bool phaseChange)
:HeatTransferElement(tag), theMaterial(0), connectedExternalNodes(2), Qp(2),
 phaseTransformation(phaseChange)
{
	pts[0] = -0.5773502691896258;
	pts[1] = 0.5773502691896258;
	

	wts[0] = 1.0;
	wts[1] = 1.0;

	//pts2[0] = -0.5773502691896258; 
	//pts2[1] =  0.5773502691896258;
	
	//wts2[0] = 1.0;
	//wts2[1] = 1.0;
  

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new HeatTransferMaterial*[2];
    
	if (theMaterial == 0) {
		opserr << "LineTwo::LineTwo - failed allocate material model pointer\n";
		exit(-1);
		}

	int i;
	for (i = 0; i < 2; i++) {

		// Get copies of the material model for each integration point
		theMaterial[i] = m.getCopy();

		// Check allocation
		if (theMaterial[i] == 0) {
			opserr << "LineTwo::LineTwo -- failed to get a copy of material model\n";
			exit(-1);
			}
		}

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;

    
    for (i=0; i<2; i++)
      theNodes[i] = 0;
}

LineTwo::~LineTwo()
{    
    for (int i = 0; i < 2; i++) {
		if (theMaterial[i])
			delete theMaterial[i];
		}

  // Delete the array of pointers to NDMaterial pointer arrays
	if (theMaterial)
		delete[] theMaterial;
}


int
LineTwo::getNumExternalNodes() const
{
    return 2;
}


const ID&
LineTwo::getExternalNodes()
{
    return connectedExternalNodes;
}


HeatTransferNode**
LineTwo::getNodePtrs(void) 
{
    return theNodes;
}


int
LineTwo::getNumDOF()
{
    return 2;
}


const ID& 
LineTwo::getNodesOnFace(int faceTag)
{
    if (faceTag < 1 || faceTag > 2) {
		opserr << "LineTwo::getNodesOnFace -- improper faceTag: " << faceTag << " for LineTwo\n";
		exit(-1);
		}
	//if (nodeNum < 1 || faceTag > 2) {
	//	opserr << "LineTwo::getNodesOnFace -- improper nodeNum: " << nodeNum << " on the face\n";
	//	exit(-1);
	//	}
	static ID fnodes(1);
	for (int i = 0; i < 1; i++) {
		int localNum = faceTag-1;
		fnodes(i) = connectedExternalNodes(localNum);
		}

	return fnodes;
}


void
LineTwo::setDomain(HeatTransferDomain* theDomain)
{
	// Check HeatTransferDomain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
		}

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
  
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
  

    if (theNodes[0] == 0 || theNodes[1] == 0 ) {
		return;
		}

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
  
    
    if (dofNd1 != 1 || dofNd2 != 1 ) {
		return;
		}

    this->HeatTransferDomainComponent::setDomain(theDomain);
}


int
LineTwo::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
	if ((retVal = this->HeatTransferElement::commitState()) != 0) {
		opserr << "LineTwo::commitState () - failed in base class";
		}    

	// Loop over the integration points and commit the material states
	for (int i = 0; i < 2; i++)
		retVal += theMaterial[i]->commitState();

    return retVal;
}


int
LineTwo::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 2; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}


int
LineTwo::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 2; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
LineTwo::update()
{
	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	
	
	static double T[2];

	T[0] = Temp1(0);
	T[1] = Temp2(0);

	double T_hat;

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < 2; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i]);

		// Interpolate nodal temperatures
		// T_hat = NT;
		T_hat = 0;
		for (int beta = 0; beta < 2; beta++) {
			T_hat += shp[beta] * T[beta];
			}

		// Set the temperature in material model
		ret += theMaterial[i]->setTrialTemperature(T_hat);
		}

	return ret;
}


const Matrix&
LineTwo::getConductionTangent()
{
	K.Zero();

	double dvol;

	// opserr << "the element tag is " << this->getTag() << " \n";
	// Loop over the integration points
	for (int i = 0; i < 2; i++) {
		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i]);
		dvol *= wts[i];
		// Get the material information--heat conductivity
		const Matrix& D = theMaterial[i]->getConductivity();

		// Perform numerical integration
		// K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;

		double D00 = D(0,0); 

		// Form elemental tangent  
		for (int m = 0; m < 2; m++) {
			for (int n = 0; n < 2; n++) {
				K(m,n) += shpd[m]* D00 * shpd[n]  * dvol;
				}
			}
		}   
	
	//opserr << this->getTag() << " conduction tang:" << endln;
	//opserr << K << endln;

    return K;
}


const Matrix&
LineTwo::getCapacityTangent()
{
	K.Zero();
	if (phaseTransformation == true){
		static double enth[2];
		static double rc[2]; // rc = rho times specific heat;

		const Vector& Temp1 = theNodes[0]->getTrialTemperature();
		const Vector& Temp2 = theNodes[1]->getTrialTemperature();

		static double T[2];

		T[0] = Temp1(0);
		T[1] = Temp2(0);


		// getting the nodal enthalpy values
		enth[0] = theMaterial[0]->getEnthalpy(T[0]);
		enth[1] = theMaterial[0]->getEnthalpy(T[1]);

		double dHdx, dTdx;	
		double rcdvol;
 
		// Loop over the integration points
		for (int i = 0; i < 2; i++) {
			// Determine Jacobian for this integration point
			rcdvol = this->shapeFunction(pts[i]);
			dHdx = 0;
			dTdx = 0;

			//double rhoi = theMaterial[i]->getRho(); 
			//rc[i] = rhoi * (theMaterial[i]->getSpecificHeat());

			for (int beta = 0; beta < 2; beta++) {
				dHdx += shpd[beta] * enth[beta];
				
				dTdx += shpd[beta] * T[beta];
				
				}

			if ((dTdx != 0)&& (dHdx != 0)) {
				//double dT = dTdx * dTdx ;

				// using Lemmon's approximation
				//double dH = dHdx * dHdx ;
				rc[i] = dHdx/dTdx;


				//// using DelGuidice's approximation
				//double dH = dHdx * dTdx + dHdy * dTdy;
				//rc[i] = dH/dT;
			} else {
				double rhoi = theMaterial[i]->getRho(); 
				rc[i] = rhoi * (theMaterial[i]->getSpecificHeat()); 
			}	

			rcdvol *= (rc[i] * wts[i]);
			//opserr<<"rcdvol"<<rc[i]<<endln;

			for (int m = 0; m < 2; m++) {
				for (int n = 0; n < 2; n++) {
					K(m,n) += shp[m] * rcdvol * shp[n];
					}
				}
			//end of m,n
			}
		//end of
		} else {
			static double rhoi[2];
			static double cpi[2];
			//opserr << this->getTag() << " capacity tangent\n";

			for (int i = 0; i < 2; i++) {
				// enable to account for temperature dependent density
				rhoi[i] = theMaterial[i]->getRho(); 
				cpi[i] = theMaterial[i]->getSpecificHeat(); // get specific heat
				}  

			double rcdvol;

			// Loop over the integration points
			for (int i = 0; i < 2; i++) {

				// Determine Jacobian for this integration point
				rcdvol = this->shapeFunction(pts[i]);
				rcdvol *= (rhoi[i] * cpi[i] * wts[i]);
				//opserr<<"rcdvol"<<rhoi[i]*cpi[i]<<endln;
				for (int m = 0; m < 2; m++) {
					for (int n = 0; n < 2; n++) {
						K(m,n) += shp[m] * rcdvol * shp[n];
						}
					}
				}
		}

//opserr << this->getTag() << "capacity tang:" << endln;
//	opserr << K << endln;

	return K;
}


const Matrix&
LineTwo::getRadiationTangent()
{
    K.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getEmissionParameter();  // get the product of epsilon and sigma

			double T_hat = (theNodes[tag-1]->getTrialTemperature())(0);
			
			double hr = eps_sgm * 4 * pow(T_hat, 3);

				
			K(tag-1, tag-1) +=  hr;
		}
	rad_iter++;
	}

	//opserr << this->getTag() << " radiation K:" << endln;
	//opserr << K << endln;
	return K;
}


const Matrix& 
LineTwo::getConvectionTangent()
{
    K.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();  // get the covective heat transfer coefficient


			K(tag-1, tag-1) +=  h ;
						
			}
		convec_iter++;
		}
	//opserr << this->getTag() << " convection tang:" << endln;
	//opserr << K << endln;

	return K;
}


void
LineTwo::zeroFlux()
{
   Qp.Zero();
   return;
}


int 
LineTwo::addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor)
{
    Qp.Zero();
    int tag = theFlux->getFaceTag();
	int the_tag = tag - 1;
    const Vector& data = theFlux->getData();
	pFlag = true;

	
	double qk = data(0);

	Qp(the_tag) = qk ;
			
	//opserr << "element: "<<this->getTag()<<"Nodal Data: "<<data<<" , Qp= " << Qp << endln;
		
	return 0;
}


const Vector&
LineTwo::get_Q_Transient()
{
    Q.Zero();

    this->getCapacityTangent();

	const Vector& Tdot1 = theNodes[0]->getTrialTdot();
	const Vector& Tdot2 = theNodes[1]->getTrialTdot();
	

	static double Td[4];

	Td[0] = Tdot1(0);
	Td[1] = Tdot2(0);

    
	for (int m = 0; m < 2; m++) {
		for (int n = 0; n < 2; n++) {
			Q(m) += (-K(m, n) * Td[n]);
			}
		}

	//opserr << this->getTag() << " transient Q:" << endln;
	//opserr << Q << endln;

	return Q;	
}


const Vector&
LineTwo::get_Q_Conduction()
{
	Q.Zero();	
	
	this->getConductionTangent();

	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	
	static double Ttrial[2];

	Ttrial[0] = Temp1(0);
	Ttrial[1] = Temp2(0);
	

	for (int m = 0; m < 2; m++) {
		for (int n = 0; n < 2; n++) {
			Q(m) += (-K(m, n) * Ttrial[n]);
			}
		}
	if (pFlag == true)
		Q.addVector(1.0, Qp, 1.0);

	//opserr << this->getTag() << " conduction Q:" << endln;
	//opserr << Q << endln;
	return Q;	
}


const Vector&
LineTwo::get_Q_Radiation()
{
    Q.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getEmissionParameter();  // get the product of epsilon and sigma

			double alpha = theRadiationBC->getAbsorptivity();
			double g = theRadiationBC->getIrradiation();  // get G

			double T_hat = (theNodes[tag-1]->getTrialTemperature())(0);
			
			double q_ir = alpha * g - eps_sgm * pow(T_hat, 4); 

			Q(tag-1) += q_ir ;
		}
	 rad_iter++;
	}

	//opserr << this->getTag() << " radiation Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


const Vector&
LineTwo::get_Q_Convection()
{
    Q.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;	
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();
			double Ta = theConvectionBC->getSurroundingTemp();
			
			double T_hat = (theNodes[tag-1]->getTrialTemperature())(0);

			double qc = h * (Ta - T_hat);

			Q(tag-1) += qc;
		}		
		convec_iter++;
	}
    //opserr << Q <<endln;

	//opserr << this->getTag() << " convection Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


double 
LineTwo::shapeFunction(double xi)
{
	const Vector& nd1Crds = theNodes[0]->getCrds();
	const Vector& nd2Crds = theNodes[1]->getCrds();

	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp[0] = 0.5 * oneMinusxi;	// N_1
	shp[1] = 0.5 * onePlusxi;		// N_2
	

	double Length;
  
  if(nd1Crds.Size()==2&&nd2Crds.Size()==2){

    Length = sqrt((nd1Crds(0) - nd2Crds(0))*(nd1Crds(0) - nd2Crds(0))+(nd1Crds(1) - nd2Crds(1))*(nd1Crds(1) - nd2Crds(1)));
  }else if(nd1Crds.Size()==3&&nd2Crds.Size()==3){
    Length = sqrt((nd1Crds(0) - nd2Crds(0))*(nd1Crds(0) - nd2Crds(0))+(nd1Crds(1) - nd2Crds(1))*(nd1Crds(1) - nd2Crds(1))+(nd1Crds(2) - nd2Crds(2))*(nd1Crds(2) - nd2Crds(2)) );
  }
  else
  {
      Length = fabs(nd1Crds(0)-nd2Crds(1));
  }

		
	
    shpd[0] = -1/Length;	// N_1,1
    shpd[1] =  1/Length;		// N_2,1


    return Length/2;
}





