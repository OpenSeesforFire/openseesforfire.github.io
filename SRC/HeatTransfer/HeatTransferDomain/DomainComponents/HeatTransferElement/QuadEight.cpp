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

#include <QuadEight.h>
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

double QuadEight::matrixData[64];
Matrix QuadEight::K(matrixData, 8, 8);
Vector QuadEight::Q(8);
double QuadEight::shp[3][8];
double QuadEight::pts[9][2];
double QuadEight::wts[9];
double QuadEight::shp2[3];
double QuadEight::pts2[3];
double QuadEight::wts2[3];
int    QuadEight::npface[4][3];

using std::list;

QuadEight::QuadEight(int tag, int nd1, int nd2, int nd3, int nd4, 
				     int nd5, int nd6, int nd7, int nd8, 
					 HeatTransferMaterial& m,
					 bool phaseChange)
:HeatTransferElement(tag), theMaterial(0), connectedExternalNodes(8), Qp(8),
 phaseTransformation(phaseChange)
{
	pts[0][0] = -0.77459666924148340;
	pts[0][1] = -0.77459666924148340;
	pts[1][0] = -0.77459666924148340;
	pts[1][1] =  0;
	pts[2][0] = -0.77459666924148340;
	pts[2][1] =  0.77459666924148340;
	pts[3][0] =  0;
	pts[3][1] = -0.77459666924148340;
	pts[4][0] =  0;
	pts[4][1] =  0;
	pts[5][0] =  0;
	pts[5][1] =  0.77459666924148340;
	pts[6][0] =  0.77459666924148340;
	pts[6][1] = -0.77459666924148340;
	pts[7][0] =  0.77459666924148340;
	pts[7][1] =  0;
	pts[8][0] =  0.77459666924148340;
	pts[8][1] =  0.77459666924148340;

	wts[0] = 0.30864197530864196;
	wts[1] = 0.49382716049382713;
	wts[2] = 0.30864197530864196;
	wts[3] = 0.49382716049382713;
	wts[4] = 0.79012345679012341;
	wts[5] = 0.49382716049382713;
	wts[6] = 0.30864197530864196;
	wts[7] = 0.49382716049382713;
	wts[8] = 0.30864197530864196;

	pts2[0] = -0.77459666924148340;
	pts2[1] =  0;
	pts2[2] =  0.77459666924148340;
	
	wts2[0] = 0.55555555555555558;
	wts2[1] = 0.88888888888888884;
	wts2[2] = 0.55555555555555558;

	npface[0][0] = 0;  //1
    npface[0][1] = 4;  //5
	npface[0][2] = 1;  //2

	npface[1][0] = 1;  //2
	npface[1][1] = 5;  //6
	npface[1][2] = 2;  //3

	npface[2][0] = 2;  //3
	npface[2][1] = 6;  //7
	npface[2][2] = 3;  //4

	npface[3][0] = 3;  //4
	npface[3][1] = 7;  //8
	npface[3][2] = 0;  //1

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new HeatTransferMaterial*[9];
    
	if (theMaterial == 0) {
		opserr << "QuadEight::QuadEight - failed allocate material model pointer\n";
		exit(-1);
		}

	int i;
	for (i = 0; i < 9; i++) {

		// Get copies of the material model for each integration point
		theMaterial[i] = m.getCopy();

		// Check allocation
		if (theMaterial[i] == 0) {
			opserr << "QuadEight::QuadEight -- failed to get a copy of material model\n";
			exit(-1);
			}
		}

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
	connectedExternalNodes(4) = nd5;
    connectedExternalNodes(5) = nd6;
    connectedExternalNodes(6) = nd7;
    connectedExternalNodes(7) = nd8;
    
    for (i = 0; i < 8; i++)
      theNodes[i] = 0;
}

QuadEight::~QuadEight()
{    
    for (int i = 0; i < 9; i++) {
		if (theMaterial[i])
			delete theMaterial[i];
		}

  // Delete the array of pointers to NDMaterial pointer arrays
	if (theMaterial)
		delete[] theMaterial;
}


int
QuadEight::getNumExternalNodes() const
{
    return 8;
}


const ID&
QuadEight::getExternalNodes()
{
    return connectedExternalNodes;
}


HeatTransferNode**
QuadEight::getNodePtrs(void) 
{
    return theNodes;
}


int
QuadEight::getNumDOF()
{
    return 8;
}


const ID& 
QuadEight::getNodesOnFace(int faceTag)
{
    if (faceTag < 1 || faceTag > 4) {
		opserr << "QuadEight::getNodesOnFace -- improper faceTag: " << faceTag << " for QuadEight\n";
		exit(-1);
		}
	//if (nodeNum < 1 || faceTag > 3) {
	//	opserr << "QuadEight::getNodesOnFace -- improper nodeNum: " << nodeNum << " on the face\n";
	//	exit(-1);
	//	}

	static ID fnodes(3);

	for (int i = 0; i < 3; i++) {
		int localNum = npface[faceTag-1][i];
		fnodes(i) = connectedExternalNodes(localNum);
		}

	return fnodes;
}


void
QuadEight::setDomain(HeatTransferDomain* theDomain)
{
	// Check HeatTransferDomain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		theNodes[4] = 0;
		theNodes[5] = 0;
		theNodes[6] = 0;
		theNodes[7] = 0;
		return;
		}

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);
	int Nd5 = connectedExternalNodes(4);
    int Nd6 = connectedExternalNodes(5);
    int Nd7 = connectedExternalNodes(6);
    int Nd8 = connectedExternalNodes(7);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);
	theNodes[4] = theDomain->getNode(Nd5);
    theNodes[5] = theDomain->getNode(Nd6);
    theNodes[6] = theDomain->getNode(Nd7);
    theNodes[7] = theDomain->getNode(Nd8);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0 ||
		theNodes[4] == 0 || theNodes[5] == 0 || theNodes[6] == 0 || theNodes[7] == 0) {
			return;
		}

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
	int dofNd5 = theNodes[4]->getNumberDOF();
    int dofNd6 = theNodes[5]->getNumberDOF();
    int dofNd7 = theNodes[6]->getNumberDOF();
    int dofNd8 = theNodes[7]->getNumberDOF();
    
    if (dofNd1 != 1 || dofNd2 != 1 || dofNd3 != 1 || dofNd4 != 1 ||
		dofNd5 != 1 || dofNd6 != 1 || dofNd7 != 1 || dofNd8 != 1) {
			return;
		}

    this->HeatTransferDomainComponent::setDomain(theDomain);
}


int
QuadEight::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
	if ((retVal = this->HeatTransferElement::commitState()) != 0) {
		opserr << "QuadEight::commitState () - failed in base class";
		}    

	// Loop over the integration points and commit the material states
	for (int i = 0; i < 9; i++)
		retVal += theMaterial[i]->commitState();

    return retVal;
}


int
QuadEight::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 9; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}


int
QuadEight::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 9; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
QuadEight::update()
{
	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	const Vector& Temp3 = theNodes[2]->getTrialTemperature();
	const Vector& Temp4 = theNodes[3]->getTrialTemperature();
	const Vector& Temp5 = theNodes[4]->getTrialTemperature();
	const Vector& Temp6 = theNodes[5]->getTrialTemperature();
	const Vector& Temp7 = theNodes[6]->getTrialTemperature();
	const Vector& Temp8 = theNodes[7]->getTrialTemperature();
	
	static double T[8];

	T[0] = Temp1(0);
	T[1] = Temp2(0);
	T[2] = Temp3(0);
	T[3] = Temp4(0);
	T[4] = Temp5(0);
	T[5] = Temp6(0);
	T[6] = Temp7(0);
	T[7] = Temp8(0);

	double T_hat;

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < 9; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate nodal temperatures
		// T_hat = NT;
		T_hat = 0;
		for (int beta = 0; beta < 8; beta++) {
			T_hat += shp[2][beta] * T[beta];
			}

		// Set the temperature in material model
		ret += theMaterial[i]->setTrialTemperature(T_hat);
		}

	return ret;
}


const Matrix&
QuadEight::getConductionTangent()
{
	K.Zero();

	double dvol;

	//opserr << "the element tag is " << this->getTag() << " \n";
	// Loop over the integration points
	for (int i = 0; i < 9; i++) {
		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= wts[i];
		// Get the material information--heat conductivity
		const Matrix& D = theMaterial[i]->getConductivity();

		// Perform numerical integration
		// K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;

		double D00 = D(0,0); double D11 = D(1,1);  // Just deal with OrthotropicMaterial 

		// Form elemental tangent  
		for (int m = 0; m < 8; m++) {
			for (int n = 0; n < 8; n++) {
				K(m,n) += (shp[0][m] * D00 * shp[0][n] + shp[1][m] * D11 * shp[1][n]) * dvol;
				}
			}
		}
	
	//opserr << this->getTag() << "conduction tang:" << endln;
	//opserr << K << endln;
    return K;
}


const Matrix&
QuadEight::getCapacityTangent()
{
	K.Zero();
	if (phaseTransformation == true){
		static double enth[8];
		static double rc[9]; // rc = rho times specific heat;
		const Vector& Temp1 = theNodes[0]->getTrialTemperature();
		const Vector& Temp2 = theNodes[1]->getTrialTemperature();
		const Vector& Temp3 = theNodes[2]->getTrialTemperature();
		const Vector& Temp4 = theNodes[3]->getTrialTemperature();
		const Vector& Temp5 = theNodes[4]->getTrialTemperature();
		const Vector& Temp6 = theNodes[5]->getTrialTemperature();
		const Vector& Temp7 = theNodes[6]->getTrialTemperature();
		const Vector& Temp8 = theNodes[7]->getTrialTemperature();

		static double T[8];

		T[0] = Temp1(0);
		T[1] = Temp2(0);
		T[2] = Temp3(0);
		T[3] = Temp4(0);
		T[4] = Temp5(0);
		T[5] = Temp6(0);
		T[6] = Temp7(0);
		T[7] = Temp8(0);

		// getting the nodal enthalpy values
		enth[0] = theMaterial[0]->getEnthalpy(T[0]);
		enth[1] = theMaterial[0]->getEnthalpy(T[1]);
		enth[2] = theMaterial[0]->getEnthalpy(T[2]);
		enth[3] = theMaterial[0]->getEnthalpy(T[3]);
		enth[4] = theMaterial[0]->getEnthalpy(T[4]);
		enth[5] = theMaterial[0]->getEnthalpy(T[5]);
		enth[6] = theMaterial[0]->getEnthalpy(T[6]);
		enth[7] = theMaterial[0]->getEnthalpy(T[7]);

		double dHdx, dHdy, dTdx, dTdy;	
		double rcdvol;

		// Loop over the integration points
		for (int i = 0; i < 9; i++) {
			// Determine Jacobian for this integration point
			rcdvol = this->shapeFunction(pts[i][0], pts[i][1]);
			dHdx = 0;
			dHdy = 0;
			dTdx = 0;
			dTdy = 0;
			for (int beta = 0; beta < 8; beta++) {
				dHdx += shp[0][beta] * enth[beta];
				dHdy += shp[1][beta] * enth[beta];
				dTdx += shp[0][beta] * T[beta];
				dTdy += shp[1][beta] * T[beta];
				}					

			if (((dTdx != 0) || (dTdy !=0)) && ((dHdx != 0) || (dHdy !=0))){
				double dT = dTdx * dTdx + dTdy * dTdy;

				//// using DelGuidice's approximation
				//double dH = dHdx * dTdx + dHdy * dTdy;
				//rc[i] = dH/dT;
				//if (rc[i] < 0)
				//	rc[i] = (theMaterial[i]->getRho()) * (theMaterial[i]->getSpecificHeat());

				//using Morgan & Lemmon's approximation
				double dH = dHdx * dHdx + dHdy * dHdy;
				rc[i] = sqrt(dH/dT);
				} else {
					double rhoi = theMaterial[i]->getRho(); 
					rc[i] = rhoi * (theMaterial[i]->getSpecificHeat()); 
				}	

			rcdvol *= (rc[i] * wts[i]);

			for (int m = 0; m < 8; m++) {
				for (int n = 0; n < 8; n++) {
					K(m,n) += shp[2][m] * rcdvol * shp[2][n];
					}
				}
			}
		} else {
			static double rhoi[9];
			static double cpi[9];

			//for (int i = 0; i < 9; i++) {
			//	rhoi[i] = theMaterial[i]->getRho(); // enable to account for temperature dependent density
			//	cpi[i] = theMaterial[i]->getSpecificHeat(); // get specific heat
			//	}

			double rcdvol;

			// Loop over the integration points
			for (int i = 0; i < 9; i++) {

				rhoi[i] = theMaterial[i]->getRho(); // enable to account for temperature dependent density
				cpi[i] = theMaterial[i]->getSpecificHeat(); // get specific heat

				// Determine Jacobian for this integration point
				rcdvol = this->shapeFunction(pts[i][0], pts[i][1]);
				rcdvol *= (rhoi[i] * cpi[i] * wts[i]);

				for (int m = 0; m < 8; m++) {
					for (int n = 0; n < 8; n++) {
						K(m,n) += shp[2][m] * rcdvol * shp[2][n];
						}
					}
				}	
		}

	//opserr << this->getTag() << "capacity tang:" << endln;
	//opserr << K << endln;
	return K;
}


const Matrix&
QuadEight::getRadiationTangent()
{
    K.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getEmissionParameter();  // get the product of epsilon and sigma
			double ds;

			for (int i = 0; i < 3; i++) {
				ds = this->shapeFunction(pts2[i],tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double hr = eps_sgm * 4 * pow(T_hat, 3);

				for (int m = 0; m < 3; m++) {
					int ftag = tag-1;
					int mm = npface[ftag][m];
					for (int n = 0; n < 3; n++) {
						int nn = npface[ftag][n];
						K(mm, nn) += shp2[m] * hr * shp2[n] * ds;
						}
					}
				}
			}
		rad_iter++;
		}

	//opserr << this->getTag() << " radiation K:" << endln;
	//opserr << K << endln;
	return K;
}


const Matrix& 
QuadEight::getConvectionTangent()
{
    K.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;
			int tag = theConvectionBC->getFaceTag();
			int ftag = tag -1;
			double h = theConvectionBC->getParameter();  // get the covective heat transfer coefficient

			double ds;

			for (int i = 0; i < 3; i++){
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				for (int m = 0; m < 3; m++) {
					int mm = npface[ftag][m];
					for (int n = 0; n < 3; n++) {
						int nn = npface[ftag][n];
						K(mm, nn) += shp2[m] * h * shp2[n] * ds;
						}
					}
				}
			}
		convec_iter++;
		}
	//opserr << this->getTag() << " convection tang:" << endln;
	//opserr << K << endln;

	return K;
}


void
QuadEight::zeroFlux()
{
   Qp.Zero();
   return;
}


int 
QuadEight::addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor)
{
    int tag = theFlux->getFaceTag();
    Vector data = theFlux->getData();
    double ds;
	pFlag = true;

	// loop over quadrature points
	for (int i = 0; i < 3; i++){
		ds = this->shapeFunction(pts2[i], tag);
		ds *= wts2[i];
		double qk = 0;

		// loop over face nodes
		for (int j = 0; j < 3; j++) {
			qk += data(j) * shp2[j];
			}

		for (int m = 0; m < 3; m++) {
			int mm = npface[tag][m];
			Qp(mm) = shp2[m] * qk * ds;
			}
		}

	return 0;
}


const Vector&
QuadEight::get_Q_Transient()
{
    Q.Zero();

    this->getCapacityTangent();

	const Vector& Tdot1 = theNodes[0]->getTrialTdot();
	const Vector& Tdot2 = theNodes[1]->getTrialTdot();
	const Vector& Tdot3 = theNodes[2]->getTrialTdot();
	const Vector& Tdot4 = theNodes[3]->getTrialTdot();
	const Vector& Tdot5 = theNodes[4]->getTrialTdot();
	const Vector& Tdot6 = theNodes[5]->getTrialTdot();
	const Vector& Tdot7 = theNodes[6]->getTrialTdot();
	const Vector& Tdot8 = theNodes[7]->getTrialTdot();

	static double Td[8];

	Td[0] = Tdot1(0);
	Td[1] = Tdot2(0);
	Td[2] = Tdot3(0);
	Td[3] = Tdot4(0);
	Td[4] = Tdot5(0);
	Td[5] = Tdot6(0);
	Td[6] = Tdot7(0);
	Td[7] = Tdot8(0);

	for (int m = 0; m < 8; m++) {
		for (int n = 0; n < 8; n++) {
			Q(m) += (-K(m, n) * Td[n]);
			}
		}
	//opserr << this->getTag() << "capacity tangent:\n" << K << "\n";
	//opserr << this->getTag() << " capacity residual:\n" << Q << "\n";

	return Q;	
}


const Vector&
QuadEight::get_Q_Conduction()
{
	Q.Zero();	
	
	this->getConductionTangent();

	const Vector& Temp1 = theNodes[0]->getTrialTemperature();
	const Vector& Temp2 = theNodes[1]->getTrialTemperature();
	const Vector& Temp3 = theNodes[2]->getTrialTemperature();
	const Vector& Temp4 = theNodes[3]->getTrialTemperature();
	const Vector& Temp5 = theNodes[4]->getTrialTemperature();
	const Vector& Temp6 = theNodes[5]->getTrialTemperature();
	const Vector& Temp7 = theNodes[6]->getTrialTemperature();
	const Vector& Temp8 = theNodes[7]->getTrialTemperature();

	static double Ttrial[8];

	Ttrial[0] = Temp1(0);
	Ttrial[1] = Temp2(0);
	Ttrial[2] = Temp3(0);
	Ttrial[3] = Temp4(0);
	Ttrial[4] = Temp5(0);
	Ttrial[5] = Temp6(0);
	Ttrial[6] = Temp7(0);
	Ttrial[7] = Temp8(0);

	for (int m = 0; m < 8; m++) {
		for (int n = 0; n < 8; n++) {
			Q(m) += (-K(m, n) * Ttrial[n]);
			}
		}

	if (pFlag == true)
		Q.addVector(1.0, Qp, 1.0); 
	
	//opserr << this->getTag() << "conduction tangent:\n" << K << "\n";
	//opserr << this->getTag() << " conduction residual:\n" << Q << "\n";
	
	return Q;	
}


const Vector&
QuadEight::get_Q_Radiation()
{
    Q.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			int ftag = tag -1;
			double eps_sgm = theRadiationBC->getEmissionParameter();  // get the product of epsilon and sigma

			double alpha = theRadiationBC->getAbsorptivity();
			double g = theRadiationBC->getIrradiation();  // get G

			double ds;

			for (int i = 0; i < 3; i++) {
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double q_ir = alpha * g - eps_sgm * pow(T_hat, 4); 

				for (int m = 0; m < 3; m++) {
					int mm = npface[ftag][m];
					Q(mm) += shp2[m] * q_ir * ds;
					}
				}
			}
		rad_iter++;
		}

	//opserr << this->getTag() << " radiation Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


const Vector&
QuadEight::get_Q_Convection()
{
    Q.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;	
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();
			double Ta = theConvectionBC->getSurroundingTemp();

			double ds;

			for (int i = 0; i < 3; i++) {
				ds = this->shapeFunction(pts2[i], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double qc = h * (Ta - T_hat);

				for (int m = 0; m < 3; m++) {
					int mm = npface[tag-1][m];
					Q(mm) += shp2[m] * qc * ds;
					}
				}
			}
		convec_iter++;
		}
    //opserr << Q <<endln;

	//opserr << this->getTag() << " convection Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


double 
QuadEight::shapeFunction(double xi, double eta)
{
	const Vector& nd1Crds = theNodes[0]->getCrds();
	const Vector& nd2Crds = theNodes[1]->getCrds();
	const Vector& nd3Crds = theNodes[2]->getCrds();
	const Vector& nd4Crds = theNodes[3]->getCrds();

	const Vector& nd5Crds = theNodes[4]->getCrds();
	const Vector& nd6Crds = theNodes[5]->getCrds();
	const Vector& nd7Crds = theNodes[6]->getCrds();
	const Vector& nd8Crds = theNodes[7]->getCrds();

	double oneMinuseta = 1.0 - eta;
	double onePluseta = 1.0 + eta;
	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;
	//double oneMinusetaSq = 1.0 - eta * eta;
	//double oneMinusxiSq = 1.0 - xi * xi;
	double twoTimesxiPluseta = 2 * xi + eta;
	double twoTimesxiMinuseta = 2 * xi - eta;
	double twoTimesetaPlusxi = 2 * eta + xi;
	double twoTimesetaMinusxi = 2 * eta - xi;

	shp[2][0] = 0.25 * oneMinusxi * oneMinuseta * (-onePlusxi - eta);	// N_1
	shp[2][1] = 0.25 * onePlusxi * oneMinuseta * (-oneMinusxi - eta);	// N_2
	shp[2][2] = 0.25 * onePlusxi * onePluseta * (-oneMinusxi + eta);	// N_3
	shp[2][3] = 0.25 * oneMinusxi * onePluseta * (-onePlusxi + eta);	// N_4
	shp[2][4] = 0.5 * oneMinusxi * onePlusxi * oneMinuseta;	// N_5
	shp[2][5] = 0.5 * onePlusxi * oneMinuseta * onePluseta;	// N_6
	shp[2][6] = 0.5 * oneMinusxi * onePlusxi * onePluseta;	// N_7
	shp[2][7] = 0.5 * oneMinusxi * oneMinuseta * onePluseta;	// N_8

	double D[2][8];
	D[0][0] = 0.25 * oneMinuseta * twoTimesxiPluseta; // dN/dxi
	D[0][1] = 0.25 * oneMinuseta * twoTimesxiMinuseta;
	D[0][2] = 0.25 * onePluseta * twoTimesxiPluseta;
	D[0][3] = 0.25 * onePluseta * twoTimesxiMinuseta;
	D[0][4] = -xi * oneMinuseta ;
	D[0][5] = 0.5 * oneMinuseta * onePluseta;
	D[0][6] = -xi * onePluseta;
	D[0][7] = -D[0][5];

	D[1][0] = 0.25 * oneMinusxi * twoTimesetaPlusxi;  // dN/deta
	D[1][1] = 0.25 * onePlusxi * twoTimesetaMinusxi;
	D[1][2] = 0.25 * onePlusxi * twoTimesetaPlusxi;
	D[1][3] = 0.25 * oneMinusxi * twoTimesetaMinusxi;
	D[1][4] = -0.5 * oneMinusxi * onePlusxi;
	D[1][5] = -eta * onePlusxi;
	D[1][6] = -D[1][4];
	D[1][7] = -eta * oneMinusxi;


	double J[2][2];

	J[0][0] = D[0][0] * nd1Crds(0) + D[0][1] * nd2Crds(0) + D[0][2] * nd3Crds(0) + D[0][3] * nd4Crds(0) +
			  D[0][4] * nd5Crds(0) + D[0][5] * nd6Crds(0) + D[0][6] * nd7Crds(0) + D[0][7] * nd8Crds(0);
	J[0][1] = D[0][0] * nd1Crds(1) + D[0][1] * nd2Crds(1) + D[0][2] * nd3Crds(1) + D[0][3] * nd4Crds(1) +
			  D[0][4] * nd5Crds(1) + D[0][5] * nd6Crds(1) + D[0][6] * nd7Crds(1) + D[0][7] * nd8Crds(1);
	J[1][0] = D[1][0] * nd1Crds(0) + D[1][1] * nd2Crds(0) + D[1][2] * nd3Crds(0) + D[1][3] * nd4Crds(0) +
			  D[1][4] * nd5Crds(0) + D[1][5] * nd6Crds(0) + D[1][6] * nd7Crds(0) + D[1][7] * nd8Crds(0);
	J[1][1] = D[1][0] * nd1Crds(1) + D[1][1] * nd2Crds(1) + D[1][2] * nd3Crds(1) + D[1][3] * nd4Crds(1) +
			  D[1][4] * nd5Crds(1) + D[1][5] * nd6Crds(1) + D[1][6] * nd7Crds(1) + D[1][7] * nd8Crds(1);

	double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	
	if (detJ < 0) {
		opserr << "QuadEight::shapeFunction(xi,eta) -- detJ < 0\n";
		exit(-1);
		}

	double oneOverdetJ = 1.0/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1] * oneOverdetJ;
	L[1][0] = -J[0][1] * oneOverdetJ;
	L[0][1] = -J[1][0] * oneOverdetJ;
	L[1][1] =  J[0][0] * oneOverdetJ;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = L[0][0] * D[0][0] + L[0][1] * D[1][0];	// N_1,1
    shp[0][1] = L[0][0] * D[0][1] + L[0][1] * D[1][1];	// N_2,1
    shp[0][2] = L[0][0] * D[0][2] + L[0][1] * D[1][2];	// N_3,1
    shp[0][3] = L[0][0] * D[0][3] + L[0][1] * D[1][3];	// N_4,1
	shp[0][4] = L[0][0] * D[0][4] + L[0][1] * D[1][4];	// N_5,1
    shp[0][5] = L[0][0] * D[0][5] + L[0][1] * D[1][5];	// N_6,1
    shp[0][6] = L[0][0] * D[0][6] + L[0][1] * D[1][6];	// N_7,1
    shp[0][7] = L[0][0] * D[0][7] + L[0][1] * D[1][7];	// N_8,1
	
    shp[1][0] = L[1][0] * D[0][0] + L[1][1] * D[1][0];	// N_1,2
    shp[1][1] = L[1][0] * D[0][1] + L[1][1] * D[1][1];	// N_2,2
    shp[1][2] = L[1][0] * D[0][2] + L[1][1] * D[1][2];	// N_3,2
    shp[1][3] = L[1][0] * D[0][3] + L[1][1] * D[1][3];	// N_4,2
	shp[1][4] = L[1][0] * D[0][4] + L[1][1] * D[1][4];	// N_5,2
    shp[1][5] = L[1][0] * D[0][5] + L[1][1] * D[1][5];	// N_6,2
    shp[1][6] = L[1][0] * D[0][6] + L[1][1] * D[1][6];	// N_7,2
    shp[1][7] = L[1][0] * D[0][7] + L[1][1] * D[1][7];	// N_8,2

    return detJ;
}


double 
QuadEight::shapeFunction(double xi, int faceTag)
{
	if (faceTag < 1 || faceTag > 4) {
		opserr << "QuadEight::shapeFunction -- improper faceTag: " << faceTag << "for QuadEight\n";
		exit(-1);
		}

	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp2[0] = -0.5 * xi * oneMinusxi;	// N_1
	shp2[1] = oneMinusxi * onePlusxi;	// N_2
	shp2[2] = 0.5 * xi * onePlusxi;	    // N_3
	
	int tag = faceTag - 1;
	int node1 = npface[tag][0];
	int node2 = npface[tag][1];
	int node3 = npface[tag][2];

	const Vector& nd1Crds = theNodes[node1]->getCrds();
	const Vector& nd2Crds = theNodes[node2]->getCrds();
	const Vector& nd3Crds = theNodes[node3]->getCrds();

	double J00, J01;
	double xiMinusHalfOne = xi - 0.5;
	double xiPlusHalfOne = xi + 0.5;

	J00 = xiMinusHalfOne * nd1Crds(0) - 2 * xi * nd2Crds(0) + xiPlusHalfOne * nd3Crds(0);  // dx/dxi
	J01 = xiMinusHalfOne * nd1Crds(1) - 2 * xi * nd2Crds(1) + xiPlusHalfOne * nd3Crds(1);  // dy/dxi

	double JGamma = sqrt(pow(J00, 2) + pow(J01, 2));

	if (JGamma < 0){
		opserr << "QuadEight::shapeFunction(xi,faceTag) -- JGamma < 0 ";
		exit(-1);
		}
  
    return JGamma;  // return the bounday Jacobian
}


double 
QuadEight::getInterpolatedTemp(int tag)
{
	int ftag = tag -1;
	int node1 = npface[ftag][0];
	int node2 = npface[ftag][1];
	int node3 = npface[ftag][2];

	const Vector& Temp1 = theNodes[node1]->getTrialTemperature();
	const Vector& Temp2 = theNodes[node2]->getTrialTemperature();
	const Vector& Temp3 = theNodes[node3]->getTrialTemperature();

	double T_hat;
	T_hat = shp2[0] * Temp1(0) + shp2[1] * Temp2(0) + shp2[2] * Temp3(0);

	return T_hat;
}


