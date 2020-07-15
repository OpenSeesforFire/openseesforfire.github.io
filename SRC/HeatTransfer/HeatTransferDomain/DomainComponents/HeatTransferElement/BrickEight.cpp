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

#include <BrickEight.h>
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

const double  BrickEight::root3 = sqrt(3.0) ;
const double  BrickEight::one_over_root3 = 1.0 / root3 ;

double BrickEight::matrixData[64];
Matrix BrickEight::K(matrixData, 8, 8);
Vector BrickEight::Q(8);
double BrickEight::shp[4][8];
double BrickEight::pts[8][3] = {
	{-one_over_root3, -one_over_root3, -one_over_root3},
	{one_over_root3, -one_over_root3, -one_over_root3},
	{one_over_root3, one_over_root3, -one_over_root3},
	{-one_over_root3, one_over_root3, -one_over_root3},
	{-one_over_root3, -one_over_root3, one_over_root3},
	{one_over_root3, -one_over_root3, one_over_root3},
	{one_over_root3, one_over_root3, one_over_root3},
	{-one_over_root3, one_over_root3, one_over_root3},
	};
const double BrickEight::wts[] = { 1.0, 1.0, 1.0, 1.0, 
                                   1.0, 1.0, 1.0, 1.0};
double BrickEight::shp2[4];
double BrickEight::pts2[4][2] = {
	{-one_over_root3, -one_over_root3},
	{one_over_root3, -one_over_root3},
	{one_over_root3, one_over_root3},
	{-one_over_root3, one_over_root3}
	};
const double BrickEight::wts2[] = {1.0, 1.0, 1.0, 1.0};
const int BrickEight::npface[6][4] = {
	{0, 1, 2, 3},
	{4, 7, 6, 5},
	{0, 4, 5, 1},
	{1, 5, 6, 2},
	{2, 6, 7, 3},
	{3, 7, 4, 0},
	};

//pts[0][0] = -one_over_root3;
//pts[0][1] = -one_over_root3;
//pts[0][2] = -one_over_root3;
//pts[1][0] =  one_over_root3;
//pts[1][1] = -one_over_root3;
//pts[1][2] = -one_over_root3;
//pts[2][0] =  one_over_root3;
//pts[2][1] =  one_over_root3;
//pts[2][2] = -one_over_root3;
//pts[3][0] = -one_over_root3;
//pts[3][1] =  one_over_root3;
//pts[3][2] = -one_over_root3;
//pts[4][0] = -one_over_root3;
//pts[4][1] = -one_over_root3;
//pts[4][2] =  one_over_root3;
//pts[5][0] =  one_over_root3;
//pts[5][1] = -one_over_root3;
//pts[5][2] =  one_over_root3;
//pts[6][0] =  one_over_root3;
//pts[6][1] =  one_over_root3;
//pts[6][2] =  one_over_root3;
//pts[7][0] = -one_over_root3;
//pts[7][1] =  one_over_root3;
//pts[7][2] =  one_over_root3;

//npface[0][0] = 0;
//npface[0][1] = 1;
//npface[0][2] = 2;
//npface[0][3] = 3;
//npface[1][0] = 4;
//npface[1][1] = 7;
//npface[1][2] = 6;
//npface[1][3] = 5;
//npface[2][0] = 0;
//npface[2][1] = 4;
//npface[2][2] = 5;
//npface[2][3] = 1;
//npface[3][0] = 1;
//npface[3][1] = 5;
//npface[3][2] = 6;
//npface[3][3] = 2;
//npface[4][0] = 2;
//npface[4][1] = 6;
//npface[4][2] = 7;
//npface[4][3] = 3;
//npface[5][0] = 3;
//npface[5][1] = 7;
//npface[5][2] = 4;
//npface[5][3] = 0;
//
//pts2[0][0] = -one_over_root3;
//pts2[0][1] = -one_over_root3;
//pts2[1][0] =  one_over_root3;
//pts2[1][1] = -one_over_root3;
//pts2[2][0] =  one_over_root3;
//pts2[2][1] =  one_over_root3;
//pts2[3][0] = -one_over_root3;
//pts2[3][1] =  one_over_root3;

using std::list;

BrickEight::BrickEight(int tag, int nd1, int nd2, int nd3, int nd4, 
					   int nd5, int nd6, int nd7, int nd8,
				       HeatTransferMaterial& m,
				       bool phaseChange)
:HeatTransferElement(tag), theMaterial(0), connectedExternalNodes(8), Qp(8),
 phaseTransformation(phaseChange)
{
    // Allocate arrays of pointers to NDMaterials
    theMaterial = new HeatTransferMaterial*[8];
    
	if (theMaterial == 0) {
		opserr << "BrickEight::BrickEight - failed allocate material model pointer\n";
		exit(-1);
		}

	int i;
	for (i = 0; i < 8; i++) {
		// Get copies of the material model for each integration point
		theMaterial[i] = m.getCopy();
		// Check allocation
		if (theMaterial[i] == 0) {
			opserr << "BrickEight::BrickEight -- failed to get a copy of material model\n";
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


BrickEight::~BrickEight()
{    
    for (int i = 0; i < 8; i++) {
		if (theMaterial[i])
			delete theMaterial[i];
		}

  // Delete the array of pointers to NDMaterial pointer arrays
	if (theMaterial)
		delete[] theMaterial;
}


int
BrickEight::getNumExternalNodes() const
{
    return 8;
}


const ID&
BrickEight::getExternalNodes()
{
    return connectedExternalNodes;
}


HeatTransferNode**
BrickEight::getNodePtrs(void) 
{
    return theNodes;
}


int
BrickEight::getNumDOF()
{
    return 8;
}


const ID& 
BrickEight::getNodesOnFace(int faceTag)
{
    if (faceTag < 1 || faceTag > 8) {
		opserr << "BrickEight::getNodesOnFace -- improper faceTag: " << faceTag << " for BrickEight\n";
		exit(-1);
		}
	//if (nodeNum < 1 || faceTag > 2) {
	//	opserr << "QuadFour::getNodesOnFace -- improper nodeNum: " << nodeNum << " on the face\n";
	//	exit(-1);
	//	}
	static ID fnodes(4);
	for (int i = 0; i < 4; i++) {
		int localNum = npface[faceTag-1][i];
		fnodes(i) = connectedExternalNodes(localNum);
		}

	return fnodes;
}


void
BrickEight::setDomain(HeatTransferDomain* theDomain)
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

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0
		|| theNodes[4] == 0 || theNodes[5] == 0 || theNodes[6] == 0 || theNodes[7] == 0
		) {
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
    
    if (dofNd1 != 1 || dofNd2 != 1 || dofNd3 != 1 || dofNd4 != 1
		|| dofNd5 != 1 || dofNd6 != 1 || dofNd7 != 1 || dofNd8 != 1
		) {
			return;
		}

    this->HeatTransferDomainComponent::setDomain(theDomain);
}


int
BrickEight::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
	if ((retVal = this->HeatTransferElement::commitState()) != 0) {
		opserr << "BrickEight::commitState () - failed in base class";
		}    

	// Loop over the integration points and commit the material states
	for (int i = 0; i < 8; i++)
		retVal += theMaterial[i]->commitState();

    return retVal;
}


int
BrickEight::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 8; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}


int
BrickEight::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 8; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
BrickEight::update()
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
	for (int i = 0; i < 8; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1], pts[i][2]);

		// Interpolate nodal temperatures
		// T_hat = NT;
		T_hat = 0;
		for (int beta = 0; beta < 8; beta++) {
			T_hat += shp[3][beta] * T[beta];
			}

		// Set the temperature in material model
		ret += theMaterial[i]->setTrialTemperature(T_hat);
		}

	return ret;
}


const Matrix&
BrickEight::getConductionTangent()
{
	K.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < 8; i++) {
		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1], pts[i][2]);
		dvol *= wts[i];
		// Get the material information--heat conductivity
		const Matrix& D = theMaterial[i]->getConductivity();

		// Perform numerical integration
		// K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;

		// Just deal with OrthotropicMaterial
		double D00 = D(0,0); 
		double D11 = D(1,1);
		double D22 = D(2,2);
 
		// Form elemental tangent  
		for (int m = 0; m < 8; m++) {
			for (int n = 0; n < 8; n++) {
				K(m,n) += (shp[0][m] * D00 * shp[0][n] + shp[1][m] * D11 * shp[1][n] 
				+ shp[2][m] * D22 * shp[2][n]) * dvol;
				}
			}
		}   
	
	//opserr << this->getTag() << "conguction tang:" << endln;
	//opserr << K << endln;
	//opserr << "capacity tangent:\n" << K << "\n";

    return K;
}


const Matrix&
BrickEight::getCapacityTangent()
{
	K.Zero();
	if (phaseTransformation == true){
		static double enth[8];
		static double rc[8]; // rc = rho times specific heat;

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

		double dHdx, dHdy, dHdz, dTdx, dTdy, dTdz;	
		double rcdvol;
 
		// Loop over the integration points
		for (int i = 0; i < 8; i++) {
			// Determine Jacobian for this integration point
			rcdvol = this->shapeFunction(pts[i][0], pts[i][1], pts[i][2]);
			dHdx = 0.0;
			dHdy = 0.0;
			dHdz = 0.0;
			dTdx = 0.0;
			dTdy = 0.0;
			dTdz = 0.0;

			//double rhoi = theMaterial[i]->getRho(); 
			//rc[i] = rhoi * (theMaterial[i]->getSpecificHeat());

			for (int beta = 0; beta < 8; beta++) {
				dHdx += shp[0][beta] * enth[beta];
				dHdy += shp[1][beta] * enth[beta];
				dHdz += shp[2][beta] * enth[beta];
				dTdx += shp[0][beta] * T[beta];
				dTdy += shp[1][beta] * T[beta];
				dTdz += shp[2][beta] * T[beta];
				}

			if (((dTdx != 0) || (dTdy !=0) || (dTdz != 0)) && ((dHdx != 0) || (dHdy !=0) || (dHdz != 0))) {
				double dT = dTdx * dTdx + dTdy * dTdy + dTdz * dTdz;
				// using Lemmon's approximation
				double dH = dHdx * dHdx + dHdy * dHdy + dHdz * dHdz;
				rc[i] = sqrt(dH/dT);

				//// using DelGuidice's approximation
				//double dH = dHdx * dTdx + dHdy * dTdy;
				//rc[i] = dH/dT;
				} else {
					double rhoi = theMaterial[i]->getRho(); 
					rc[i] = rhoi * (theMaterial[i]->getSpecificHeat()); 
				}	

			rcdvol *= (rc[i] * wts[i]);

			for (int m = 0; m < 8; m++) {
				for (int n = 0; n < 8; n++) {
					K(m,n) += shp[3][m] * rcdvol * shp[3][n];
					}
				}
			}
		} else {
			static double rhoi[8];
			static double cpi[8];
			//opserr << this->getTag() << " capacity tangent\n";

			for (int i = 0; i < 8; i++) {
				// enable to account for temperature dependent density
				rhoi[i] = theMaterial[i]->getRho(); 
				cpi[i] = theMaterial[i]->getSpecificHeat(); // get specific heat
				}

			double rcdvol;

			// Loop over the integration points
			for (int i = 0; i < 8; i++) {

				// Determine Jacobian for this integration point
				rcdvol = this->shapeFunction(pts[i][0], pts[i][1], pts[i][2]);
				rcdvol *= (rhoi[i] * cpi[i] * wts[i]);

				for (int m = 0; m < 8; m++) {
					for (int n = 0; n < 8; n++) {
						K(m,n) += shp[3][m] * rcdvol * shp[3][n];
						}
					}
				}
		}

	//opserr << "capacity tangent:\n" << K << "\n";

	return K;
}


const Matrix&
BrickEight::getRadiationTangent()
{
    K.Zero();
	list<Radiation*>::iterator rad_iter = theRadiationBCs.begin();

	while (rad_iter != theRadiationBCs.end()) {
		if (*rad_iter != 0){
			Radiation* theRadiationBC = *rad_iter;
			int tag = theRadiationBC->getFaceTag();
			double eps_sgm = theRadiationBC->getEmissionParameter();  // get the product of epsilon and sigma
			double ds;

			for (int i = 0; i < 4; i++) {
				ds = this->shapeFunction(pts2[i][0], pts2[i][1], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double hr = eps_sgm * 4 * pow(T_hat, 3);

				for (int m = 0; m < 4; m++) {
					int mm = npface[tag-1][m];
					for (int n = 0; n < 4; n++) {
						int nn = npface[tag-1][n];
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
BrickEight::getConvectionTangent()
{
    K.Zero();
	list<Convection*>::iterator convec_iter = theConvectionBCs.begin();

	while (convec_iter != theConvectionBCs.end()) {
		if (*convec_iter != 0) {
			Convection* theConvectionBC = *convec_iter;
			int tag = theConvectionBC->getFaceTag();
			double h = theConvectionBC->getParameter();  // get the covective heat transfer coefficient

			double ds;

			for (int i = 0; i < 4; i++){
				ds = this->shapeFunction(pts2[i][0], pts2[i][1], tag);
				ds *= wts2[i];

				for (int m = 0; m < 4; m++) {
					int mm = npface[tag-1][m];
					for (int n = 0; n < 4; n++) {
						int nn = npface[tag-1][n];
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
BrickEight::zeroFlux()
{ 
   Qp.Zero();
   return;
}


int 
BrickEight::addPrecribedSurfFlux(PrescribedSurfFlux* theFlux, double factor)
{
    Qp.Zero();
    int tag = theFlux->getFaceTag();
	int the_tag = tag - 1;
    const Vector& data = theFlux->getData();
    double ds;
	pFlag = true;

	// loop over quadrature points
	for (int i = 0; i < 4; i++){
		ds = this->shapeFunction(pts2[i][0], pts2[i][1], tag);
		ds *= wts2[i];
		double qk = 0;

		// loop over face nodes
		for (int j = 0; j < 4; j++) {
			qk += data(j) * shp2[j];
			}

		for (int m = 0; m < 4; m++) {
			int mm = npface[the_tag][m];
			Qp(mm) += shp2[m] * qk * ds;
			}
		}
   // opserr << "element: "<<this->getTag()<<"Nodal Data: "<<data<<" , Qp= " << Qp << endln;
	return 0;
}


const Vector&
BrickEight::get_Q_Transient()
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
	//opserr << "capacity tangent:\n" << K << "\n";
	//opserr <<"Element: "<<  this->getTag() << " capacity residual:\n" << Q << "\n";

	return Q;	
}


const Vector&
BrickEight::get_Q_Conduction()
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

	//opserr << "conduction tangent:\n" << K << "\n";
	//opserr <<"Element: "<<  this->getTag() << " conduction residual:\n" << Q << "\n";
	return Q;	
}


const Vector&
BrickEight::get_Q_Radiation()
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

			double ds;

			for (int i = 0; i < 4; i++) {
				ds = this->shapeFunction(pts2[i][0], pts2[i][1], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double q_ir = alpha * g - eps_sgm * pow(T_hat, 4); 

				for (int m = 0; m < 4; m++) {
					int mm = npface[tag-1][m];
					Q(mm) += shp2[m] * q_ir * ds;
					}
				}
			}
		rad_iter++;
		}

	//opserr <<"Element: "<< this->getTag() << " radiation Q:" << endln;
	//opserr << Q << endln;
	return Q;
}


const Vector&
BrickEight::get_Q_Convection()
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

			for (int i = 0; i < 4; i++) {
				ds = this->shapeFunction(pts2[i][0], pts2[i][1], tag);
				ds *= wts2[i];

				double T_hat = this->getInterpolatedTemp(tag);
				double qc = h * (Ta - T_hat);

				for (int m = 0; m < 4; m++) {
					int mm = npface[tag-1][m];
					Q(mm) += shp2[m] * qc * ds;
					}
				}
			}
		convec_iter++;
		}
    //opserr << Q <<endln;

	//opserr <<"Element: "<<  this->getTag() << " convection Q:" <<Q<< endln;
	//opserr << Q << endln;
	return Q;
}


double 
BrickEight::shapeFunction(double xi, double eta, double zeta)
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
	double oneMinuszeta = 1.0 - zeta;
	double onePluszeta = 1.0 + zeta;

	shp[3][0] = 0.125 * oneMinusxi * oneMinuseta * oneMinuszeta;	// N1
	shp[3][1] = 0.125 * onePlusxi * oneMinuseta * oneMinuszeta;		// N2
	shp[3][2] = 0.125 * onePlusxi * onePluseta * oneMinuszeta;		// N3
	shp[3][3] = 0.125 * oneMinusxi * onePluseta * oneMinuszeta;		// N4
	shp[3][4] = 0.125 * oneMinusxi * oneMinuseta * onePluszeta;	    // N5
	shp[3][5] = 0.125 * onePlusxi * oneMinuseta * onePluszeta;		// N6
	shp[3][6] = 0.125 * onePlusxi * onePluseta * onePluszeta;		// N7
	shp[3][7] = 0.125 * oneMinusxi * onePluseta * onePluszeta;		// N8
    
	// Derivative matrix in reference coordinate system
	double D[3][8];
	D[0][0] = -0.125 * oneMinuseta * oneMinuszeta;
	D[0][1] = 0.125 * oneMinuseta * oneMinuszeta;
	D[0][2] = 0.125 * onePluseta * oneMinuszeta;
	D[0][3] = -0.125 * onePluseta * oneMinuszeta;
	D[0][4] = -0.125 * oneMinuseta * onePluszeta;
	D[0][5] = 0.125 * oneMinuseta * onePluszeta;
	D[0][6] = 0.125 * onePluseta * onePluszeta;
	D[0][7] = -0.125 * onePluseta * onePluszeta;
	D[1][0] = -0.125 * oneMinusxi * oneMinuszeta;
	D[1][1] = -0.125 * onePlusxi * oneMinuszeta;
	D[1][2] = 0.125 * onePlusxi * oneMinuszeta;
	D[1][3] = 0.125 * oneMinusxi * oneMinuszeta;
	D[1][4] = -0.125 * oneMinusxi * onePluszeta;
	D[1][5] = -0.125 * onePlusxi * onePluszeta;
	D[1][6] = 0.125 * onePlusxi * onePluszeta;
	D[1][7] = 0.125 * oneMinusxi * onePluszeta;
	D[2][0] = -0.125 * oneMinusxi * oneMinuseta;
	D[2][1] = -0.125 * onePlusxi * oneMinuseta;
	D[2][2] = -0.125 * onePlusxi * onePluseta;
	D[2][3] = -0.125 * oneMinusxi * onePluseta;
	D[2][4] = 0.125 * oneMinusxi * oneMinuseta;
	D[2][5] = 0.125 * onePlusxi * oneMinuseta;
	D[2][6] = 0.125 * onePlusxi * onePluseta;
	D[2][7] = 0.125 * oneMinusxi * onePluseta;
    
	// calculate Jacobian
	double J[3][3];
    // first column
	J[0][0] =   nd1Crds(0) * D[0][0] + nd2Crds(0) * D[0][1] +
				nd3Crds(0) * D[0][2] + nd4Crds(0) * D[0][3] +
				nd5Crds(0) * D[0][4] + nd6Crds(0) * D[0][5] +
				nd7Crds(0) * D[0][6] + nd8Crds(0) * D[0][7];

	J[1][0] =   nd1Crds(0) * D[1][0] + nd2Crds(0) * D[1][1] +
				nd3Crds(0) * D[1][2] + nd4Crds(0) * D[1][3] +
				nd5Crds(0) * D[1][4] + nd6Crds(0) * D[1][5] +
				nd7Crds(0) * D[1][6] + nd8Crds(0) * D[1][7];

	J[2][0] =   nd1Crds(0) * D[2][0] + nd2Crds(0) * D[2][1] +
				nd3Crds(0) * D[2][2] + nd4Crds(0) * D[2][3] +
				nd5Crds(0) * D[2][4] + nd6Crds(0) * D[2][5] +
				nd7Crds(0) * D[2][6] + nd8Crds(0) * D[2][7];
    // second column
	J[0][1] =   nd1Crds(1) * D[0][0] + nd2Crds(1) * D[0][1] +
				nd3Crds(1) * D[0][2] + nd4Crds(1) * D[0][3] +
				nd5Crds(1) * D[0][4] + nd6Crds(1) * D[0][5] +
				nd7Crds(1) * D[0][6] + nd8Crds(1) * D[0][7];

	J[1][1] =   nd1Crds(1) * D[1][0] + nd2Crds(1) * D[1][1] +
				nd3Crds(1) * D[1][2] + nd4Crds(1) * D[1][3] +
				nd5Crds(1) * D[1][4] + nd6Crds(1) * D[1][5] +
				nd7Crds(1) * D[1][6] + nd8Crds(1) * D[1][7];

	J[2][1] =   nd1Crds(1) * D[2][0] + nd2Crds(1) * D[2][1] +
				nd3Crds(1) * D[2][2] + nd4Crds(1) * D[2][3] +
				nd5Crds(1) * D[2][4] + nd6Crds(1) * D[2][5] +
				nd7Crds(1) * D[2][6] + nd8Crds(1) * D[2][7];
    // third column
	J[0][2] =   nd1Crds(2) * D[0][0] + nd2Crds(2) * D[0][1] +
				nd3Crds(2) * D[0][2] + nd4Crds(2) * D[0][3] +
				nd5Crds(2) * D[0][4] + nd6Crds(2) * D[0][5] +
				nd7Crds(2) * D[0][6] + nd8Crds(2) * D[0][7];

	J[1][2] =   nd1Crds(2) * D[1][0] + nd2Crds(2) * D[1][1] +
				nd3Crds(2) * D[1][2] + nd4Crds(2) * D[1][3] +
				nd5Crds(2) * D[1][4] + nd6Crds(2) * D[1][5] +
				nd7Crds(2) * D[1][6] + nd8Crds(2) * D[1][7];

	J[2][2] =   nd1Crds(2) * D[2][0] + nd2Crds(2) * D[2][1] +
				nd3Crds(2) * D[2][2] + nd4Crds(2) * D[2][3] +
				nd5Crds(2) * D[2][4] + nd6Crds(2) * D[2][5] +
				nd7Crds(2) * D[2][6] + nd8Crds(2) * D[2][7];

	double detJ = J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + 
		          J[0][2] * J[1][0] * J[2][1] - J[0][0] * J[1][2] * J[2][1] - 
				  J[0][1] * J[1][0] * J[2][2] - J[0][2] * J[1][1] * J[2][0];
	
	if (detJ < 0) {
		opserr << "BrickEight::shapeFunction(xi,eta) -- detJ < 0 ";
		exit(-1);
		}

    // calculate the inverse of Jacobian
	double oneOverdetJ = 1.0 / detJ;
	double L[3][3];

	L[0][0] = oneOverdetJ * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
	L[1][0] = -oneOverdetJ * (J[1][0] * J[2][2] - J[1][2] * J[2][0]);
	L[2][0] = oneOverdetJ * (J[1][0] * J[2][1] - J[2][0] * J[1][1]);

	L[0][1] = -oneOverdetJ * (J[0][1] * J[2][2] - J[0][2] * J[2][1]);
	L[1][1] = oneOverdetJ * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
	L[2][1] = -oneOverdetJ * (J[0][0] * J[2][1] - J[2][0] * J[0][1]);

	L[0][2] = oneOverdetJ * (J[0][1] * J[1][2] - J[0][2] * J[1][1]);
	L[1][2] = -oneOverdetJ * (J[0][0] * J[1][2] - J[0][2] * J[1][0]);
	L[2][2] = oneOverdetJ * (J[0][0] * J[1][1] - J[1][0] * J[0][1]);

	// calculate B matrix
    shp[0][0] = L[0][0] * D[0][0] + L[0][1] * D[1][0] + L[0][2] * D[2][0];	// N_1,1
    shp[0][1] = L[0][0] * D[0][1] + L[0][1] * D[1][1] + L[0][2] * D[2][1];	// N_2,1
    shp[0][2] = L[0][0] * D[0][2] + L[0][1] * D[1][2] + L[0][2] * D[2][2];	// N_3,1
    shp[0][3] = L[0][0] * D[0][3] + L[0][1] * D[1][3] + L[0][2] * D[2][3];	// N_4,1
	shp[0][4] = L[0][0] * D[0][4] + L[0][1] * D[1][4] + L[0][2] * D[2][4];	// N_5,1
    shp[0][5] = L[0][0] * D[0][5] + L[0][1] * D[1][5] + L[0][2] * D[2][5];	// N_6,1
    shp[0][6] = L[0][0] * D[0][6] + L[0][1] * D[1][6] + L[0][2] * D[2][6];	// N_7,1
    shp[0][7] = L[0][0] * D[0][7] + L[0][1] * D[1][7] + L[0][2] * D[2][7];	// N_8,1
	
    shp[1][0] = L[1][0] * D[0][0] + L[1][1] * D[1][0] + L[1][2] * D[2][0];	// N_1,2
    shp[1][1] = L[1][0] * D[0][1] + L[1][1] * D[1][1] + L[1][2] * D[2][1];	// N_2,2
    shp[1][2] = L[1][0] * D[0][2] + L[1][1] * D[1][2] + L[1][2] * D[2][2];	// N_3,2
    shp[1][3] = L[1][0] * D[0][3] + L[1][1] * D[1][3] + L[1][2] * D[2][3];	// N_4,2
	shp[1][4] = L[1][0] * D[0][4] + L[1][1] * D[1][4] + L[1][2] * D[2][4];	// N_5,2
    shp[1][5] = L[1][0] * D[0][5] + L[1][1] * D[1][5] + L[1][2] * D[2][5];	// N_6,2
    shp[1][6] = L[1][0] * D[0][6] + L[1][1] * D[1][6] + L[1][2] * D[2][6];	// N_7,2
    shp[1][7] = L[1][0] * D[0][7] + L[1][1] * D[1][7] + L[1][2] * D[2][7];	// N_8,2

	shp[2][0] = L[2][0] * D[0][0] + L[2][1] * D[1][0] + L[2][2] * D[2][0];	// N_1,2
    shp[2][1] = L[2][0] * D[0][1] + L[2][1] * D[1][1] + L[2][2] * D[2][1];	// N_2,2
    shp[2][2] = L[2][0] * D[0][2] + L[2][1] * D[1][2] + L[2][2] * D[2][2];	// N_3,2
    shp[2][3] = L[2][0] * D[0][3] + L[2][1] * D[1][3] + L[2][2] * D[2][3];	// N_4,2
	shp[2][4] = L[2][0] * D[0][4] + L[2][1] * D[1][4] + L[2][2] * D[2][4];	// N_5,2
    shp[2][5] = L[2][0] * D[0][5] + L[2][1] * D[1][5] + L[2][2] * D[2][5];	// N_6,2
    shp[2][6] = L[2][0] * D[0][6] + L[2][1] * D[1][6] + L[2][2] * D[2][6];	// N_7,2
    shp[2][7] = L[2][0] * D[0][7] + L[2][1] * D[1][7] + L[2][2] * D[2][7];	// N_8,2

    return detJ;
}


double 
BrickEight::shapeFunction(double xi, double eta, int faceTag)
{
	if (faceTag < 1 || faceTag > 8) {
		opserr << "BrickEight::shapeFunction -- improper faceTag: " << faceTag << "for BrickEight\n";
		exit(-1);
		}

	int fNum = faceTag-1;

	int node1 = npface[fNum][0];
	int node2 = npface[fNum][1];
	int node3 = npface[fNum][2];
	int node4 = npface[fNum][3];

	const Vector& nd1Crds = theNodes[node1]->getCrds();
	const Vector& nd2Crds = theNodes[node2]->getCrds();
	const Vector& nd3Crds = theNodes[node3]->getCrds();
	const Vector& nd4Crds = theNodes[node4]->getCrds();

	double oneMinuseta = 1.0 - eta;
	double onePluseta = 1.0 + eta;
	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp2[0] = 0.25 * oneMinusxi * oneMinuseta;	// N_1
	shp2[1] = 0.25 * onePlusxi * oneMinuseta;	// N_2
	shp2[2] = 0.25 * onePlusxi * onePluseta;	// N_3
	shp2[3] = 0.25 * oneMinusxi * onePluseta;	// N_4

	// Derivative matrix in reference coordinate system
	double D[4][2];
	//dNdxi
	D[0][0] = -0.25 * oneMinuseta;
	D[1][0] =  0.25 * oneMinuseta;
	D[2][0] =  0.25 * onePluseta;
	D[3][0] = -0.25 * onePluseta;
    //dNdeta
	D[0][1] = -0.25 * oneMinusxi;
	D[1][1] = -0.25 * onePlusxi;
	D[2][1] =  0.25 * onePlusxi;
	D[3][1] =  0.25 * oneMinusxi;

	double E1[3], E2[3];
	//e1
	E1[0] = nd1Crds(0) * D[0][0] + nd2Crds(0) * D[1][0] + nd3Crds(0) * D[2][0] + nd4Crds(0) * D[3][0];
	E1[1] = nd1Crds(1) * D[0][0] + nd2Crds(1) * D[1][0] + nd3Crds(1) * D[2][0] + nd4Crds(1) * D[3][0];
	E1[2] = nd1Crds(2) * D[0][0] + nd2Crds(2) * D[1][0] + nd3Crds(2) * D[2][0] + nd4Crds(2) * D[3][0];
    //e2
	E2[0] = nd1Crds(0) * D[0][1] + nd2Crds(0) * D[1][1] + nd3Crds(0) * D[2][1] + nd4Crds(0) * D[3][1];
	E2[1] = nd1Crds(1) * D[0][1] + nd2Crds(1) * D[1][1] + nd3Crds(1) * D[2][1] + nd4Crds(1) * D[3][1];
	E2[2] = nd1Crds(2) * D[0][1] + nd2Crds(2) * D[1][1] + nd3Crds(2) * D[2][1] + nd4Crds(2) * D[3][1];
	double x1 = E1[1] * E2[2] - E1[2] * E2[1];
	double x2 = E1[2] * E2[0] - E1[0] * E2[2];
	double x3 = E1[0] * E2[1] - E1[1] * E2[0];
	double X = x1 * x1 + x2 * x2 + x3 * x3;
	double JGamma = sqrt(X);

	return JGamma;
}


double 
BrickEight::getInterpolatedTemp(int tag)
{
	int fNum = tag -1;
    int node1 = npface[fNum][0];
	int node2 = npface[fNum][1];
	int node3 = npface[fNum][2];
	int node4 = npface[fNum][3];

	const Vector& Temp1 = theNodes[node1]->getTrialTemperature();
	const Vector& Temp2 = theNodes[node2]->getTrialTemperature();
	const Vector& Temp3 = theNodes[node3]->getTrialTemperature();
	const Vector& Temp4 = theNodes[node4]->getTrialTemperature();

	double T_hat;
	T_hat = shp2[0] * Temp1(0) + shp2[1] * Temp2(0) + shp2[2] * Temp3(0)
		    + shp2[3] * Temp4(0);

	return T_hat;
}


