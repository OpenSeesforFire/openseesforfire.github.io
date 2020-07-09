/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.6 $
// $Date: 2007-07-27 19:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint3dThermal.cpp,v $
                                                                        
// Written: NM (nmitra@u.washington.edu)
// Created: April 2002
// Last Revised: January 2007
//
// Description: This file contains the class implementation for beam-column joint.
// This element is a 4 noded 12 dof (3 dof at each node) finite area super-element.
// The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint.
// Details about this element can be found in 
// Mitra and Lowes, J. Struct. Eng. ASCE, Volume 133, Number 1 (January 2007), pp. 105-120
// Updates: Several concerning Joint formulation (presently a revised formulation for joints)


#include <BeamColumnJoint3dThermal.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>
#include <ElementalLoad.h>

#include <UniaxialMaterial.h>
#include <string.h>
#include <math.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <string>

//static Matrix CoupledZeroLengthM6(6, 6);   // class wide matrix for 6*6
static Matrix TwoNodeM6(6,6);   // class wide matrix for 6*6
static Vector TwoNodeV6(6);   // class wide Vector for size 6
int  numMaterials1d = 3;

void* OPS_BeamColumnJoint3dThermal()
{
    if (OPS_GetNumRemainingInputArgs() < 6) {
	opserr << "WARNING insufficient arguments\n";
	  opserr << "Want: element beamColumnJoint eleTag? node1? node2? \n";
	  opserr << "matTag1? matTag2? matTag3?\n";
	  //opserr << "<ElementHeightFactor? ElementWidthFactor?>\n";
	  return 0;
    }

    int idata[6];
    int numdata = 6;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    double data[2] = {1.0, 1.0};
    numdata = 2;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr<<"WARNING: invalid double inputs\n";
	    return 0;
	}
    }

    UniaxialMaterial* mats[3];
    for (int i = 0; i < 3; i++) {
	mats[i] = OPS_getUniaxialMaterial(idata[3+i]);
	if (mats[i] == 0) {
	    opserr<<"WARNING: material "<<idata[3+i]<<" is not defined\n";
	    return 0;
	}
    }

    return new BeamColumnJoint3dThermal(idata[0],idata[1],idata[2],
				 *mats[0],*mats[1],*mats[2],data[0],data[1]);

}

// full constructors:
BeamColumnJoint3dThermal::BeamColumnJoint3dThermal(int tag,int Nd1, int Nd2,
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3):
  Element(tag,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(2),
 theMatrix(0), d0(0), v0(0), theVector(0), theLoad(0), elemType(1),
	dimension(2), nDOF(0), transformation(3, 3)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 2)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJointThermal " << tag << "failed to create an ID of size 2" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;


	MaterialPtr = new UniaxialMaterial*[3];

	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }



	nodePtr[0] = 0;
	nodePtr[1] = 0;

// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
  MaterialPtr[1] = theMat2.getCopy();
  if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2"<< endln;}
  MaterialPtr[2] = theMat3.getCopy();
  if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3"<< endln;}
  
  mInitialize = 1;
}

// full constructors:
BeamColumnJoint3dThermal::BeamColumnJoint3dThermal(int tag,int Nd1, int Nd2, 
				     UniaxialMaterial& theMat1,
				     UniaxialMaterial& theMat2,
				     UniaxialMaterial& theMat3,
					 double elHgtFac, double elWdtFac):
  Element(tag,ELE_TAG_BeamColumnJoint3d), connectedExternalNodes(2),
	d0(0), v0(0), theMatrix(0), theVector(0), theLoad(0), elemType(1),
	dimension(2), nDOF(0), transformation(3, 3)
{
	// ensure the connectedExternalNode ID is of correct size & set values
 
	if (connectedExternalNodes.Size() != 2)
      opserr << "ERROR : BeamColumnJoint::BeamColumnJointThermal " << tag << "failed to create an ID of size 2" << endln;

	connectedExternalNodes(0) = Nd1 ;
    connectedExternalNodes(1) = Nd2 ;
    

	MaterialPtr = new UniaxialMaterial*[3];
	deform = Vector(3);

	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }


	// added
	nodePtr[0] = 0;
	nodePtr[1] = 0;


	// get a copy of the material and check we obtained a valid copy
  MaterialPtr[0] = theMat1.getCopy();
  if (!MaterialPtr[0]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 1" << endln;}
  MaterialPtr[1] = theMat2.getCopy();
  if (!MaterialPtr[1]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 2"<< endln;}
  MaterialPtr[2] = theMat3.getCopy();
  if (!MaterialPtr[2]){
		opserr << "ERROR : BeamColumnJoint::Constructor failed to get a copy of material 3"<< endln;}
  mInitialize = 1;
}

// default constructor:
BeamColumnJoint3dThermal::BeamColumnJoint3dThermal():
  Element(0,ELE_TAG_BeamColumnJoint2d), connectedExternalNodes(2),
 d0(0), v0(0), theMatrix(0), theVector(0), theLoad(0),elemType(1)
{
	nodePtr[0] = 0;
	nodePtr[1] = 0;
	for (int x = 0; x <3; x++)
	{	MaterialPtr[x] = 0; }
	mInitialize = 0;
   // does nothing (invoked by FEM_ObjectBroker)
}

//  destructor:
BeamColumnJoint3dThermal::~BeamColumnJoint3dThermal()
{
	for (int i =0; i<3; i++)
	{
	    if (MaterialPtr[i] != 0)
		delete MaterialPtr[i];
	}

	if (MaterialPtr)
		 delete [] MaterialPtr;
}

// public methods
int
BeamColumnJoint3dThermal::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
BeamColumnJoint3dThermal::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **BeamColumnJoint3dThermal::getNodePtrs(void)
{
	return nodePtr;
}

int
BeamColumnJoint3dThermal::getNumDOF(void) 
{
    return 6;
}

void
BeamColumnJoint3dThermal::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	opserr << "ERROR : BeamColumnJoint::setDomain -- Domain is null" << endln;
	
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nodePtr[0] = 0;
	nodePtr[1] = 0;
    }

	//node pointers set
    for (int i = 0; i < 2; i++ ) {
		nodePtr[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR : BeamColumnJoint::setDomain -- node pointer is null"<< endln;
			exit(-1); // donot go any further - otherwise segmentation fault
		}
	}

	// ensure connected nodes have correct dof's
	int dofNd1 = nodePtr[0]->getNumberDOF();
	int dofNd2 = nodePtr[1]->getNumberDOF();


	// if differing dof at the ends - print a warning message
	if (dofNd1 != dofNd2) {
		opserr << "WARNING ZeroLength::setDomain(): nodes " << connectedExternalNodes <<
			"have differing dof at ends for ZeroLength " << this->getTag() << endln;
		return;
	}


	// Check that length is zero within tolerance
	const Vector& end1Crd = nodePtr[0]->getCrds();
	const Vector& end2Crd = nodePtr[1]->getCrds();
	Vector diff = end1Crd - end2Crd;
	double L = diff.Norm();
	double v1 = end1Crd.Norm();
	double v2 = end2Crd.Norm();
	double vm;

	vm = (v1 < v2) ? v2 : v1;

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	Vector Node1(end1Crd);
	Vector Node2(end2Crd);

	theMatrix = &TwoNodeM6;
	theVector = &TwoNodeV6;

	// get trial displacements and take difference
	const Vector& disp1 = nodePtr[0]->getTrialDisp();
	const Vector& disp2 = nodePtr[1]->getTrialDisp();
	Vector  diffD = disp2 - disp1;
	const Vector& vel1 = nodePtr[0]->getTrialVel();
	const Vector& vel2 = nodePtr[1]->getTrialVel();
	Vector  diffV = vel2 - vel1;

	// to avoid incorrect results, do not set initial disp/vel upon call of null constructor
	// when using database commands
	if (mInitialize == 1) {
		if (diffD != 0.0)
			d0 = new Vector(diffD);

		if (diffV != 0)
			v0 = new Vector(diffV);
	}

}   

int
BeamColumnJoint3dThermal::commitState(void)
{

	int errCode = 0;

	// commit material models
	for (int i = 0; i < 3; i++)
		errCode += MaterialPtr[i]->commitState();

	// commit the base class
	errCode += this->Element::commitState();

	return errCode;
}

int
BeamColumnJoint3dThermal::revertToLastCommit(void)
{
	int errCode = 0;
	// revert material models
	for (int i = 0; i < 3; i++)
		errCode += MaterialPtr[i]->revertToLastCommit();

	return errCode;
}

int
BeamColumnJoint3dThermal::revertToStart(void)
{
	int mcs = 0;
	for (int j=0; j<3; j++)
	{
		if (MaterialPtr[j] != 0) mcs = MaterialPtr[j]->revertToStart();
		if (mcs != 0) break;
	}
	
	return mcs;
}

int
BeamColumnJoint3dThermal::update(void)
{	

	double strain;
	double strainRate;

	// get trial displacements and take difference
	const Vector& disp1 = nodePtr[0]->getTrialDisp();
	const Vector& disp2 = nodePtr[1]->getTrialDisp();
	Vector  diff = disp2 - disp1;
	const Vector& vel1 = nodePtr[0]->getTrialVel();
	const Vector& vel2 = nodePtr[1]->getTrialVel();
	Vector  diffv = vel2 - vel1;

	if (d0 != 0)
		diff -= *d0;

	if (v0 != 0)
		diffv -= *v0;

	// loop over 1d materials

	//    Matrix& tran = *t1d;
	int ret = 0;
	for (int mat = 0; mat < numMaterials1d; mat++) {
		// compute strain and rate; set as current trial for material
		strain = this->computeCurrentStrain1d(mat, diff);
		deform(mat) = strain;
		//opserr << "  strain " << mat << " " << strain;
		//strainRate = this->computeCurrentStrain1d(mat, diffv);
		ret += MaterialPtr[mat]->setTrialStrain(strain);
	} 
	if (abs(deform(2)) > 0.2) { 
		Domain* theDomain = OPS_GetDomain();
		opserr << "WARNing::Connection Failure! Element removed";
		theDomain->removeElement(this->getTag());
	}

	return ret;
}

const Matrix &
BeamColumnJoint3dThermal::getTangentStiff(void)
{
	double E;


	// stiff is a reference to the matrix holding the stiffness matrix
	Matrix& stiff = *theMatrix;

	// zero stiffness matrix
	stiff.Zero();

	// loop over 1d materials

	//Matrix& tran = *t1d;;
	for (int mat = 0; mat < numMaterials1d; mat++) {

		// get tangent for material
		E = MaterialPtr[mat]->getTangent();

		// compute contribution of material to tangent matrix
		if (mat == 0) {
			stiff(mat, mat) = E;
			stiff(mat, mat+3) = -E;
			stiff(mat+3, mat) = -E;
			stiff(mat + 3, mat + 3) = E;
		}
		else if (mat == 1) {
			stiff(mat, mat) = E;
			stiff(mat, mat + 3) = -E;
			stiff(mat + 3, mat) = -E;
			stiff(mat + 3, mat + 3) = E;
		}
		else if (mat == 2) {
			stiff(mat, mat) = E;
			stiff(mat, mat + 3) = -E;
			stiff(mat + 3, mat) = -E;
			stiff(mat + 3, mat + 3) = E;
		}
		

	}

	// end loop over 1d materials 

	//complete symmetric stiffness matrix
	//for (int i = 0; i < 6; i++)
	//	for (int j = 0; j < i; j++)
	//		stiff(j, i) = stiff(i, j);
	//opserr << "stiff " << stiff << endln;
	return stiff;

}

const Matrix &
BeamColumnJoint3dThermal::getInitialStiff(void)
{
	return getTangentStiff();
}

const Vector &
BeamColumnJoint3dThermal::getResistingForce(void)
{
	//Global forces at two nodes (Fx, Fy, M)
	double force;
	// zero the residual
	theVector->Zero();

	// loop over 1d materials
	for (int mat = 0; mat < numMaterials1d; mat++) {

		// get resisting force for material
		force = MaterialPtr[mat]->getStress();

		// compute residual due to resisting force
		if (mat == 0) {
			(*theVector)(mat) = -force;
			(*theVector)(mat + 3) = force;
		}
		else if (mat == 1) {
			(*theVector)(mat) = -force;
			(*theVector)(mat + 3) = force;
		}
		else if (mat == 2) {
			(*theVector)(mat) = -force;
			(*theVector)(mat + 3) = force;
		}

	} // end loop over 1d materials 

	//opserr << "   Resisting F: " << *theVector << endln;
	return *theVector;
}


    
const Matrix &
BeamColumnJoint3dThermal::getDamp(void)
{
	//not applicable (stiffness being returned)
	return *theMatrix;
}

const Matrix &
BeamColumnJoint3dThermal::getMass(void)
{ 
	//not applicable  (stiffness being returned)
	return *theMatrix;
}

void 
BeamColumnJoint3dThermal::zeroLoad(void)
{
	//not applicable  
	return;
}

int 
BeamColumnJoint3dThermal::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	//aDD THERMAL ACTION
	int type;
	const Vector& data = theLoad->getData(type, loadFactor);
	//opserr << "Tmp Data " << data << endln;
	if (type == LOAD_TAG_Beam2dThermalAction) {
		// load not inside fire load pattern
		//const Vector &data = theLoad->getData(type, loadFactor);

		/// This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[END]
		double tangent = 0.0;
		double ThermalElongation = 0.0;
		double FiberTemperature = 0;
		double FiberTempMax = 0;
		/*
		Vector dataMixV(27);
			for (int m = 0; m < 9; m++) {
				dataMixV(2 * m) = data(2 * m); //Linear temperature interpolation
				dataMixV(2 * m + 1) = data(2 * m + 1);
				dataMixV(18 + m) = 1000;
			}
		*/
			
			FiberTemperature = data(0);
			static Vector tData(4);
			static Information iData(tData);
			tData(0) = FiberTemperature;
			tData(1) = tangent;
			tData(2) = ThermalElongation;
			tData(3) = FiberTempMax;
			iData.setVector(tData);
			

			for (int mat = 0; mat < numMaterials1d; mat++) {
				MaterialPtr[mat]->getVariable("ElongTangent", iData);
				tData = iData.getData();
			}
			//const Vector& s = MaterialPtr[i]->data(1);    //contribuited by ThermalElongation
#ifdef _BDEBUG
			if (this->getTag() == 1)
				opserr << "Thermal Stress " << s << endln;
#endif
		


	}


	return 0;
}

int 
BeamColumnJoint3dThermal::addInertiaLoadToUnbalance(const Vector &accel)
{
	//not applicable
	return 0;
}


const Vector &
BeamColumnJoint3dThermal::getResistingForceIncInertia()
{	
  //not applicable (residual being returned)
	return *theVector;
}

int
BeamColumnJoint3dThermal::sendSelf(int commitTag, Channel &theChannel)
{
	// yet to do.
	return -1;
}

int
BeamColumnJoint3dThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// yet to do.
	return -1;
}

int
BeamColumnJoint3dThermal::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	const Vector &node1Crd = nodePtr[0]->getCrds();
	const Vector &node2Crd = nodePtr[1]->getCrds();	


	const Vector &node1Disp = nodePtr[0]->getDisp();
	const Vector &node2Disp = nodePtr[1]->getDisp();    
 

	static Vector v1(3);
	static Vector v2(3);

	
	// calculate the current coordinates of four external nodes
	for (int i=0; i<2; i++) 
    {
		v1(i) = node1Crd(i)+node1Disp(i)*fact;
		v2(i) = node2Crd(i)+node2Disp(i)*fact;

	}
	

	return 0;

}

void
BeamColumnJoint3dThermal::Print(OPS_Stream &s, int flag)
{
	s << "Element: " << this->getTag() << " Type: Beam Column Joint Thermal" << endln;
	for (int i = 0; i<2; i++)
	{
		s << "Node :" << connectedExternalNodes(i);
		s << "DOF :" << nodePtr[i]->getNumberDOF();
	}
	return;
}

Response*
BeamColumnJoint3dThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  // we will compare argv[0] to determine the type of response required
  
 if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"Deformation") == 0 || strcmp(argv[0], "deform") == 0 || strcmp(argv[0], "Deform") == 0)
    return new ElementResponse(this,1,Vector(3));
 else if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "Force") == 0 || strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "Forces") == 0)
	 return new ElementResponse(this, 2, Vector(6));
  
  else
    return 0;
}

int 
BeamColumnJoint3dThermal::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1:  // deform
		return eleInfo.setVector(deform);
	case 2:  // global forces
		return eleInfo.setVector(this->getResistingForce());

	default:
		return -1;
	}
}

int
BeamColumnJoint3dThermal::setParameter (char **argv, int argc, Information &info)
{
  return -1;
}
    
int
BeamColumnJoint3dThermal::updateParameter (int parameterID, Information &info)
{
  return -1;
}


// Compute current strain for 1d material mat
// dispDiff are the displacements of node 2 minus those
// of node 1
double
BeamColumnJoint3dThermal::computeCurrentStrain1d(int mat,
	const Vector& dispDiff) const
{
	double strain = 0.0;

	//for (int i = 0; i < 6 / 2; i++) {
		strain = dispDiff(mat) ;
	//}

	return strain;
}



void
BeamColumnJoint3dThermal::updateDir(const Vector& x, const Vector& y)
{
	this->setUp(connectedExternalNodes(0), connectedExternalNodes(1), x, y);
	this->setTran1d(elemType, numMaterials1d);
}


// Set basic deformation-displacement transformation matrix for 1d
// uniaxial materials
void
BeamColumnJoint3dThermal::setTran1d(int elemType,
	int   numMat)
{
	int nDoF;
	enum Dtype { TRANS, ROTATE };

	int   indx, dir;
	Dtype dirType;

	// Create 1d transformation matrix
	if (elemType == 1)
		nDoF = 6;
	else if(elemType == 2)
		nDoF = 6;
	else if (elemType == 3)
		nDoF = 12;

	t1d = new Matrix(numMat, nDoF);

	if (t1d == 0)
		opserr << "FATAL ZeroLength::setTran1d - can't allocate 1d transformation matrix\n";

	// Use reference for convenience and zero matrix.
	Matrix& tran = *t1d;
	tran.Zero();

	// loop over materials, setting row in tran for each material depending on dimensionality of element

	for (int i = 0; i < numMat; i++) {

		dir = (*dir1d)(i);	// direction 0 to 5;
		indx = dir % 3;		// direction 0, 1, 2 for axis of translation or rotation

		// set direction type to translation or rotation
		dirType = (dir < 3) ? TRANS : ROTATE;

		// now switch on dimensionality of element

		switch (elemType) {
//D2N6
		case 1:
			if (dirType == TRANS) {
				tran(i, 3) = transformation(indx, 0);
				tran(i, 4) = transformation(indx, 1);
				tran(i, 5) = 0.0;
			}
			else if (dirType == ROTATE) {
				tran(i, 3) = 0.0;
				tran(i, 4) = 0.0;
				tran(i, 5) = transformation(indx, 2);
			}
			break;
//D3N6
		case 2:
			if (dirType == TRANS) {
				tran(i, 3) = transformation(indx, 0);
				tran(i, 4) = transformation(indx, 1);
				tran(i, 5) = transformation(indx, 2);
			}
			break;
//D3N12
		case 3:
			if (dirType == TRANS) {
				tran(i, 6) = transformation(indx, 0);
				tran(i, 7) = transformation(indx, 1);
				tran(i, 8) = transformation(indx, 2);
				tran(i, 9) = 0.0;
				tran(i, 10) = 0.0;
				tran(i, 11) = 0.0;
			}
			else if (dirType == ROTATE) {
				tran(i, 6) = 0.0;
				tran(i, 7) = 0.0;
				tran(i, 8) = 0.0;
				tran(i, 9) = transformation(indx, 0);
				tran(i, 10) = transformation(indx, 1);
				tran(i, 11) = transformation(indx, 2);
			}
			break;

		} // end switch

		// fill in first half of transformation matrix with
		// negative sign

		for (int j = 0; j < nDOF / 2; j++)
			tran(i, j) = -tran(i, j + nDOF / 2);

	} // end loop over 1d materials
}


void
BeamColumnJoint3dThermal::setUp(int Nd1, int Nd2,
	const Vector& x,
	const Vector& yp)
{
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2)
		opserr << "FATAL ZeroLength::setUp - failed to create an ID of correct size\n";

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	int i;
	for (i = 0; i < 2; i++)
		nodePtr[i] = 0;

	// check that vectors for orientation are correct size
	if (x.Size() != 3 || yp.Size() != 3)
		opserr << "FATAL ZeroLength::setUp - incorrect dimension of orientation vectors\n";

	// establish orientation of element for the transformation matrix
	// z = x cross yp
	Vector z(3);
	z(0) = x(1) * yp(2) - x(2) * yp(1);
	z(1) = x(2) * yp(0) - x(0) * yp(2);
	z(2) = x(0) * yp(1) - x(1) * yp(0);

	// y = z cross x
	Vector y(3);
	y(0) = z(1) * x(2) - z(2) * x(1);
	y(1) = z(2) * x(0) - z(0) * x(2);
	y(2) = z(0) * x(1) - z(1) * x(0);

	// compute length(norm) of vectors
	double xn = x.Norm();
	double yn = y.Norm();
	double zn = z.Norm();

	// check valid x and y vectors, i.e. not parallel and of zero length
	if (xn == 0 || yn == 0 || zn == 0) {
		opserr << "FATAL ZeroLength::setUp - invalid vectors to constructor\n";
	}

	// create transformation matrix of direction cosines
	for (i = 0; i < 3; i++) {
		transformation(0, i) = x(i) / xn;
		transformation(1, i) = y(i) / yn;
		transformation(2, i) = z(i) / zn;
	}

}