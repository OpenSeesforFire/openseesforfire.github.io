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
                                                                        
// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/IGAQuad/IGAQuad.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for IGAQuad.

//detect memory leak
//#include "stdafx.h"
//#include "vld.h"

//int _tmain(int argc, _TCHAR* argv[])
//{
//    char* pBuf = new char[200];
//
//    return 0;
//}


#include <IGAQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>



//Definiton of single element from the script input
void* OPS_IGAQuad()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 2) {
	opserr << "WARNING -- IGAquad element expects ndm=2 and ndf=2\n";
	return 0;
    }
    
 /*   if (OPS_GetNumRemainingInputArgs() < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element IGAQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
	return 0;
    }*/


    //int num = OPS_GetNumRemainingInputArgs()-3;
    // IGAQuadId, iNode, jNode, kNode, lNode: Now they are IGA control points
    int num = 1;
    int ex = 1, nex = 1, ey = 1, ney = 1;
    OPS_GetIntInput(&num, &ex);
    OPS_GetIntInput(&num, &nex);
    OPS_GetIntInput(&num, &ey);
    OPS_GetIntInput(&num, &ney);
    num = 2;
    int* obfs1 = new int[2];
    OPS_GetIntInput(&num, obfs1);

    num = (obfs1[0] + 1) * (obfs1[1] + 1);
    int* CPIds = new int[num];
    OPS_GetIntInput(&num, CPIds);
    
    num = 1;
    int numKnotVect_x=0;
    OPS_GetIntInput(&num, &numKnotVect_x);
    num = numKnotVect_x;
    double* datax = new double[num];
    OPS_GetDoubleInput(&num, datax);
    Vector KnotVect_x1(numKnotVect_x);
    for (int i = 1; i <= numKnotVect_x; i++) { KnotVect_x1(i - 1) = datax[i - 1]; }
    delete[] datax;
    
    int numKnotVect_y = 0;
    OPS_GetIntInput(&num, &numKnotVect_y);
    num = numKnotVect_y;
    double* datay = new double[num];
    OPS_GetDoubleInput(&num, datay);
    Vector KnotVect_y1(numKnotVect_y);
    for (int i = 1; i <= numKnotVect_y; i++) { KnotVect_y1(i - 1) = datay[i - 1]; }
    delete[] datay;
    
    num = 2;
    int numMults[2];
    OPS_GetIntInput(&num, numMults);


    //int num = 5;

    double thk = 1.0;
    num = 1;
    if (OPS_GetDoubleInput(&num,&thk) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    const char* type = OPS_GetString();

    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << matTag;
	//opserr << "\nIGAQuad element: " << idata[0] << endln;
	return 0;
    }

    // p, rho, b1, b2
    double data[4] = {0,0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 4) {
	num = 4;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,data) < 0) {
	    opserr<<"WARNING: invalid integer data\n";
	    return 0;
	}	
    }
    int tag = ex + ney * (ey - 1);
    int numCPs1 = (obfs1[0] + 1) * (obfs1[1] + 1);
    int ndof1 = numCPs1 * 2;
    return new IGAQuad(tag,numCPs1,ex,nex,ey,ney,obfs1,ndof1,CPIds,KnotVect_x1,KnotVect_y1,numMults,
			                *mat,type,thk,data[0],data[1],data[2],data[3]);
	// elementTag, cp1, cp2, cp3, cp4. material pointer, type, thickness, load, pressure, density, Ki
}


double IGAQuad::matrixData[324];
Matrix IGAQuad::K(18,18);
Vector IGAQuad::P(18);
double IGAQuad::shp[3][9];
double IGAQuad::pts[9][2];
double IGAQuad::wts[9];

//Matrix* IGAQuad::N0n;

IGAQuad::IGAQuad(int tag, int numCPs1, int ex, int nex, int ey, int ney, int* obfs1, int ndof1,
    int* CPIds, Vector KnotVect_x1, Vector KnotVect_y1, int* numMults,
    NDMaterial& m, const char* type, double t,
    double p, double r, double b1, double b2)
    :Element(tag, ELE_TAG_IGAQuad),
    theMaterial(0), connectedExternalNodes(4),
    Q(8), pressureLoad(8), thickness(t), applyLoad(0), pressure(p), rho(r), Ki(0)
{
    //int NGPss[2];
    //for (int i = 1; i <= 2; i++) { NGPss[i - 1] = obfs[i - 1] + 1; }
    numCPs = numCPs1;
    //const int nCPs = numCPs;
    ndof = ndof1;
    //const int numDof = ndof;
    obfs[0] = obfs1[0];
    obfs[1] = obfs1[1];
    delete[] obfs1;
    KnotVect_y = KnotVect_y1;
    KnotVect_x = KnotVect_x1;
    //const int numMatrixData = ndof * ndof;
    //static double matrixData[324];
    //static Matrix K(numDof, numDof);
    //static Vector P(numDof);
    //static double shp[3][9];
    //static double pts[9][2];
    //static double wts[9];

    eleIdInfo.resize(4);
    eleIdInfo(0) = ex;
    eleIdInfo(1) = nex;
    eleIdInfo(2) = ey;
    eleIdInfo(3) = ney;

    N0nx = new Matrix[obfs[0] + 1];
    N0ny = new Matrix[obfs[1] + 1];


	if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
	    && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
	  opserr << "IGAQuad::IGAQuad -- improper material type: " << type << "for IGAQuad\n";
	  exit(-1);
	}

	// Body forces
	b[0] = b1;
	b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[numCPs];
    //theNodes = new Node *[numCPs];
    
    if (theMaterial == 0) {
      opserr << "IGAQuad::IGAQuad - failed allocate material model pointer\n";
      exit(-1);
    }

    for (int i = 0; i < numCPs; i++) {
      
      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);
			
      // Check allocation
      if (theMaterial[i] == 0) {
          opserr << "IGAQuad::IGAQuad -- failed to get a copy of material model\n";
          exit(-1);
      }
    }

    // Set connected external node IDs
    for (int i = 1; i <= numCPs; i++) {
        connectedExternalNodes(i - 1) = CPIds[i - 1];
        theNodes[i - 1] = 0;
    }
    delete[] CPIds;
    //connectedExternalNodes(0) = nd1;
    //connectedExternalNodes(1) = nd2;
    //connectedExternalNodes(2) = nd3;
    //connectedExternalNodes(3) = nd4;
    
    //for (int i=1; i<=numCPs; i++)
    //  theNodes[i-1] = 0;
}


IGAQuad::IGAQuad()
:Element (0,ELE_TAG_IGAQuad),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), pressureLoad(8), thickness(0.0), applyLoad(0), pressure(0.0), Ki(0)
{
  //pts(0,0) = -0.577350269189626;
  //pts[0][1] = -0.577350269189626;
  //pts[1][0] =  0.577350269189626;
  //pts[1][1] = -0.577350269189626;
  //pts[2][0] =  0.577350269189626;
  //pts[2][1] =  0.577350269189626;
  //pts[3][0] = -0.577350269189626;
  //pts[3][1] =  0.577350269189626;
  //
  //wts[0] = 1.0;
  //wts[1] = 1.0;
  //wts[2] = 1.0;
  //wts[3] = 1.0;

    for (int i=0; i<4; i++)
      theNodes[i] = 0;
}

IGAQuad::~IGAQuad()
{    
  for (int i = 0; i < numCPs; i++) {
      if (theMaterial[i]) {
          delete theMaterial[i];
      }
  }
  delete[] N0nx;
  delete[] N0ny;
  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial) {
      delete[] theMaterial;
  }

  if (Ki != 0) {
      delete Ki;
  }
}

int
IGAQuad::getNumExternalNodes() const
{
    return numCPs;
}

const ID&
IGAQuad::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
IGAQuad::getNodePtrs(void) 
{
  return theNodes;
}

int
IGAQuad::getNumDOF()
{
    return ndof;
}

void
IGAQuad::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        for (int i = 1; i <= numCPs; i++) {
            theNodes[i-1] = 0;
        }
        return;
    }
    Vector Nd(numCPs);
    for (int i = 1; i <= numCPs; i++) {
        Nd(i - 1) = connectedExternalNodes(i - 1);
        theNodes[i - 1] = theDomain->getNode(Nd(i - 1));
    }

 //   if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
	////opserr << "FATAL ERROR IGAQuad (tag: %d), node not found in domain",
	////	this->getTag());
	//
	//return;
 //   }
    //Vector dofNd(numCPs);
    //for (int i = 1; i <= numCPs; i++) {
    //    dofNd(i - 1) = theNodes[i - 1]->getNumberDOF();
    //}
    //int dofNd1 = theNodes[0]->getNumberDOF();
    //int dofNd2 = theNodes[1]->getNumberDOF();
    //int dofNd3 = theNodes[2]->getNumberDOF();
    //int dofNd4 = theNodes[3]->getNumberDOF();
    
    //if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
	//opserr << "FATAL ERROR IGAQuad (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	//return;
 //   }
    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
}

int
IGAQuad::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "IGAQuad::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numCPs; i++)
      retVal += theMaterial[i]->commitState();

    return retVal;
}

int
IGAQuad::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numCPs; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
IGAQuad::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numCPs; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}

//update the deformation at intergation points( Gauss). IGA should be considered 
int
IGAQuad::update()
{
    opserr << "Element "<<this->getTag()<<" update is done" << endln;

    Matrix u(2,numCPs);
    //Vector disp_temp = new Vector & [numCPs];
    Vector* disp_temp = new Vector[numCPs];
    for (int i = 1; i <= numCPs; i++) {
        disp_temp[i-1] = theNodes[i - 1]->getTrialDisp();
        for (int j = 1; j <= 2; j++) {
            u(j - 1, i - 1) = (disp_temp[i-1])(j - 1);
        }
        
    }
    delete[] disp_temp;
    //double uTest[4][2];
	//const Vector &disp1 = theNodes[0]->getTrialDisp();
	//const Vector &disp2 = theNodes[1]->getTrialDisp();
	//const Vector &disp3 = theNodes[2]->getTrialDisp();
	//const Vector &disp4 = theNodes[3]->getTrialDisp();
	//
	//static double u[2][4];

	//u[0][0] = disp1(0);
	//u[1][0] = disp1(1);
	//u[0][1] = disp2(0);
	//u[1][1] = disp2(1);
	//u[0][2] = disp3(0);
	//u[1][2] = disp3(1);
	//u[0][3] = disp4(0);
	//u[1][3] = disp4(1);

	static Vector eps(3);

	int ret = 0;

	//creation of knots may be written here
	// Loop over the integration points FEMID:1,2,3,4, IGAID:1,2,4,3 
	for (int j = 1; j <= obfs[1]+1; j++) { // loop x-directional GPs
        for (int i = 1; i <= obfs[0]+1; i++) {// loop y-directional GPs
            // Determine Jacobian for this integration point
            this->shapeFunction(i, j, eleIdInfo, KnotVect_x, KnotVect_y);
            // Interpolate strains
		    //eps = B*u;
		    //eps.addMatrixVector(0.0, B, u, 1.0);
		    eps.Zero();
		    for (int beta = 0; beta < numCPs; beta++) {
			    eps(0) += shp[0][beta]*u(0,beta);
			    eps(1) += shp[1][beta]*u(1,beta);
			    eps(2) += shp[0][beta]*u(1,beta) + shp[1][beta]*u(0,beta);
		    }
            int CPId = (obfs[1]+1) * (j - 1) + i;// for temperaly use
            ret += theMaterial[CPId - 1]->setTrialStrain(eps);
        }
	}


	return ret;
}

//This is key function for stiffness matrix, IGA should be included here
const Matrix&
IGAQuad::getTangentStiff()
{

	K.Zero();

	double dvol;
	double DB[3][2];
    Vector DataJ(2);//J1,J2 for IGA
    double KK[18][18];
    // Loop over the integration points FEMID:1,2,3,4, IGAID:1,2,4,3 
    for (int j = 1; j <= obfs[0] + 1; j++) { // loop x-directional GPs
        for (int i = 1; i <= obfs[1] + 1; i++) {// loop y-directional GPs
            // Determine Jacobian for this integration point
            DataJ = this->shapeFunction(i, j, eleIdInfo, KnotVect_x, KnotVect_y);
            int CPId = (obfs[0] + 1) * (j - 1) + i;// for temperaly use
            dvol = DataJ(0) * DataJ(1) * thickness * wts[CPId - 1];// J1*J2*t*W

            // Get the material tangent
            const Matrix& D = theMaterial[CPId - 1]->getTangent();

            // Perform numerical integration
            //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
            //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

            double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
            double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
            double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

            for (int alpha = 0, ia = 0; alpha < numCPs; alpha++, ia += 2) {
                for (int beta = 0, ib = 0; beta < numCPs; beta++, ib += 2) {

                    DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
                    DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
                    DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
                    DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
                    DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
                    DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);


                    K(ia, ib) += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
                    K(ia, ib + 1) += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
                    K(ia + 1, ib) += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
                    K(ia + 1, ib + 1) += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];
                    //	      matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
                    //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
                    //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
                    //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
                }
            }
            for (int ii = 1; ii <= ndof;ii++) {
                //KK[i - 1] = new double[8];
                for (int jj = 1; jj <= ndof;jj++) {
                    KK[ii - 1][jj - 1] = K(ii - 1, jj - 1);
                }
            }
        }
    }

	return K;
}

//This will be similar to getTangentStiff()
const Matrix&
IGAQuad::getInitialStiff()
{
    if (Ki != 0)
        return *Ki;

    K.Zero();

    double dvol;
    double DB[3][2];
    Vector DataJ(2);

    // Loop over the integration points FEMID:1,2,3,4, IGAID:1,2,4,3 
    for (int j = 1; j <= obfs[0]+1; j++) { // loop x-directional GPs
        for (int i = 1; i <= obfs[1]+1; i++) {// loop y-directional GPs
            //if (j == 1) { ii = i_temp; }
            //else { ii = 3 - i_temp; }
            int CPId = (obfs[0] + 1) * (j - 1) + i;// for temperaly use
            // Determine Jacobian for this integration point
            DataJ = this->shapeFunction(i, j, eleIdInfo, KnotVect_x, KnotVect_y);
            dvol = DataJ(0) * DataJ(1) * thickness * wts[CPId - 1];// J1*J2*t*W

            // Get the material tangent
            const Matrix& D = theMaterial[CPId - 1]->getInitialTangent();

            double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
            double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
            double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

            // Perform numerical integration
            //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
            //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
            for (int beta = 0, ib = 0, colIb = 0, colIbP1 = ndof;
                beta < ndof;
                beta++, ib += 2, colIb += ndof*2, colIbP1 += ndof*2) {

                for (int alpha = 0, ia = 0; alpha < numCPs; alpha++, ia += 2) {

                    DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
                    DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
                    DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
                    DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
                    DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
                    DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

                    matrixData[colIb + ia] += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
                    matrixData[colIbP1 + ia] += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
                    matrixData[colIb + ia + 1] += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
                    matrixData[colIbP1 + ia + 1] += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];
                }
            }
        }
    }

    Ki = new Matrix(K);
    return K;
}
// getMass is currently not modified for IGA, unchanged
//getMass is a function for dynamic analysis, currently not considered for IGA
const Matrix&
IGAQuad::getMass()
{
	K.Zero();

	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
	  if (rho == 0)
	    rhoi[i] = theMaterial[i]->getRho();
	  else
	    rhoi[i] = rho;	    
	  sum += rhoi[i];
	}

	if (sum == 0.0)
	  return K;

	double rhodvol, Nrho;
    Vector DataJ(2);

    // Loop over the integration points FEMID:1,2,3,4, IGAID:1,2,4,3 
    for (int j = 1; j <= obfs[1]+1; j++) { // loop x-directional GPs
        for (int i = 1; i <= obfs[1]+2; i++) {// loop y-directional GPs
            int CPId = 2 * (j - 1) + i;// for temperaly use
            // Determine Jacobian for this integration point
            DataJ = this->shapeFunction(i, j, eleIdInfo, KnotVect_x, KnotVect_y);

        // Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
            rhodvol *= DataJ(0) * DataJ(1) * rhoi[CPId - 1] * thickness * wts[CPId - 1];

            for (int alpha = 0, ia = 0; alpha < numCPs; alpha++, ia++) {
                Nrho = shp[2][alpha] * rhodvol;
                K(ia, ia) += Nrho;
                ia++;
                K(ia, ia) += Nrho;
            }
        }
	}

	return K;
}

void
IGAQuad::zeroLoad(void)
{
	Q.Zero();

	applyLoad = 0;

	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

  	return;
}

//adding load to the quad element, currently may be not considered
int 
IGAQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	// Added option for applying body forces in load pattern: C.McGann, U.Washington
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		applyLoad = 1;
		appliedB[0] += loadFactor*data(0)*b[0];
		appliedB[1] += loadFactor*data(1)*b[1];
		return 0;
	} else {
		opserr << "IGAQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	} 

	return -1;
}

int 
IGAQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < numCPs; i++) {
    rhoi[i] = theMaterial[i]->getRho();
    sum += rhoi[i];
  }
  
  if (sum == 0.0)
    return 0;
  
  // Get R * accel from the nodes
  static Vector ra(ndof);
  for (i = 1; i <= numCPs; i++) {
      const Vector& Raccel_temp = theNodes[i - 1]->getRV(accel);
      ra[i - 1] = Raccel_temp(0);
      ra[2 * i - 1] = Raccel_temp(1);
  }
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);
  
  //if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
  //    2 != Raccel4.Size()) {
  //  opserr << "IGAQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
  //  return -1;
  //}
  
  //static double ra[8];
  //
  //ra[0] = Raccel1(0);
  //ra[1] = Raccel1(1);
  //ra[2] = Raccel2(0);
  //ra[3] = Raccel2(1);
  //ra[4] = Raccel3(0);
  //ra[5] = Raccel3(1);
  //ra[6] = Raccel4(0);
  //ra[7] = Raccel4(1);
  
  // Compute mass matrix
  this->getMass();
  
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (i = 0; i < ndof; i++)
    Q(i) += -K(i,i)*ra[i];
  
  return 0;
}


//resisting force at the nodes, should be revised for IGA
const Vector&
IGAQuad::getResistingForce()
{
	P.Zero();

	double dvol;
    Vector DataJ(2);//J1,J2 for IGA


    // Loop over the integration points FEMID:1,2,3,4, IGAID:1,2,4,3 
    for (int j = 1; j <= obfs[1]+1; j++) { // loop x-directional GPs
        for (int i = 1; i <= obfs[0]+1; i++) {// loop y-directional GPs
            // Determine Jacobian for this integration point
            DataJ = this->shapeFunction(i, j, eleIdInfo, KnotVect_x, KnotVect_y);
            int CPId = 2 * (j - 1) + i;// for temperaly use
            dvol = DataJ(0) * DataJ(1) * thickness * wts[CPId - 1];// J1*J2*t*W

            // Get material stress response
            const Vector& sigma = theMaterial[CPId - 1]->getStress();

            // Perform numerical integration on internal force
            //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
            //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
            for (int alpha = 0, ia = 0; alpha < numCPs; alpha++, ia += 2) {

                P(ia) += dvol * (shp[0][alpha] * sigma(0) + shp[1][alpha] * sigma(2));

                P(ia + 1) += dvol * (shp[1][alpha] * sigma(1) + shp[0][alpha] * sigma(2));

                // Subtract equiv. body forces from the nodes
                //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
                //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
                if (applyLoad == 0) {
                    P(ia) -= dvol * (shp[2][alpha] * b[0]);
                    P(ia + 1) -= dvol * (shp[2][alpha] * b[1]);
                }
                else {
                    P(ia) -= dvol * (shp[2][alpha] * appliedB[0]);
                    P(ia + 1) -= dvol * (shp[2][alpha] * appliedB[1]);
                }
            }
        }
	}

	// Subtract pressure loading from resisting force
	if (pressure != 0.0) {
		//P = P - pressureLoad;
		P.addVector(1.0, pressureLoad, -1.0);
	}
	
	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;
}

const Vector&
IGAQuad::getResistingForceIncInertia()
{
	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
	  rhoi[i] = theMaterial[i]->getRho();
	  sum += rhoi[i];
	}

	// if no mass terms .. just add damping terms
	if (sum == 0.0) {
	  this->getResistingForce();

	  // add the damping forces if rayleigh damping
	  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	    P += this->getRayleighDampingForces();

	  return P;
	}

	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	const Vector &accel4 = theNodes[3]->getTrialAccel();
	
	static double a[8];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 8; i++)
		P(i) += K(i,i)*a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	  P += this->getRayleighDampingForces();

	return P;
}

int
IGAQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = pressure;

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING IGAQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(3*numCPs);
  
  int i;
  for (i = 0; i < numCPs; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+numCPs) = matDbTag;
    idData(i + 2 * numCPs) = connectedExternalNodes(i - 1);
  }
  
  //idData(8) = connectedExternalNodes(0);
  //idData(9) = connectedExternalNodes(1);
  //idData(10) = connectedExternalNodes(2);
  //idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING IGAQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < numCPs; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING IGAQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
IGAQuad::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING IGAQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  pressure = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  static ID idData(3*numCPs);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING IGAQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }
  for (int i = 1; i <= numCPs; i++) {
      connectedExternalNodes(i - 1) = idData(2*numCPs + i);
  }
  //connectedExternalNodes(0) = idData(8);
  //connectedExternalNodes(1) = idData(9);
  //connectedExternalNodes(2) = idData(10);
  //connectedExternalNodes(3) = idData(11);
  
  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[numCPs];
    if (theMaterial == 0) {
      opserr << "IGAQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < numCPs; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+numCPs);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	    opserr << "IGAQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	    return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "IGAQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < numCPs; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+numCPs);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	    delete theMaterial[i];
	    theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	    if (theMaterial[i] == 0) {
          opserr << "IGAQuad::recvSelf() - material " << i << "failed to create\n";
	      return -1;
	    }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "IGAQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }
  
  return res;
}

void
IGAQuad::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#IGAQuad\n";
    
    int i;
    int numNodes = numCPs;
    const int nstress = 3 ;
    
    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = numCPs;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;
  }
  
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nIGAQuad, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	s << "\tsurface pressure:  " << pressure << endln;
    s << "\tmass density:  " << rho << endln;
    s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < numCPs; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"IGANURBS\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
      for (int i =1 ; i<numCPs;i++){ s << connectedExternalNodes(i-1) << ", "; }
      //s << connectedExternalNodes(1) << ", ";
      //s << connectedExternalNodes(2) << ", ";
      s << connectedExternalNodes(numCPs-1) << "], ";
      s << "\"thickness\": " << thickness << ", ";
      s << "\"surfacePressure\": " << pressure << ", ";
      s << "\"masspervolume\": " << rho << ", ";
      s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
      s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int
IGAQuad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

    // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

    static Vector values(numCPs);

    for (int j=0; j<numCPs; j++)
	   values(j) = 0.0;

    if (displayMode < 4 && displayMode > 0) {
	for (int i=0; i<numCPs; i++) {
	  const Vector &stress = theMaterial[i]->getStress();
	  values(i) = stress(displayMode-1);
	}
    }

    // now  determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    static Matrix coords(numCPs, 3);

    Vector* endCrd_temp = new Vector[numCPs];
    Vector* endDisp_temp = new Vector[numCPs];
    Matrix* eigen_temp = new Matrix[numCPs];
    for (int i = 1; i <= numCPs;i++) {
        endCrd_temp[i-1] = theNodes[i - 1]->getCrds();
        if (displayMode >= 0) {
            endDisp_temp[i-1] = theNodes[i - 1]->getDisp();
            for (int j = 1; j <= 2;j++) {
                coords(i - 1, j - 1) = (endCrd_temp[i-1])(j - 1) + (endDisp_temp[i-1])(j - 1) * fact;
            }
        }
        else {
            int mode = displayMode * -1;
            eigen_temp[i-1] = theNodes[i - 1]->getEigenvectors();
            if ((eigen_temp[i-1]).noCols() >= mode) {
                for (int j = 1; j <= 2;j++) {
                    coords(i - 1, j - 1) = (endCrd_temp[i-1])(j - 1) + (eigen_temp[i-1])(j - 1, mode - 1) * fact;
                }
            }
            else {
                for (int j = 1; j <= 2; j++) {
                    coords(i - 1, j - 1) = (endCrd_temp[i-1])(j - 1);
                }
            }
        }
    }
    delete[] endCrd_temp;
    delete[] endDisp_temp;
    delete[] eigen_temp;
 //   const Vector &end1Crd = theNodes[0]->getCrds();
 //   const Vector &end2Crd = theNodes[1]->getCrds();	
 //   const Vector &end3Crd = theNodes[2]->getCrds();	
 //   const Vector &end4Crd = theNodes[3]->getCrds();	

 //   static Matrix coords(numCPs,3);

 //   if (displayMode >= 0) {    
 //     
 //     const Vector &end1Disp = theNodes[0]->getDisp();
 //     const Vector &end2Disp = theNodes[1]->getDisp();
 //     const Vector &end3Disp = theNodes[2]->getDisp();
 //     const Vector &end4Disp = theNodes[3]->getDisp();

 //     for (int i = 0; i < 2; i++) {
	//coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	//coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
	//coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
	//coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
 //     }
 //   } else {
 //     int mode = displayMode * -1;
 //     const Matrix &eigen1 = theNodes[0]->getEigenvectors();
 //     const Matrix &eigen2 = theNodes[1]->getEigenvectors();
 //     const Matrix &eigen3 = theNodes[2]->getEigenvectors();
 //     const Matrix &eigen4 = theNodes[3]->getEigenvectors();
 //     if (eigen1.noCols() >= mode) {
	//for (int i = 0; i < 2; i++) {
	//  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	//  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	//  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	//  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	//}    
 //     } else {
	//for (int i = 0; i < 2; i++) {
	//  coords(0,i) = end1Crd(i);
	//  coords(1,i) = end2Crd(i);
	//  coords(2,i) = end3Crd(i);
	//  coords(3,i) = end4Crd(i);
	//}    
 //     }
 //   }
    
    int error = 0;

    // finally we  the element using drawPolygon
    error += theViewer.drawPolygon (coords, values, this->getTag());

    return error;
}

Response*
IGAQuad::setResponse(const char **argv, int argc, 
			  OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","IGANURBS");
  output.attr("eleTag",this->getTag());
  for (int i = 1; i <= numCPs; i++) {
      output.attr("nodei", connectedExternalNodes[i-1]);
  }
  //output.attr("node1",connectedExternalNodes[0]);
  //output.attr("node2",connectedExternalNodes[1]);
  //output.attr("node3",connectedExternalNodes[2]);
  //output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=numCPs; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, P);
  }   

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= numCPs) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<numCPs; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(3*numCPs));
  }

  else if ((strcmp(argv[0],"strain") ==0) || (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<numCPs; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(3*numCPs));//4 is the response id
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int 
IGAQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(3*numCPs);
    int cnt = 0;
    for (int i = 0; i < numCPs; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    
    return eleInfo.setVector(stresses);
      
  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(3*numCPs);
    int cnt = 0;
    for (int i = 0; i < numCPs; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}

int
IGAQuad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // quad pressure loading
  if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);
  }
  // a material parameter
  else if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= numCPs)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {

    int matRes = res;
    for (int i=0; i<numCPs; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
	res = matRes;
    }
  }
  
  return res;
}
    
int
IGAQuad::updateParameter(int parameterID, Information &info)
{
	int res = -1;
		int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

	case 1:
		
		for (int i = 0; i<numCPs; i++) {
		matRes = theMaterial[i]->updateParameter(parameterID, info);
		}
		if (matRes != -1) {
			res = matRes;
		}
		return res;
      
	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;

	default: 
	  /*	  
	  if (parameterID >= 100) { // material parameter
	    int pointNum = parameterID/100;
	    if (pointNum > 0 && pointNum <= 4)
	      return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
	    else
	      return -1;
	  } else // unknown
	  */
	    return -1;
  }
}

void IGAQuad::GaussRule(int NGPs1,Vector* Xg, Vector* WXg)
{
    double NGPs = (double)NGPs1;
    Vector X(NGPs), W(NGPs);
    Matrix result(NGPs, 2);
    //Matrix ptss(NGPs, 2);
    //Vector wtss(NGPs);
    //this->GaussRule(NGPs);
    int iter = 2;
    double e1 = NGPs * (NGPs + 1);
    int mm = trunc((NGPs + 1) / 2);
    Vector tt(mm);
    int kk = 1;
    double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
    for (int i = 3; i <= 4 * mm - 1;i += 4) {
        tt(kk - 1) = i * PI / (4 * NGPs + 2);
        kk++;
    }
    Vector x0(tt.Size());
    Vector h(tt.Size()), d1(tt.Size()), pk(tt.Size()), dpn(tt.Size()), d2pn(tt.Size()), d3pn(tt.Size()), d4pn(tt.Size());
    for (int i = 1; i <= tt.Size();i++) {
        x0(i - 1) = (1 - (1 - 1 / NGPs) / (8 * NGPs * NGPs)) * cos(tt(i - 1));
    }
    for (int j = 1; j <= iter;j++) {
        Vector pkm1(x0.Size());
        pkm1.Zero();
        pkm1 += 1;
        pk = x0;
        for (int k = 2; k <= NGPs; k++) {
            Vector pkp1 = 2 * dotProduct(x0, pk) - pkm1 - (dotProduct(x0, pk) - pkm1) / k;
            pkm1 = pk;
            pk = pkp1;
        }
        //Vector Test = dotProduct(x0, x0);
        Vector den = -1 * dotProduct(x0, x0) + 1;
        d1 = NGPs * (pkm1 - dotProduct(x0, pk));
        dpn = dotDivide(d1, den);
        d2pn = dotDivide((2 * dotProduct(x0, dpn) - e1 * pk), den);
        d3pn = dotDivide((4 * dotProduct(x0, d2pn) + (2 - e1) * dpn), den);
        d4pn = dotDivide((6 * dotProduct(x0, d3pn) + (6 - e1) * d2pn), den);
        Vector uu = dotDivide(pk, dpn);
        Vector vv = dotDivide(d2pn, dpn);
        //Initial approximation H
        h = dotProduct(-1 * uu, 0.5 * dotProduct(uu, vv + dotProduct(uu, dotDivide(dotProduct(vv, vv) - dotProduct(uu, d3pn), 3 * dpn))) + 1);
        //Refine H using one step of Newton's method
        Vector p = pk + dotProduct(h, dpn + 0.5 * dotProduct(h, d2pn + dotProduct(h / 3, d3pn + 0.25 * dotProduct(h, d4pn))));
        x0 += h;
    }
    for (int i = 1; i <= x0.Size();i++) { X(i - 1) = -1 * x0(i - 1) - h(i - 1); }
    Vector fx = d1 - dotProduct(h, e1 * (pk + 0.5 * dotProduct(h, dpn + dotProduct(h / 3, d2pn + 0.25 * dotProduct(h, d3pn + 0.2 * dotProduct(h, d4pn))))));
    for (int i = 1;i <= x0.Size();i++) { W(i - 1) = 2 * (1 - X(i - 1) * X(i - 1)) / (fx(i - 1) * fx(i - 1)); }
    if (mm + mm > NGPs) { X(mm - 1) = 0; }
    if (mm + mm != NGPs) { mm = mm - 1; }
    for (int i = 1;i <= mm;i++) {
        X(NGPs - i) = -X(i - 1);
        W(NGPs - i) = W(i - 1);
    }
    for (int i = 1; i <= NGPs; i++) {
        (*Xg)(i - 1) = X(i - 1);
        (*WXg)(i - 1) = W(i - 1);
    }
}
Vector IGAQuad::dotProduct(Vector VecA, Vector VecB)
{
    Vector result(VecA.Size());
    for (int i = 1; i <= VecA.Size();i++) {
        result(i - 1) = VecA(i - 1) * VecB(i - 1);
    }
    return result;
}

Vector IGAQuad::dotDivide(Vector VecA, Vector VecB)
{
    Vector result(VecA.Size());
    for (int i = 1; i <= VecA.Size();i++) {
        result(i - 1) = VecA(i - 1) / VecB(i - 1);
    }
    return result;
}
//The shape function should be revised for IGA, this is the 1d code
void IGAQuad::DerBasisFuns(double Idx, Vector Pts, int obf, int n, Vector KnotVect, Matrix** N0n)
{
    // Idx: span index of the parametric location; Pts: parametric points; obf: order of basis function;
    // n: maximum order of derivatives; KnotVect: knot vector
    int PtsSize = obf+1;
    double iXi;
    double saved;
    double temp;
    int s1, s2;
    double d;
    int j1, j2;
    //N0n.resize(PtsSize,PtsSize);
    //N0n.Zero();
    Matrix Ni(obf + 1, n + 1);
    Matrix ndu(obf + 1, obf + 1);
    Vector left(obf+1); // more 1 space to avoid error
    Vector right(obf+1);
    Matrix a(2, obf + 1);
    for (int i = 1; i <= PtsSize; i++) {
        iXi = Pts(i - 1);
        ndu(0, 0) = 1.0;
        for (int j = 1; j <= obf; j++) {
            left(j) = iXi - KnotVect(Idx - j);
            right(j) = KnotVect(Idx + j - 1) - iXi;
            saved = 0;
            for (int r = 0; r <= j - 1; r++) {
                //lower triangle
                ndu(j, r) = right(r + 1) + left(j - r);
                temp = ndu(r, j - 1) / ndu(j, r);
                //upper triangle
                ndu(r, j) = saved + right(r + 1) * temp;
                saved = left(j - r) * temp;
            }
            ndu(j, j) = saved;
        }
        for (int r = 1; r <= obf + 1;r++) {
            Ni(r - 1, 0) = ndu(r - 1, obf);
        }
        //Ni.setData(ndu(:, p),[],0);//the syntax is wrong
        for (int r = 0; r <= obf; r++) {
            s1 = 0; s2 = 1;
            a(0, 0) = 1;
            for (int k = 1; k <= n; k++) {
                d = 0;
                int rk = r - k;
                int pk = obf - k;
                if (r >= k) {
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                    d = a(s2, 0) * ndu(rk, pk);
                }
                if (rk >= -1) {
                    j1 = 1;
                }
                else {
                    j1 = -rk;
                }
                if ((r - 1) <= pk) {
                    j2 = k - 1;
                }
                else {
                    j2 = obf - r;
                }
                for (int j = j1; j <= j2; j++) {
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                    d += a(s2, j) * ndu(rk + j, pk);
                }
                if (r <= pk) {
                    a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                    d += a(s2, k) * ndu(r, pk);
                }
                Ni(r, k) = d;
                int j = s1;
                s1 = s2;
                s2 = j;//switch rows
            }
        }
        int r = obf;
        for (int k = 1;k <= n;k++) {
            for (int rr = 1; rr <= obf + 1; rr++) {
                Ni(rr - 1, k) *= r;// the syntax is also wrong
            }
            r = r * (obf - k);
        }
        (*N0n)[i-1] = Ni;
    }
}
// to calculate the shape functions and derivatives of Gaussian points
// 需要修改成单个单元的版本！
void IGAQuad::calcDersBasisFunsAtGPs(int obf, int ncp, Vector KnotVect, int d, int NGPs, int Idx, double* J2, Vector* WXg, Matrix** N0n)
{
    // obf: order of basis function, ncp: number of control points;
    // KnotVect: knot vector; d: degree of derivative; NGPs: number of gauss points; Idx: span index of the element;
    Vector Xg(NGPs);
    GaussRule(NGPs,&Xg,WXg);
    //Vector Xg(NGPs);
    
    //for (int i = 1; i <= NGPs; i++) {
    //    Xg(i - 1) = resultGauss(i - 1, 0);
    //    WXg(i - 1) = resultGauss(i - 1, 1);
    //}
    //double J2;
    Vector Xi_e(NGPs);
    //Matrix* N = new Matrix [NGPs];
    //N = new Matrix * [NEl];
    //for (int i = 1; i <= NEl;i++) {
    //    N[i-1] = new Matrix [NGPs];// N = zeros(NEl,NGPs,obf+1,d+1)
    //}
    //int i = Idx[ex - 1];
    *J2 = (KnotVect(Idx) - KnotVect(Idx - 1)) / 2;
    //Vector Xi_e(NGPs);
    for (int ii = 1;ii <= NGPs;ii++) {//有问题！！！
        Xi_e(ii - 1) = (Xg(ii - 1) + 1) * (*J2) + KnotVect(Idx - 1);//Xi(ex,:)=(Xg+1)*J2(ex)+KnotVect(i);
        //Xi_ex(ii - 1) = Xi(ex - 1, ii - 1);
    }
    this->DerBasisFuns(Idx, Xi_e, obf, d, KnotVect, N0n);//N(ex,:,:,:) = DersBasisFuns(i,Xi(ex,:),obf,d,KnotVect)
    //resultDersBasisFunsAtGPS result;// return a struct
    //result.J2 = J2;
    //result.WXg = WXg;
    //result.Xi = Xi_e;
    //result.N = N0n;


}
void IGAQuad::Rationalize(Vector WeightsCP, Vector N0, Matrix N1, Vector* R0, Matrix* R1)
{
    // convert B-spline to NURBS basis functions for 1D
    Vector N0W(N0.Size());
    double W0 = 0;
    for (int i = 1; i <= N0.Size(); i++) {
        N0W(i - 1) = N0(i - 1) * WeightsCP(i - 1);
        W0 += N0W(i - 1);
    }
    *R0 = N0W / W0;
    //for (int i = 1; i <= R0.Size();i++) {
    //    (RatFunc[0])(i - 1, 0) = R0(i - 1);
    //}
    // First derivatives of NURBS basis functions
    Matrix N1W(N1.noRows(), N1.noCols());
    for (int i = 1; i <= N1.noRows(); i++) {
        for (int j = 1;j <= N1.noCols();j++) {
            N1W(i - 1, j - 1) = N1(i - 1, j - 1) * WeightsCP(j - 1);
        }
    }
    Vector W1(N1.noRows());
    Matrix Temp(N1W.noRows(), N1W.noCols());
    for (int i = 1; i <= N1W.noRows();i++) {
        for (int j = 1; j <= N1W.noCols();j++) {
            W1(i - 1) += N1W(i - 1, j - 1);
        }
    }
    for (int i = 1; i <= N1W.noRows();i++) {
        for (int j = 1; j <= N1W.noCols();j++) {
            Temp(i - 1, j - 1) = (*R0)(j - 1) * W1(i - 1);
        }
    }
    Matrix N1WR0W1 = N1W - Temp;
    *R1 = N1WR0W1 / W0;

    //result.R0 = R0;
    //result.R1 = R1;
}
//double IGAQuad::shapeFunction1(int obfs, int ncps, Vector* KnotVects, int nds, int NGPss, int Idxs)
//{
//    //calculate the shape function and derivatives of x and y directions
//    int obf_x = obfs[0]; int obf_y = obfs[1];
//    int ncp_x = ncps[0]; int ncp_y = ncps[1];
//    Vector KnotVect_x = KnotVects[0]; Vector KnotVect_y = KnotVects[1];
//    int nd_x = nds[0]; int nd_y = nds[1];
//    int NGPs_x = NGPss[0]; int NGPs_y = NGPss[1];
//    int Idx_x = Idxs[0]; int Idx_y = Idxs[1];
//
//    resultDersBasisFunsAtGPS ValueBFG_x = calcDersBasisFunsAtGPs(obf_x, ncp_x, KnotVect_x, nd_x, NGPs_x, Idx_x);// return basis function and derivatives of GP
//    resultDersBasisFunsAtGPS ValueBFG_y = calcDersBasisFunsAtGPs(obf_y, ncp_y, KnotVect_y, nd_y, NGPs_y, Idx_y);
//    double Jx = ValueBFG_x.J2; Vector Wx = ValueBFG_x.WXg; Matrix* Nx = ValueBFG_X.N;
//    double Jy = ValueBFG_y.J2; Vector Wy = ValueBFG_y.WXg; Matrix* Ny = ValueBFG_X.N;
//    int NEN = (obf_x + 1) * (obf_y + 1); // number of local basic functions
//    Vector N0(NEN);
//    Matrix N1(2, NEN);
//    double WeightsCP1[] = { 1,1 };
//    Vector WeightsCP(WeightsCP1, 2);
//    const Vector& nd1Crds = theNodes[0]->getCrds();
//    const Vector& nd2Crds = theNodes[1]->getCrds();
//    const Vector& nd3Crds = theNodes[2]->getCrds();
//    const Vector& nd4Crds = theNodes[3]->getCrds();
//    Matrix CtrlPts(4, 3);
//    for (int i = 1; i <= 4; i++) {
//        for (int j = 1; j <= 3;j++) {
//            Vector& ndCrds = theNodes[i - 1]->getCrds();
//            CtrlPts(i - 1, j - 1) = *ndCrds(j - 1);
//        }
//    }
//    Matrix Ke(NEN * 2, NEN * 2);
//    for (int qy = 1; qy <= NGPs_y;qy++) {
//        for (int qx = 1;qx <= NGPs_x;qx++) {
//            int k = 1;
//            for (int j = 1;j <= obf_y + 1;j++) {
//                for (int i = 1; i <= obf_x + 1;i++) {
//                    N0(k - 1) = (Nx[qx - 1])(i - 1, 0) * (Ny[qy - 1])(j - 1, 0);//shape function
//                    N1(0, k - 1) = (Nx[qx - 1])(i - 1, 1) * (Ny[qy - 1])(j - 1, 0);// 1st derivatives
//                    N1(1, k - 1) = (Nx[qx - 1])(i - 1, 0) * (Ny[qy - 1])(j - 1, 1);
//                    k++;
//                }
//            }
//        }
//    }
//    //R1 = Rationalize()
//}
// calculate the shape function and derivatives of the quadrature point
int IGAQuad::findSpan(int order, int ncp, int eleId, int nele,Vector KnotVect) {
    ID Idxs(nele);
    Idxs.Zero();
    int iE = 1;
    for (int i = order; i <= ncp; i++) {
        if (abs(KnotVect(i - 1) - KnotVect(i)) > 1.4901e-08) {
            Idxs(iE - 1) = i;
            iE = iE + 1;
        }
    }
    int Idx = Idxs(eleId - 1);
    return Idx;
}
Vector IGAQuad::shapeFunction(int qx, int qy, ID eleIdInfo, Vector KnotVect_x, Vector KnotVect_y){
    //需要的输入：int obfs[2]: x,y方向的基函数阶次；int mults[2]：x,y方向重叠节点个数，目前可以假定为1
    // Vector KnotVect_x/y: x/y方向的KnotVector; Vector Weights_x/y: x/y方向的控制点权重；int NumIntegPoints[2]: x,y方向的高斯点数目
    //const Vector& nd1Crds = theNodes[0]->getCrds();
    //const Vector& nd2Crds = theNodes[1]->getCrds();
    //const Vector& nd3Crds = theNodes[2]->getCrds();
    //const Vector& nd4Crds = theNodes[3]->getCrds();
    int ex = eleIdInfo(0); int nex = eleIdInfo(1);
    int ey = eleIdInfo(2); int ney = eleIdInfo(3);
    Matrix CtrlPts(numCPs, 2);
    Vector* ndCrds_temp = new Vector[numCPs];
    for (int i = 1; i <= numCPs; i++) {
        ndCrds_temp[i-1] = theNodes[i - 1]->getCrds();
        for (int j = 1; j <= 2;j++) {
            CtrlPts(i - 1, j - 1) = (ndCrds_temp[i-1])(j - 1);
        }
    }
    delete[] ndCrds_temp;

    int nds[] = { 1,1 }; // maximum order of derivatives
    //int obfs[] = { 1,1 }; // order of basis function
    int ncps[] = { obfs[0]+1,obfs[1]+1 };//number of control points
    //double KnotVect1[] = { 0,0,1,1 };
    //Vector KnotVect(KnotVect1, 4);
    //Vector* KnotVects = new Vector[4];
    //for (int i = 1; i <= 2; i++) {
    //    KnotVects[i - 1] = KnotVect;// assume knot vector in x-y direction are same
    //}
    int NGPss[] = { obfs[0]+1,obfs[1]+1 };
    //calculate the shape function and derivatives of x and y directions
    int obf_x = obfs[0]; int obf_y = obfs[1];
    int ncp_x = ncps[0]; int ncp_y = ncps[1];
    //Vector KnotVect_x = KnotVects[0]; Vector KnotVect_y = KnotVects[1];
    int nd_x = nds[0]; int nd_y = nds[1];
    int NGPs_x = NGPss[0]; int NGPs_y = NGPss[1];
    int Idx_x = findSpan(obfs[0], ncp_x, ex, nex, KnotVect_x);
    int Idx_y = findSpan(obfs[1], ncp_y, ey, ney, KnotVect_y);// find the knot span

    double Jx=0., Jy=0.;
    Vector Wx(NGPs_x), Wy(NGPs_y);
    //Matrix* N0nx, N0ny;
    calcDersBasisFunsAtGPs(obf_x, ncp_x, KnotVect_x, nd_x, NGPs_x, Idx_x,&Jx,&Wx,&N0nx);// return basis function and derivatives of GP
    calcDersBasisFunsAtGPs(obf_y, ncp_y, KnotVect_y, nd_y, NGPs_y, Idx_y,&Jy,&Wy,&N0ny);
    //double Jx = ValBFG_x.J2; Vector Wx = ValBFG_x.WXg; Matrix* Nx = ValBFG_x.N;
    //double Jy = ValBFG_y.J2; Vector Wy = ValBFG_y.WXg; Matrix* Ny = ValBFG_y.N;
    //int NEN = (obf_x + 1) * (obf_y + 1); // number of local basic functions
    Vector N0(numCPs);
    Matrix N1(2, numCPs);
    Vector WeightsCP(ndof+1);
    for (int i = 1; i <= ndof;i++) { WeightsCP(i - 1) = 1.0; }
    //Matrix CtrlPts(4, 2);
    //CtrlPts(0, 0) = nd1Crds(0); CtrlPts(0, 1) = nd1Crds(1);//FEM  43   IGA  34
    //CtrlPts(1, 0) = nd2Crds(0); CtrlPts(1, 1) = nd2Crds(1);//     12        12
    //CtrlPts(2, 0) = nd4Crds(0); CtrlPts(2, 1) = nd4Crds(1);//
    //CtrlPts(3, 0) = nd3Crds(0); CtrlPts(3, 1) = nd3Crds(1);    
    //Matrix Ke(, NEN * 2);
    Vector DataJ(2+1);//J1 and J2 of the quadrature point
    int k = 1;
    for (int j = 1;j <= obf_y + 1;j++) {
        for (int i = 1; i <= obf_x + 1;i++) {
            N0(k - 1) = (N0nx[qx - 1])(i - 1, 0) * (N0ny[qy - 1])(j - 1, 0);//shape function
            N1(0, k - 1) = (N0nx[qx - 1])(i - 1, 1) * (N0ny[qy - 1])(j - 1, 0);// 1st derivatives
            N1(1, k - 1) = (N0nx[qx - 1])(i - 1, 0) * (N0ny[qy - 1])(j - 1, 1);
            k++;
        }
    }
    //delete[] Nx;
    //delete[] Ny;
    Vector R0;
    Matrix R1;
    Rationalize(WeightsCP, N0, N1, &R0, &R1);//obtaining the rationalized shape function

    double J2 = Jx * Jy; // J2
    DataJ(1) = J2;
    double W = Wx(qx-1) * Wy(qy-1);
    wts[qx - 1 + (qy - 1) * (obfs[1] + 1)] = W;
    //gradient of mapping from parametrical space to physical space
    Matrix dxdxi(2, 2);
    dxdxi.Zero();
    for (int i = 1; i <= R1.noRows(); i++) {
        for (int j = 1; j <= CtrlPts.noCols(); j++) {
            for (int k = 1; k <= R1.noCols();k++) {
                dxdxi(i-1, j-1) += R1(i-1, k-1) * CtrlPts(k-1, j-1);
            }
        }
    }
    double J1 = abs(dxdxi(0, 0) * dxdxi(1, 1) - dxdxi(0, 1) * dxdxi(1, 0));// J1
    DataJ(0) = J1;
    //static Matrix dxdxi1 = dxdxi;
    Matrix dxdxiInv;
    dxdxi.Invert(dxdxiInv);// inverse of jacobian matrix

    Matrix dRdx(dxdxi.noRows(), R1.noCols());
    dRdx.Zero();
    for (int i = 1; i <= dxdxiInv.noRows();i++) {
        for (int j = 1; j <= R1.noCols();j++) {
            for (int k = 1; k <= R1.noRows();k++) {
                dRdx(i-1,j-1) += dxdxiInv(i - 1, k - 1)* R1(k - 1, j - 1);
            }
        }
    }
    //int IGAId[] = { 1,2,4,3 };
    for (int i = 1;i <= R0.Size();i++) {
        //int ii = IGAId[i - 1];
        shp[2][i - 1] = R0(i - 1);
        shp[0][i - 1] = dRdx(0, i - 1);
        shp[1][i - 1] = dRdx(1, i - 1);
    }
    return DataJ;
}
void 
IGAQuad::setPressureLoadAtNodes(void)
{
        pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();
	const Vector &node4 = theNodes[3]->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);
	double x4 = node4(0);
	double y4 = node4(1);

	double dx12 = x2-x1;
	double dy12 = y2-y1;
	double dx23 = x3-x2;
	double dy23 = y3-y2;
	double dx34 = x4-x3;
	double dy34 = y4-y3;
	double dx41 = x1-x4;
	double dy41 = y1-y4;

	double pressureOver2 = pressure/2.0;

	// Contribution from side 12
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(2) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(2) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*dy23;
	pressureLoad(3) += pressureOver2*-dx23;
	pressureLoad(5) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(4) += pressureOver2*dy34;
	pressureLoad(6) += pressureOver2*dy34;
	pressureLoad(5) += pressureOver2*-dx34;
	pressureLoad(7) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(6) += pressureOver2*dy41;
	pressureLoad(0) += pressureOver2*dy41;
	pressureLoad(7) += pressureOver2*-dx41;
	pressureLoad(1) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}








