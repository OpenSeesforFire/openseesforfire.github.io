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

#ifndef GmshBuilder_h
#define GmshBuilder_h

#include <HTModelBuilder.h>

class HeatTransferElement;
class HeatTransferNode;
class HeatTransferDomain;
class ParametricFireEC1;
class UserDefinedFire;
class AlpertCeilingJetModel;
class BoundaryPattern;
class FireModel;

/*
Read the mesh infomation generated Gmsh. Element type tag: 1 for 2-node line element,
2 for 3-node triangle, 3 for 4-node quadrangle, 5 for 8-node hexahedron(brick), 8 for
3-node second order line, 16 for 8-node second order quadrangle.
*/


class GmshBuilder : public HTModelBuilder
{
  public:
    GmshBuilder(HeatTransferDomain& theDomain, const char* filename);
    ~GmshBuilder();    

    int BuildModel(void);

	// impose the same convection boundary conditions for all the boundaries
	int setConvectionBC(double h, double Tf);
	// impose a specific convection boundary condition for bounaries with
	// the same tag
	int setConvectionBC(double h, double Tf, int tag);

	// impose the same convection boundary conditions for all the boundaries
	int setRadiationBC(double epsilon, double sigma, 
					   double alpha, double qir);
    // impose a specific radiation boundary condition for bounaries with
	// the same tag
	int setRadiationBC(double epsilon, double sigma, 
					   double alpha, double qir, int tag);

	int setDirichletBC(double T, int tag);

	int setParametricFireConvcBC(ParametricFireEC1*  themodel, double h, int tag);

	int setParametricFireRadBC(ParametricFireEC1*  themodel, double epsilon, 
		                       double sigma, double alpha,
							   int tag);


	//int setParametricFireConvcBC(double I, double Av, double H, double At,
	//	double Af, double Qf, double tlim,
	//	double h, int tag);

	//int setParametricFireRadBC(double I, double Av, double H, double At,
	//	double Af, double Qf, double tlim,
	//	double epsilon, double sigma, double alpha,
	//	int tag);

	int setLocalisedFireBC(double crd1, double crd2, double crd3, 
		                   double D, double Q, double H, 
						   int centerLineTag, int tag);

	int setUDFFireConvcBC(UserDefinedFire* themodel, double h, int tag);
	int setUDFFireRadBC(UserDefinedFire* themodel, double epsilon, 
		                double sigma, double alpha, int tag);
	int setUDFFlux(UserDefinedFire* themodel, int tag);

	int setAlpertFireConvcBC(AlpertCeilingJetModel* themodel, double h, int tag);
	int setAlpertFireRadBC(AlpertCeilingJetModel* themodel, double epsilon, 
		                   double sigma, double alpha, int tag);

    // Tf is the flame temperautre
	int setLinearTravellingFireBC(AlpertCeilingJetModel* themodel, double pos1, double pos2,
		                          double Lmax, double h, double Ta, double Tf, double epsilon, 
								  double sigma, double alpha, int DirectionFlag, int tag);

	int createBoundaryPattern(void); 
	int createFireImposedPattern(FireModel* themodel);
	void removePattern(int tag);

  protected:

  private:

	int** PhysicalID;
	int** ConnectivityBC;
	int** ConnectivityDomain;
	int** BoundaryEleID;
	int numFluxBC;
	int numPatterns, numNearFieldPattern, numFarFieldPattern;
	int TotalElements, NumBoundaryEle, NumDomainEle, NUMPhysicalNames;
	const char* fileName;
	bool travellingStamp, TravelStop;
};

#endif