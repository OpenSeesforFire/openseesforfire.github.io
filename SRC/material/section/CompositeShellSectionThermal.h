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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-21 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/CompositeShellSectionThermal.h,v $

// To define shell section for composite floor slabs
//
// Composite Shell Section
//
// Added for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>



class CompositeShellSectionThermal : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    CompositeShellSectionThermal( ) ;

    //full constructor
    CompositeShellSectionThermal(   int tag, SectionForceDeformation* theSec1, SectionForceDeformation* theSec2,
                                double ratio1,double ratio2, double ribAngle, int distribution);

    const char *getClassType(void) const {return "CompositeShellSectionThermal";};

    //destructor
    virtual ~CompositeShellSectionThermal( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;

    //mass per unit area
    double getRho() ;

    //send back order of strain in vector form
    int getOrder( ) const ;

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    //send back order of strain in vector form
    const ID& getType( ) ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //send back the initial tangent 
    const Matrix& getInitialTangent( ) {return this->getSectionTangent();}

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	const Vector &getTemperatureStress(const Vector&); //Added by LMJ
	double determineFiberTemperature(const Vector& DataMixed, double fiberLoc) ; //Added by LMJ



  private :

    //quadrature data
    double Ratio1;
    double Ratio2;

    double ribAng;
    int distributionMode;

    SectionForceDeformation* theSection1;  //pointers to the first section
    SectionForceDeformation* theSection2;  //pointers to the second section

    Vector strainResultant ;

    static Vector stressResultant ;

    static Matrix tangent ;


    Vector* sT;
    double* ThermalElongation;  // thermal elongation

    Vector thermalForce1; //added by QJ
    Vector thermalForce2; //added by QJ

    Vector thermalDeformation1; //added by QJ
    Vector thermalDeformation2; //added by QJ

    double axialStrain; //added by QJ
    double curvatureStrain; //added by QJ

    Matrix tangentSection1; //added by QJ
    Matrix tangentSection2; //added by QJ

    Matrix inverseTangentSection1; //added by QJ
    Matrix inverseTangentSection2; //added by QJ
    
    double commitAxial1; //added by QJ
    double commitCurve1; //added by QJ
    double commitAxial2; //added by QJ
    double commitCurve2; //added by QJ
    double commitAxial; //added by QJ
    double commitCurvature; //added by QJ
    double axialStrain1;
    double curvatureStrain1;
    double axialStrain2;
    double curvatureStrain2;
 
} ; //end of CompositeShellSectionThermal declarations





