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
// $Date: 2006-08-03 23:42:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/JointEPMaterialThermal.h,v $
                                                                        
#ifndef JointEPMaterialThermal_h
#define JointEPMaterialThermal_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// JointEPMaterialThermal. JointEPMaterialThermal provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) JointEPMaterialThermal.h, revA"

#include <UniaxialMaterial.h>

class JointEPMaterialThermal : public UniaxialMaterial
{
  public:
    JointEPMaterialThermal(int tag, double E, double eyp,int soft=0);    
    JointEPMaterialThermal(int tag, double E, double eyp, double eyn, double ezero=0.0, int soft =0);    
    JointEPMaterialThermal();    

    ~JointEPMaterialThermal();

    const char *getClassType(void) const {return "JointEPMaterialThermal";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return E;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    double getThermalElongation(void);//returning the Thermal Elongation
    double getElongTangent(double, double&, double&, double);
    int getVariable(const char* variable, Information&);  //For providing the universal function
    
  protected:
    
  private:
    double fyp, fyn;	// positive and negative yield stress
    double ezero;	// initial strain
    double E;		// elastic modulus
    double E0;     //initial elastic modulus
    double ep;		// plastic strain at last commit

    int softIndex;
    double ThermalElongation;
    double Temp;


    double trialStrain;	     // current trial strain
    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
    double commitStrain;     // last commited strain
    double commitStress;     // last commited stress
    double commitTangent;    // last committed  tangent
};


#endif



