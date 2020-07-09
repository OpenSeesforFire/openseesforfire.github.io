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
// $Date: 2003-02-14 23:01:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/JointEPMaterialThermal.cpp,v $
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) JointEPMaterialThermal.C, revA"


#include <JointEPMaterialThermal.h>
#include <Vector.h>
#include <Channel.h>
#include <Parameter.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_JointEPMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3 || numArgs > 5) {
    opserr << "Invalid #args,  want: JointEPMaterialThermal $tag $E $epsP <$epsN $eps0>\n";
    return 0;
  }
  
  int iData[1];
  double dData[4];
  dData[3] = 0.0; // setting default eps0 to 0.

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for JointEPMaterialThermal ElasticPP" << endln;
    return 0;
  }

  numData = numArgs-1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for JointEPMaterialThermal " << iData[0] << endln;
    return 0;	
  }
        
  if (numData == 2) 
    dData[2] = -dData[1];

  int softindex = 1;
  // Parsing was successful, allocate the material
    theMaterial = new JointEPMaterialThermal(iData[0], dData[0], dData[1], dData[2], dData[3], softindex);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type JointEPMaterialThermal\n";
    return 0;
  }

  return theMaterial;
}


JointEPMaterialThermal::JointEPMaterialThermal(int tag, double e, double eyp, int soft)
:UniaxialMaterial(tag,MAT_TAG_ElasticPPMaterial),
 ezero(0.0), E(e), ep(0.0), E0(e), softIndex(soft),
 trialStrain(0.0), trialStress(0.0), trialTangent(E),
 commitStrain(0.0), commitStress(0.0), commitTangent(E)
{
  fyp = E*eyp;
  fyn = -fyp;
}

JointEPMaterialThermal::JointEPMaterialThermal(int tag, double e, double eyp,
				     double eyn, double ez, int soft  )
:UniaxialMaterial(tag,MAT_TAG_ElasticPPMaterial),
 ezero(ez), E(e), E0(e),  ep(0.0), softIndex(soft),
 trialStrain(0.0), trialStress(0.0), trialTangent(E),
 commitStrain(0.0), commitStress(0.0), commitTangent(E)
{
    if (eyp < 0) {
	opserr << "JointEPMaterialThermal::JointEPMaterialThermal() - eyp < 0, setting > 0\n";
	eyp *= -1.;
    }
    if (eyn > 0) {
	opserr << "JointEPMaterialThermal::JointEPMaterialThermal() - eyn > 0, setting < 0\n";
	eyn *= -1.;
    }    
    
    fyp = E*eyp;
    fyn = E*eyn;
}

JointEPMaterialThermal::JointEPMaterialThermal()
:UniaxialMaterial(0,MAT_TAG_ElasticPPMaterial),
 fyp(0.0), fyn(0.0), ezero(0.0), E(0.0), E0(0.0), ep(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
 commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

JointEPMaterialThermal::~JointEPMaterialThermal()
{
  // does nothing
}

int 
JointEPMaterialThermal::setTrialStrain(double strain, double strainRate)
{
  /*
    if (fabs(trialStrain - strain) < DBL_EPSILON)
      return 0;
  */
    trialStrain = strain;

    double sigtrial;	// trial stress
    double f;		// yield function

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    //sigtrial  = E * trialStrain;
    //sigtrial -= E * ezero;
    //sigtrial -= E *  ep;

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f <= fYieldSurface ) {

      // elastic
      trialStress = sigtrial;
      trialTangent = E;

    } else {

      // plastic
      if ( sigtrial > 0.0 ) {
	trialStress = fyp;
      } else {
	trialStress = fyn;
      }

      trialTangent = 0.0;
    }

    return 0;
}

double 
JointEPMaterialThermal::getStrain(void)
{
  return trialStrain;
}

double 
JointEPMaterialThermal::getStress(void)
{
  return trialStress;
}


double 
JointEPMaterialThermal::getTangent(void)
{
  return trialTangent;
}

int 
JointEPMaterialThermal::commitState(void)
{
    double sigtrial;	// trial stress
    double f;		// yield function

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f > fYieldSurface ) {
      // plastic
      if ( sigtrial > 0.0 ) {
	ep += f / E;
      } else {
	ep -= f / E;
      }
    }

    commitStrain = trialStrain;
    commitTangent=trialTangent;
    commitStress = trialStress;

    return 0;
}	


int 
JointEPMaterialThermal::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;

  return 0;
}


int 
JointEPMaterialThermal::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialTangent = commitTangent = E;
  trialStress = commitStress = 0.0;

  ep = 0.0;
  
  return 0;
}


UniaxialMaterial *
JointEPMaterialThermal::getCopy(void)
{
  JointEPMaterialThermal *theCopy =
    new JointEPMaterialThermal(this->getTag(),E,fyp/E,fyn/E,ezero);
  theCopy->ep = this->ep;
  
  return theCopy;
}


int 
JointEPMaterialThermal::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = ep;
  data(2) = E;
  data(3) = ezero;
  data(4) = fyp;
  data(5) = fyn;
  data(6) = commitStrain;
  data(7) = commitStress;
  data(8) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "JointEPMaterialThermal::sendSelf() - failed to send data\n";

  return res;
}

int 
JointEPMaterialThermal::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "JointEPMaterialThermal::recvSelf() - failed to recv data\n";
  else {
    this->setTag(int(data(0)));
    ep    = data(1);
    E     = data(2);
    ezero = data(3);
    fyp   = data(4);
    fyn   = data(5);  
    commitStrain=data(6);
    commitStress=data(7);
    commitTangent=data(8);
    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
  }

  return res;
}

void 
JointEPMaterialThermal::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "JointEPMaterialThermal tag: " << this->getTag() << endln;
		s << "  E: " << E << endln;
		s << "  ep: " << ep << endln;
		s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
	}
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"JointEPMaterialThermal\", ";
		s << "\"E\": " << E << ", ";
		s << "\"epsp\": " << ep << "}";
	}
}


int
JointEPMaterialThermal::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(fyp);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"epsP") == 0 || strcmp(argv[0],"ep") == 0) {
    param.setValue(ep);
    return param.addObject(3, this);
  }

  return -1;
}

int
JointEPMaterialThermal::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->fyp = info.theDouble;
    this->fyn = -fyp;
    break;
  case 2:
    this->E = info.theDouble;
    trialTangent = E;
    break;
  case 3:
    this->ep = info.theDouble;
    break;
  default:
    return -1;
  }
  
  return 0;
}



//Liming for updating the reduction factors////////////start
double
JointEPMaterialThermal::getElongTangent(double TempT, double& ET, double& Elong, double TempTmax)
{
    double ThermalElongation;
    Temp = TempT;
    double* redfactors;

    if (softIndex == 0) {
        //doNothing
    }
    else if (softIndex == 1) {
        redfactors = new double[12];
        double SteelRedfactors[12] = { 1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225 , 0.0 };
        for (int i = 0; i < 12; i++)
            redfactors[i] = SteelRedfactors[i];
    }
    else if (softIndex == 2) {
        redfactors = new double[12];
        double ConcreteRedfactors[12] = { 0.625, 0.4318 ,0.3036, 0.1875 ,0.1, 0.045, 0.03, 0.015, 0.008 , 0.004,0.001,0.0 };
        for (int i = 0; i < 12; i++)
            redfactors[i] = ConcreteRedfactors[i];
    }
    else {
        opserr << "ElasticIsotropic3DThermal " << this->getTag() << " recieves an invalid softening index" << endln;
    }

    if (softIndex ==1) {

        for (int i = 0; i < 13; i++) {
            if (Temp <= 80 + 100 * i)
            {
                if (i == 0) {
                    E = E0 * (1.0 - Temp * (1.0 - redfactors[0]) / 80);
                    E = E0 * (1.0 - Temp * (1.0 - redfactors[0]) / 80);
                }
                else if (i == 12) {
                    opserr << "Warning:The temperature " << Temp << " for SteelECthermal is out of range\n";
                    return -1;
                }
                else {
                    E = E0 * (redfactors[i - 1] - (Temp + 20 - 100 * i) * (redfactors[i - 1] - redfactors[i]) / 100);
                    E = E0 * (redfactors[i - 1] - (Temp + 20 - 100 * i) * (redfactors[i - 1] - redfactors[i]) / 100);
                }
                break;
            }

        }

        ThermalElongation = 0;
        Elong = ThermalElongation;

    }
    else {
        ET = E0;
        ThermalElongation = 0;
    }
    ET = E;
    Elong = ThermalElongation;
    return 0;
}

double
JointEPMaterialThermal::getThermalElongation(void)
{
    return ThermalElongation;
}



int
JointEPMaterialThermal::getVariable(const char* variable, Information& info)
{
    if (strcmp(variable, "ThermalElongation") == 0) {
        info.theDouble = ThermalElongation;
        return 0;
    }
    else if (strcmp(variable, "ElongTangent") == 0) {
        Vector* theVector = info.theVector;
        if (theVector != 0) {
            double tempT, ET, Elong, TempTmax;
            tempT = (*theVector)(0);
            ET = (*theVector)(1);
            Elong = (*theVector)(2);
            TempTmax = (*theVector)(3);
            this->getElongTangent(tempT, ET, Elong, TempTmax);
            (*theVector)(0) = tempT;
            (*theVector)(1) = ET;
            (*theVector)(2) = Elong;
            (*theVector)(3) = TempTmax;
        }
        return 0;
    }
    else if (strcmp(variable, "TempAndElong") == 0) {
        Vector* theVector = info.theVector;
        if (theVector != 0) {
            (*theVector)(0) = Temp;
            (*theVector)(1) = ThermalElongation;
        }
        else {
            opserr << "null Vector in EC" << endln;
        }

        return 0;
    }
    return -1;
}