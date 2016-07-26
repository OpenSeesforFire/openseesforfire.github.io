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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterialThermal.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterialThermal. 
//
// What: "@(#) ElasticMaterialThermal.C, revA"

#include <ElasticMaterialThermal.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_NewElasticMaterialThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? alpha?<eta?> ... " << endln;
    return 0;
  }
  
  int iData[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData >= 3) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;	
    }
  } 
  else if (numData ==2) {
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
	}
	dData[2]= 0.0;
 }
 else {
	 opserr << "Invalid #args, want: uniaxialMaterial Elastic " << iData[0] <<"E? alpha?" << endln;
  }

 //opserr<< "receieved alpha"<<dData[1]<<endln;

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMaterialThermal(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticMaterialThermal\n";
    return 0;
  }

  return theMaterial;
}


ElasticMaterialThermal::ElasticMaterialThermal(int tag, double e, double alpha, double et)
:UniaxialMaterial(tag,MAT_TAG_ElasticMaterialThermal),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(e),Alpha(alpha),eta(et), parameterID(0),ThermalElongation(0),TempT(0)
{

}

ElasticMaterialThermal::ElasticMaterialThermal()
:UniaxialMaterial(0,MAT_TAG_ElasticMaterialThermal),
 trialStrain(0.0),  trialStrainRate(0.0),
 E(0.0),Alpha(0.0), eta(0.0), parameterID(0),ThermalElongation(0),TempT(0)
{

}

ElasticMaterialThermal::~ElasticMaterialThermal()
{
  // does nothing
}

int 
ElasticMaterialThermal::setTrialStrain(double strain, double strainRate)
{
    opserr << "ElasticMaterial::setTrialStrain (double strain, double strainRate) - should never be called\n";
  return 0;
}

int 
ElasticMaterialThermal::setTrialStrain(double strain, double FiberTemperature, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;

    return 0;
}


int 
ElasticMaterialThermal::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;

    stress = E*trialStrain + eta*trialStrainRate;
    tangent = E;

    return 0;
}



double 
ElasticMaterialThermal::getStress(void)
{
    return E*trialStrain + eta*trialStrainRate;
}


int 
ElasticMaterialThermal::commitState(void)
{
    return 0;
}


int 
ElasticMaterialThermal::revertToLastCommit(void)
{
    return 0;
}


int 
ElasticMaterialThermal::revertToStart(void)
{
    trialStrain      = 0.0;
    trialStrainRate  = 0.0;
    return 0;
}

UniaxialMaterial *
ElasticMaterialThermal::getCopy(void)
{
    ElasticMaterialThermal *theCopy = new ElasticMaterialThermal(this->getTag(),E,Alpha,eta);
    theCopy->trialStrain     = trialStrain;
    theCopy->trialStrainRate = trialStrainRate;
    return theCopy;
}

int 
ElasticMaterialThermal::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = Alpha;
  data(3) = eta;
  

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMaterialThermal::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticMaterialThermal::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(4);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ElasticMaterialThermal::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag(data(0));
    E   = data(1);
	Alpha = data(2);
    eta = data(3);
  }
    
  return res;
}

void 
ElasticMaterialThermal::Print(OPS_Stream &s, int flag)
{
    s << "Elastic tag: " << this->getTag() << endln;
	s << "  E: " << E << " eta: " << eta << "alpha: "<< Alpha<< endln;
}

int
ElasticMaterialThermal::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);
  
  else if (strcmp(argv[0],"eta") == 0)
    return param.addObject(2, this);

  return -1;
}

int 
ElasticMaterialThermal::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    eta = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
ElasticMaterialThermal::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
ElasticMaterialThermal::getStressSensitivity(int gradIndex, bool conditional)
{
  if (parameterID == 1)
    return trialStrain;
  else if (parameterID == 2)
    return trialStrainRate;
  else
    return 0.0;
}

double
ElasticMaterialThermal::getInitialTangentSensitivity(int gradIndex)
{
  return 0.0;
}

int
ElasticMaterialThermal::commitSensitivity(double strainGradient,
				   int gradIndex, int numGrads)
{
  // Nothing to commit ... path independent
  return 0.0;
}


//Liming for updating the reduction factors////////////start
double 
ElasticMaterialThermal::getElongTangent(double tempT, double &ET, double &Elong, double TempTmax) 
{
  ThermalElongation = tempT*Alpha;
  TempT = tempT;
  //ThermalElongation = 0 ;   //debug  Liming
  ET = E;   
  Elong = ThermalElongation;

  return 0;
}

double 
ElasticMaterialThermal::getThermalElongation(void) 
{
  return ThermalElongation;
}

int 
ElasticMaterialThermal::getVariable(const char *variable, Information &info)
{
  if (strcmp(variable,"ThermalElongation") == 0) {
    info.theDouble = ThermalElongation;    
    return 0;
  } 
  else if (strcmp(variable,"ElongTangent") == 0) {
    Vector *theVector = info.theVector;
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
  else if (strcmp(variable,"TempAndElong") == 0) {
    Vector *theVector = info.theVector;
	if (theVector!= 0) {
		(*theVector)(0) = TempT;
        (*theVector)(1) = ThermalElongation;
	}else{
		opserr<<"null Vector in EC"<<endln;
	}

	return 0;
  }
  return -1;
}