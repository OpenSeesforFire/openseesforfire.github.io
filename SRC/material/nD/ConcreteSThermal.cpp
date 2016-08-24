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
// $Date: 2012-05-24 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ConcreteSThermal.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Simple Concrete Model, plane stress
// 


#include <ConcreteSThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>

//null constructor
ConcreteSThermal::ConcreteSThermal( ) : 
NDMaterial(0, ND_TAG_ConcreteSThermal ), 
strain0(3), strain(3), stress0(3), stress(3), stressd(3),
tangent(3,3),eTangent(3,3),
cStrain0(0),cStrain(0),
TempAndElong(2)
{ }


//full constructor
ConcreteSThermal::ConcreteSThermal(int tag, double rE, double rnu, double rfc, double rft, double rEs) :
NDMaterial( tag, ND_TAG_ConcreteSThermal ),
strain0(3), strain(3), stress0(3), stress(3), stressd(3),
tangent(3,3),eTangent(3,3),
cStrain0(0),cStrain(0),
E0(rE), nu(rnu), fc0(fabs(rfc)), ft0(rft), Es0(fabs(rEs)),
TempAndElong(2)
{
  fc = fc0;
  ft = ft0;
  E=E0;
  Es=Es0;

  setInitials();
}


//destructor
ConcreteSThermal::~ConcreteSThermal( ) 
{ 

} 

void ConcreteSThermal::setInitials()
{
  double fac;
    strain0.Zero();
  strain.Zero();
  stress0.Zero();
  stress.Zero();
  eTangent.Zero();
  eTangent(0,0) = 1.0;
  eTangent(0,1) = nu;
  eTangent(1,0) = nu;
  eTangent(1,1) = 1.0;
  eTangent(2,2) = 0.5 * (1.0 - nu);
  fac = E / (1.0 - nu * nu);
  eTangent *= fac;
  tangent = eTangent;
  Ep = 0.5* E;
  EmEp1 = 1.0 / (E + Ep);
  beta = - (Ep -2* Es) / ft;

//  Ep = E * Es / (E + Es)
//  EmEp1 = 1.0 / (E - Ep) = (E + Es) / E^2;
//  EmEp1 = (E + Es) / (E * E);
//  beta = -(E + 2 * Es) / ft;
//  Es = -E/2;
//  EmEp1 = 0.5 / E;
  ftmin = 1.0e-2 * ft;
}


void ConcreteSThermal::setTempInitials()
{
  double fac;
  eTangent.Zero();
  eTangent(0,0) = 1.0;
  eTangent(0,1) = nu;
  eTangent(1,0) = nu;
  eTangent(1,1) = 1.0;
  eTangent(2,2) = 0.5 * (1.0 - nu);
  fac = E / (1.0 - nu * nu);
  eTangent *= fac;
  tangent = eTangent;
  Ep = 0.5* E;
  EmEp1 = 1.0 / (E + Ep);
  beta = - (Ep -2* Es) / ft;

  //Ep = 10 * E;
  //EmEp1 = 1.0 / (E + Ep);
  //beta = - (Ep + 11 * Es) / ft;


//  Ep = E * Es / (E + Es)
//  EmEp1 = 1.0 / (E - Ep) = (E + Es) / E^2;
//  EmEp1 = (E + Es) / (E * E);
//  beta = -(E + 2 * Es) / ft;
//  Es = -E/2;
//  EmEp1 = 0.5 / E;
  ftmin = 1.0e-2 * ft;
}

//make a clone of this material
NDMaterial*
ConcreteSThermal::getCopy( ) 
{
  ConcreteSThermal *clone ;   //new instance of this class

  clone = new ConcreteSThermal( this->getTag(), E, nu, fc, ft, Es);

  return clone ;
}


//make a clone of this material
NDMaterial* 
ConcreteSThermal::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
ConcreteSThermal::getOrder( ) const
{
  return 3 ;
}


const char*
ConcreteSThermal::getType( ) const 
{
  return "ConcreteSThermal" ; 
}



//swap history variables
int 
ConcreteSThermal::commitState( ) 
{
  stress0 = stress;
  strain0 = strain;
  cStrain0 = cStrain;
  return 0;
}


//revert to last saved state
int 
ConcreteSThermal::revertToLastCommit( )
{
  return 0;
}


//revert to start
int
ConcreteSThermal::revertToStart( )
{
  strain0.Zero();
  strain.Zero();
  stress0.Zero();
  stress.Zero();
  cStrain0 = 0.0;
  cStrain = 0.0;
  return 0;
}

//Set TemperatureAndElongation

double
ConcreteSThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{

  double Temp = tempT;
  double epsc0 = -0.0025;
  if (Temp <= 80) {
    ft = ft0;
  }
  else if (Temp <= 580) {
    ft = (1.0 - 1.0*(Temp -80)/500)*ft0;
    //Ets = (1.0 - 1.0*(Temp -80)/500)*fc0 * 1.5 / epsc0T;
    //Ets = (1.0 - 1.0*(Temp -80)/500)*EtsT;
	if(ft<0.01*ft0)
		ft = 0.01*ft0;
  }
  else {
    ft = 0.01*ft0;
    //Ets = 1.0e-10;
    //ft = 0;
    //Ets = 0;
  }
  
  // compression strength, at elevated temperature
  //   strain at compression strength, at elevated temperature
  //   ultimate (crushing) strain, at elevated temperature
  if (Temp <= 0) {
    fc = fc0;
    epsc0 = -0.0025;
    //fcu = fcuT;
    //epscu = -0.02;
    //Ets = EtsT;  jz what is there the statement?
  }
  else if (Temp <= 80) {
    fc = fc0;
    epsc0 = -(0.0025 + (0.004-0.0025)*(Temp - 0)/(80 - 0));
  }
  else if (Temp <= 180) {
    fc = fc0*(1 - (Temp - 80)*0.05/100);
    epsc0 = -(0.0040 + (0.0055-0.0040)*(Temp - 80)/100);
  }
  else if (Temp <= 280) {
    fc = fc0*(0.95 - (Temp - 180)*0.1/100);
    epsc0 = -(0.0055 + (0.0070-0.0055)*(Temp - 180)/100);
   
  }
  else if (Temp <= 380) {
    fc = fc0*(0.85 - (Temp - 280)*0.1/100);
    epsc0 = -(0.0070 + (0.0100-0.0070)*(Temp - 280)/100);
  }
  else if (Temp <= 480) {
    fc = fc0*(0.75 - (Temp - 380)*0.15/100);
    epsc0 = -(0.0100 + (0.0150-0.0100)*(Temp - 380)/100);

  }
  else if (Temp <= 580) {
    fc = fc0*(0.60 - (Temp - 480)*0.15/100);
    epsc0 = -(0.0150 + (0.0250-0.0150)*(Temp - 480)/100);
  }
  else if (Temp <= 680) {
    fc = fc0*(0.45 - (Temp - 580)*0.15/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 780) {
    fc = fc0*(0.30 - (Temp - 680)*0.15/100);
    epsc0 = -0.0250;

  }
  else if (Temp <= 880) {
    fc = fc0*(0.15 - (Temp - 780)*0.07/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 980) {
    fc = fc0*(0.08 - (Temp - 880)*0.04/100);
    epsc0 = -0.0250;
  }
  else if (Temp <= 1080) {
    fc = fc0*(0.04 - (Temp - 980)*0.03/100);
    epsc0 = -0.0250;
  }
  else  {
    opserr << "the temperature is invalid\n";
    
  }
  
  double ThermalElongation=0;
  if (Temp <= 1) {
		  ThermalElongation = Temp  * 9.213e-6;
  }
  else if (Temp <= 680) {
    ThermalElongation = -1.8e-4 + 9e-6 *(Temp+20) + 2.3e-11 *(Temp+20)*(Temp+20)*(Temp+20);
  }
  else if (Temp <= 1180) {
    ThermalElongation = 14e-3;
  }
  else {
    opserr << "the temperature is invalid\n";
  }
  
  E = -fc*0.0025/epsc0/fc0*E0;
  Es = -fc*0.0025/epsc0/fc0*Es0;
  Elong = ThermalElongation;
  this->setTempInitials();
 // E = E0;
  //Es = Es0;
  TempAndElong(0) = Temp;
  TempAndElong(1) = E;
  ET =E;

return 0;
}



//receive the strain
int
ConcreteSThermal::setTrialStrain( const Vector &strainFromElement )
{
  static Matrix CfCf(3,3);
  static Vector flow(3), Cf(3);
  double vStress, vStress1, yieldFunc;
  double fCf, sigm, sigd, theta;
  double ps1, ps2, psmax, tStrain, eps;
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;
//#ifdef _DEBUG
  //if(TempAndElong(0)>10)
	 // opserr<<"strain  "<<strain<<endln<<"strain0  "<<strain0<<endln;
//#endif
  stress = stress0 + eTangent * (strain - strain0);
  
  tangent = eTangent;
  
  vStress = sqrt(  stress(0) * stress(0)
                 - stress(0) * stress(1)
                 + stress(1) * stress(1)
                 + stress(2) * stress(2) * 3.0);
  
  yieldFunc = vStress - fc;
  sigm = 0.5 * (stress(0) + stress(1));
  sigd = 0.5 * (stress(0) - stress(1));
  ps2 = sqrt(sigd * sigd + stress(2) * stress(2));
  ps1 = sigm + ps2;
  ps2 = sigm - ps2;
  if (yieldFunc > 0.0 && ps1 <= ft && ps2 <= ft)
  {
    vStress1 = 1.0 / vStress;
    flow(0) = (stress(0) - 0.5 * stress(1)) * vStress1;
    flow(1) = (stress(1) - 0.5 * stress(0)) * vStress1;
    flow(2) = 3.0 * stress(2) * vStress1;
    Cf = eTangent * flow;
    fCf = flow^Cf;
    eps = yieldFunc / fCf;
    stress -= eTangent * (eps * flow);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        CfCf(i,j) = Cf(i) * Cf(j);
    tangent -= CfCf / (fCf + 1.e-3 * E);
  }
  sigm = 0.5 * (stress(0) + stress(1));
  sigd = 0.5 * (stress(0) - stress(1));
  theta = atan2(stress(2), sigd);
  ps2 = sqrt(sigd * sigd + stress(2) * stress(2));
  ps1 = sigm + ps2;
  ps2 = sigm - ps2;
  
  if (ps2 > 0.0 && ps1 < -fc) ps1 = -fc;
  if (ps1 > 0.0 && ps2 < -fc) ps2 = -fc;
  psmax = ft + Ep * cStrain0;
//  if (psmax < ftmin) psmax = ftmin;
  yieldFunc = ps1 - psmax;

  if (yieldFunc > 0.0)
  {
    cStrain = cStrain0 + yieldFunc * EmEp1;
    ps1 = ft + Ep * cStrain;
  }

  yieldFunc = ps2 - psmax;

  if (yieldFunc > 0.0)
  {
    tStrain = cStrain0 + yieldFunc * EmEp1;
    ps2 = ft + Ep * tStrain;
    if (cStrain < tStrain) cStrain = tStrain;
  }
  
  sigm = 0.5 * (ps1 + ps2);
  sigd = 0.5 * fabs(ps1 - ps2);
  
  stress(1) = sigd * cos(theta);
  stress(0) = sigm + stress(1);
  stress(1) = sigm - stress(1);
  stress(2) = sigd * sin(theta);

  stressd = stress;   
  
  double damage = exp(beta * cStrain);
  if (ps1 > 0.0){
	  ps1 *= damage;
	  //tangent(0,0) *= (-damage);
	  //tangent *= 0.5;
	  }
  if (ps2 > 0.0) {
	  ps2 *= damage;
	  //tangent(1,1) *= (-damage);
	  //tangent *= 0.5;
	  }

  sigm = 0.5 * (ps1 + ps2);
  sigd = 0.5 * fabs(ps1 - ps2);
  
stressd(1) = sigd * cos(theta);
  stressd(0) = sigm + stressd(1);
  stressd(1) = sigm - stressd(1);
  stressd(2) = sigd * sin(theta);

//if (ps1 > 0.0 || ps2 > 0.0)
//{
//    double damage = exp(beta * cStrain);
//    stressd *= damage;
  
  // tangent(1,1) *= damage;
//}
#ifdef _DEBUG
  //if(TempAndElong(0)>10)
	//  opserr<<stressd<<endln;
#endif
  return 0;
}


//send back the strain
const Vector& 
ConcreteSThermal::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
ConcreteSThermal::getStress( )
{
  return stressd ;
}


//send back the tangent 
const Matrix&  
ConcreteSThermal::getTangent( )
{
  return tangent ;
}

const Matrix&  
ConcreteSThermal::getInitialTangent
( )
{
  return eTangent ;
}


//print out data
void  
ConcreteSThermal::Print( OPS_Stream &s, int flag )
{
  s << "ConcreteSThermal Material tag: " << this->getTag() << endln ; 
  s << "  E:  " << E  << " ";
  s << "  nu: " << nu << " ";
  s << "  fc: " << fc << " ";
  s << "  ft: " << ft << " ";
  s << "  Es: " << Es << " ";

}


int 
ConcreteSThermal::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0, cnt = 0;

  static Vector data(13);

  data(cnt++) = this->getTag();
  data(cnt++) = E;
  data(cnt++) = nu;
  data(cnt++) = fc;
  data(cnt++) = ft;
  data(cnt++) = Es;
  data(cnt++) = cStrain0;

  int i;
  for (i = 0; i < 3; i++) 
    data(cnt++) = strain0(i);

  for (i = 0; i < 3; i++) 
    data(cnt++) = stress0(i);

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "ConcreteSThermal::sendSelf() - failed to send data" << endln;

   return res;
}

int 
ConcreteSThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0, cnt = 0;

  static Vector data(13);

  res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "ConcreteSThermal::recvSelf -- could not recv Vector" << endln;
   return res;
  }

  this->setTag(int(data(cnt++)));
  E = data(cnt++);
  nu = data(cnt++);
  fc = data(cnt++);
  ft = data(cnt++);
  Es = data(cnt++);
  cStrain0 = data(cnt++);

  setInitials();

  int i;
  for (i = 0; i < 3; i++)
    strain0(i) = data(cnt++);

  for (i = 0; i < 3; i++)
    stress0(i) = data(cnt++);

  return res;
}



//send back TempAndElong(Liming,UoE)
const Vector&
ConcreteSThermal::getTempAndElong( void)
{
    return TempAndElong;
}
 
