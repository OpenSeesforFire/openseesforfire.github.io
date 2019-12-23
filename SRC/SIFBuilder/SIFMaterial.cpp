
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
                                                                   
/**********************************************************************
** This project is aiming to provide a Tcl interface to define models  **
** for simulating structural behaviours under fire action.           **
** Developed by:                  `                                   **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $
// This file constructs the class SIFMaterial for storing and generating material models
// Created by Liming Jiang @UoE

#include <SIFMaterial.h>
#include <CarbonSteelEC3.h>
#include <ConcreteEC2.h>
#include <Concrete02Thermal.h>
#include <SimpleMaterial.h>
#include <HeatTransferMaterial.h>
#include <UniaxialMaterial.h>
#include <SteelECThermal.h>
#include <Steel01Thermal.h>
#include <StainlessSteelEC.h>
#include <StainlessECThermal.h>
#include <ConcreteECThermal.h>
#include <ElasticMaterialThermal.h>

SIFMaterial::SIFMaterial(int tag, int MaterialTypeTag, double fy, double E0):TaggedObject(tag), 
MaterialTypeTag(MaterialTypeTag),MaterialPars(2),theHTMaterial(0)
{
	MaterialPars(0)= fy;
	MaterialPars(1)= E0;
#ifdef _DEBUG
	opserr<<"Successfully added material "<<tag<<endln
		<<"Material type: Steel "<<MaterialTypeTag <<" ,fy: "<<fy<<" , E0: "<<E0<<endln;
#endif
}

SIFMaterial::SIFMaterial(int tag, int MaterialTypeTag, double fy, double fu, double E0) :TaggedObject(tag),
MaterialTypeTag(MaterialTypeTag), MaterialPars(3), theHTMaterial(0)
{
	MaterialPars(0) = fy;
	MaterialPars(1) = E0;
	MaterialPars(2) = fu;
#ifdef _DEBUG
	opserr << "Successfully added material " << tag << endln
		<< "Material type: Steel " << MaterialTypeTag << " ,fy: " << fy << " , E0: " << E0 << endln;
#endif
}

SIFMaterial::SIFMaterial(int tag, int MaterialTypeTag, 
						 double fc, double epsc0, double fcu, double epscu, double rat, double ft,
						 double Ets, double moisture):TaggedObject(tag),
						 MaterialTypeTag(MaterialTypeTag),MaterialPars(8),theHTMaterial(0)
{
  MaterialPars(0)= fc;
  MaterialPars(1)= epsc0;
  MaterialPars(2)= fcu;
  MaterialPars(3)= epscu;
  MaterialPars(4)= rat;
  MaterialPars(5)= ft;
  MaterialPars(6)= Ets;
  MaterialPars(7)= moisture;
}

SIFMaterial::~SIFMaterial()
{
	if (theHTMaterial!=0)
		delete theHTMaterial;

}

int 
SIFMaterial::getMaterialTypeTag()
{
	return MaterialTypeTag;
}


const Vector& 
SIFMaterial::getMaterialPars(void)
{
		return 0;
}


HeatTransferMaterial*
SIFMaterial::getHeatTransferMaterial()
{
	int tag = this->getTag();
	if(MaterialTypeTag==130){
    //Steel material, default EC thermal properties will be defined
    theHTMaterial = new CarbonSteelEC3(tag);
  }
	else if (MaterialTypeTag == 301) {
		//Steel material, default EC thermal properties will be defined
		//theHTMaterial = new CarbonSteelEC3(tag);
		theHTMaterial = new StainlessSteelEC(tag);
	}
	else if(MaterialTypeTag==220){
		if (MaterialPars.Size()!=8) {
		opserr<<"WARNING:: SIFMaterial fialed to create heat transfer material, moisture doesn't exsist"<<endln;
		return 0;
    }
    theHTMaterial = new ConcreteEC2(tag, MaterialPars(7));
  }
  
  return theHTMaterial;
}

UniaxialMaterial*
SIFMaterial::getUniaxialMaterial(bool isElastic)
{
	//isElastic = true;
	int tag = this->getTag();
	UniaxialMaterial* theuniMaterial= 0;
  if(MaterialTypeTag==130){
    //Steel material, default EC thermal properties will be defined
	double fy = MaterialPars(0);
	double E0= MaterialPars(1);
	if (isElastic)
		theuniMaterial = new ElasticMaterialThermal(tag, E0, 1.2e-5);
	else
		//theuniMaterial = new SteelECThermal(tag,3,fy, E0);
		theuniMaterial = new Steel01Thermal(tag, fy, E0, 0.05);
		 // Steel01Thermal(int tag, double fy, double E0, double b,
  }
  else if (MaterialTypeTag == 301) {
	  //StainlessECThermal(int tag, int grade, double Fy, double E, double Fu,double sigInit)
	  double fy = MaterialPars(0);
	  double fu = MaterialPars(2);
	  double E0 = MaterialPars(1);
	  if (isElastic)
		  theuniMaterial = new ElasticMaterialThermal(tag, E0, 1.2e-5);
	  else
		  theuniMaterial = new StainlessECThermal(tag, 3, fy, E0, fu);

	  //Grade14301:1, Grade14401:2, Grade14571:3, Grade14003:4, Grade14462:5
  }
  else if(MaterialTypeTag==220){
   //ConcreteECThermal::ConcreteECThermal(int tag, double _fc, double _epsc0, double _fcu,
				     //double _epscu, double _rat, double _ft, double _Ets):
	  //isElastic = true;
	    double fc= MaterialPars(0);
		double epsc0 = MaterialPars(1) ;
		double fcu = MaterialPars(2) ;
		double epsc = MaterialPars(3) ;
		double rat =MaterialPars(4) ;
		double ft = MaterialPars(5) ;
		double Ets = MaterialPars(6) ;
		double moisture =MaterialPars(7) ;
		if(isElastic)
		theuniMaterial = new ElasticMaterialThermal(tag, 3e10,1.4e-5, true);
		else
		theuniMaterial = new ConcreteECThermal(tag,fc, epsc0, fcu, epsc, rat, ft, Ets);
		//theuniMaterial = new Concrete02Thermal(tag, fc, epsc0, fcu, epsc, rat, ft, Ets);
		// ElasticMaterialThermal(int tag, double E, double alpha, double eta = 0.0);
		//theuniMaterial = new ElasticMaterialThermal(tag,3e10, 1.4e-5);
		
  }


  if(theuniMaterial ==0){
	opserr<<"SIFMaterial failed to generatre uniaxailMaterial with tag: "<<tag<<endln;
  }
  return theuniMaterial;
}


double 
SIFMaterial::getInitialModulus(){
	double E0;

	if(MaterialTypeTag==130){
		E0 = MaterialPars(1);
		}

	return E0;
	
}