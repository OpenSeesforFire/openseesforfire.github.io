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
                                                                        

  //Based on JZ, JJ @Edinburgh University                                                                    
 //Modified by Liming Jiang [http://openseesforfire.github.io]


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <ShellThermalAction.h>
#include <Vector.h>

Vector ShellThermalAction::data(18);
ShellThermalAction::ShellThermalAction(int tag, 
                         double t1, double locY1, double t2, double locY2,
                         double t3, double locY3, double t4, double locY4,
                         double t5, double locY5, double t6, double locY6,
                         double t7, double locY7, double t8, double locY8,
                         double t9, double locY9, 
			       int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_ShellThermalAction, theElementTag), 
	ThermalActionType(LOAD_TAG_ShellThermalAction),theSeries(0),theSeries1(0), Temp(9), TempApp(9), Loc(9)
{
    Temp(0)=t1; Temp(1) = t2; Temp(2) = t3; Temp(3) = t4; Temp(4) = t5;
	Temp(5)=t6; Temp(6) = t7; Temp(7) = t8; Temp(8) = t9; 

    Loc(0)=locY1; Loc(1) = locY2; Loc(2) = locY3; Loc(3) = locY4; Loc(4) = locY5;
    Loc(5)=locY6; Loc(6) = locY7; Loc(7) = locY8; Loc(8) = locY9;

  TempApp.Zero();
  indicator=1; //without path timeseries defined;
}



ShellThermalAction::ShellThermalAction(int tag, 
                         double t1, double locY1, double t2, double locY2, 
			       int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_ShellThermalAction, theElementTag), 
ThermalActionType(LOAD_TAG_ShellThermalAction),theSeries(0), theSeries1(0), Temp(9), TempApp(9),Loc(9)
   
{
	Temp = Vector(9);
	Temp(0) = t1;
	Temp(8) = t2;
	Loc = Vector(9);
	Loc(0)=locY1;
	Loc(8)=locY2;
	for(int i= 1; i<8;i++){
		Temp(i)=Temp(0)-i*(Temp(0)-Temp(8))/8;
		Loc(i)=Loc(0)-i*(Loc(0)-Loc(8))/8;
	}

	TempApp.Zero();

	indicator=1; //without path timeseries defined;
}


ShellThermalAction::ShellThermalAction(int tag,
					 double locY1, double locY2,
					 TimeSeries* theSeries,int theElementTag
					 )
:ElementalLoad(tag, LOAD_TAG_ShellThermalAction, theElementTag),theSeries(theSeries),  theSeries1(0),
ThermalActionType(LOAD_TAG_ShellThermalAction),Loc(9),Temp(9),TempApp(9)
{
  Loc(0)=locY1;
  Loc(8)=locY2;
  
  for(int i= 1; i<8;i++){
		Loc(i)=Loc(0)-i*(Loc(0)-Loc(8))/8;
	}
  Temp.Zero();
  TempApp.Zero();

	indicator=2 ;// Independent timeseries were created;

}


//for composite section using two time series
ShellThermalAction::ShellThermalAction(int tag,
	double locY1, double locY2, double locY3, double locY4,
	TimeSeries* theSeries, TimeSeries* theSeries1, int theElementTag
)
	:ElementalLoad(tag, LOAD_TAG_ShellThermalAction, theElementTag), theSeries(theSeries), theSeries1(theSeries1),
	ThermalActionType(LOAD_TAG_ShellThermalAction), Temp(0), TempApp(18), Loc(18)
{
	Loc(0) = locY1;
	Loc(8) = locY2;
	Loc(9) = locY3;
	Loc(17) = locY4;

	for (int i = 1; i < 8; i++) {
		Loc(i) = Loc(0) - i * (Loc(0) - Loc(8)) / 8;
		Loc(i+9) = Loc(9) - i * (Loc(9) - Loc(17)) / 8;
	}

	data.resize(36);
	TempApp.Zero();

	indicator = 4;// Two timeseries were created for composite section;

}

ShellThermalAction::ShellThermalAction(int tag,  
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_ShellThermalAction, theElementTag),
  ThermalActionType(LOAD_TAG_NodalThermalAction),theSeries(0), theSeries1(0)
{
	Loc.Zero();
	Temp.Zero();
	TempApp.Zero();

    indicator=3 ;// USing Nodal Thermal Action;
}

ShellThermalAction::~ShellThermalAction()
{
     if(theSeries!=0)
		 theSeries=0;
}

const Vector &
ShellThermalAction::getData(int& type, double loadFactor)
{
	type = ThermalActionType;
	//0,2,4,6,8...16 storing temperature;
	//1,3,5,7,9...17 storing local y location;
	// opserr << "TempPP" << TempApp << endln;
	if (indicator != 4) {
		for (int i = 0; i < 9; i++) {
			data(2 * i) = TempApp(i);
			data(2 * i + 1) = Loc(i);
		}
	}
	else {
		//two time series
		for (int i = 0; i < 18; i++) {
			data(2 * i) = TempApp(i);
			data(2 * i + 1) = Loc(i);
		}
	}
	
	return data;
}

void 
ShellThermalAction::applyLoad(const Vector &factors) 
{
	for(int i=0; i<9 ;i++) {
		TempApp(i)= Temp(i)*factors(i);
	}
	if (theElement != 0)
    theElement->addLoad(this, factors(0));
}

void 
ShellThermalAction::applyLoad(double loadfactor)
{
	// first determine the load factor
	if (indicator==2) {
		TempApp =((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		// opserr << "TempPP" << TempApp << endln;
		   //PathTimeSeriesThermal returns absolute temperature;
	}
	else if (indicator == 4) {
		Vector vec1 = ((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		Vector vec2 = ((PathTimeSeriesThermal*)theSeries1)->getFactors(loadfactor);
		//PathTimeSeriesThermal returns absolute temperature;
		for (int i = 0; i < 9; i++) {
			TempApp(i) = vec1(i);
			TempApp(i + 9) = vec2(i);
		}
	}
	else{
		  TempApp=Temp*loadfactor;
	}
	if (theElement != 0)
    theElement->addLoad(this, loadfactor);
}



int 
ShellThermalAction::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
ShellThermalAction::recvSelf(int commitTag, Channel &theChannel,  
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}

// do it later
void 
ShellThermalAction::Print(OPS_Stream &s, int flag)
{
  s << "ShellThermalAction - reference load : " <<Temp[0] <<" change  temp of bot\n";
  s <<  Temp[8] << " change  temp at top\n";
  s << "  element acted on: " << eleTag << endln;
}

