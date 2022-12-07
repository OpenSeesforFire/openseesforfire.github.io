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

                                                                        
#ifndef BrickThermalAction_h
#define BrickThermalAction_h
                                                                  
 //Modified by Liming Jiang [http://openseesforfire.github.io]

#include <ElementalLoad.h>
#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>
#include <Element.h>

class BrickThermalAction : public ElementalLoad
{
  public:
  // Constructors based on 9, 2 or 0 temperature changes given
  BrickThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                double t3, double locY3, double t4, double locY4,
                double t5, double locY5, double t6, double locY6,
                double t7, double locY7, double t8, double locY8,
                double t9, double locY9, 
		        int theElementTag);


    BrickThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                int theElementTag);
	BrickThermalAction(int tag, 
					 double locY1, double locY2,
					 TimeSeries* theSeries,int theElementTag
					 );

	BrickThermalAction(int tag, int theElementTag);

  BrickThermalAction();    

  ~BrickThermalAction();
  
  const Vector &getData(int &type, double loadFactor);
  virtual void applyLoad(const Vector &factors);
  virtual void applyLoad(double loadFactor);
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double Temp[9]; //Initial Temperature 
  double TempApp[9]; // Temperature applied
  double Loc[9]; // Location through the depth of section
  static Vector data; // data for temperature and locations
  int ThermalActionType;

  int indicator; //indicator if fireloadpattern was called
  Vector Factors;
  TimeSeries* theSeries;
  //TimeSeries* theSeries1;
};

#endif

