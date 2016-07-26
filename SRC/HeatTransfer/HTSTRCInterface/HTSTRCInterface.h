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


// Written: Yaqiang Jiang (y.jiang@ed.ac.uk)
// Created: March 2012
// Revision: A

#ifndef HTSTRCInterface_h
#define HTSTRCInterface_h

#include <vector>

class KDTree;
class HeatTransferDomain;

class HTSTRCInterface {

public:

  HTSTRCInterface(int numNodes, int numSteps, int numDim, const char* filename);

  HTSTRCInterface(int numNodes, int numDim, HeatTransferDomain* theDomain);
  
  virtual ~ HTSTRCInterface();

  //void getTemperature(double x, double y, double y);
  //void getTemperature(double x, double y);
  
  // return temperature for a fiber, perform searching in 2D
  const std::vector<double>& getFiberTemperature(double x1, double x2, 
						                         double y1, double y2);

  //return temperature for a fiber, perform searching in 3D
  const std::vector<double>& getFiberTemperature(double x1, double x2, 
						                         double y1, double y2,
						                         double z1, double z2);

  const std::vector<int>& getPointsInRegion2D(double x1, double x2, 
						                         double y1, double y2);
  const std::vector<int>& getPointsInRegion3D(double x1, double x2, 
						                         double y1, double y2,
						                         double z1, double z2);


  int getClosestPoint2D(double x, double y);

  int getClosestPoint3D(double x, double y, double z);
  
  //save the temperature data to a file for a specified region
  void getTemperature2DRegion(double x1, double x2, double y1, double y2);
  void getTemperature3DRegion(double x1, double x2, double y1, double y2,
	                          double z1, double z2);

  //return 9 temperatures for the I section and 20 for the slab;
  //however you can get as many as you can depending on the meshing 
  //density in heat transfer analysis
  void getTemperatureForCompositeSection(double h, double t, double s, double b,
	                                     double slab_w, double slab_d);
  //return 9 temperatures for the I section
  void getTemperatureForISection(double h, double t, double s, double b);


  //return 9 temperatures for the I section and 20 for slab with 
  //TimeSeries format
  void getTemperatureSeriesForCompositeSection(double h, double t, double s, double b,
	                                           double slab_w, double slab_d);
  //return 9 temperatures for the I section with TimeSeries format
  void getTemperatureSeriesForISection(double h, double t, double s, double b);

  //extract the temperature dataset at time step ntsp
  void getTemperatureSnapShot(int ntsp, double dtime);


private:
    int regionCnt2d;
	int regionCnt3d;
	const char* fileName;
	double** datastored;
	std::vector<double>* Temp_history;
	KDTree* mytree;
	int numnodes, numsteps,numdimensions;
};

#endif
