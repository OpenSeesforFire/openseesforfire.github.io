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

#include <HTSTRCInterface.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HT_NodeIter.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <KDTree.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;


HTSTRCInterface::HTSTRCInterface(int numNodes, int numSteps, int numDim, const char* filename)
:numnodes(numNodes), numsteps(numSteps), numdimensions(numDim), datastored(0), 
 fileName(filename), mytree(0), Temp_history(0),regionCnt2d(0),regionCnt3d(0)
{
    if(datastored)
		delete[]datastored;

	int col = numnodes * numsteps;
	datastored = new double*[col];

	
	ifstream inputFile(fileName, ios::in);
	if (!inputFile) {
		cout << "FATAL:HTSTRCInterface::HTSTRCInterfacec( ) -";
		cout << " could not open file: " << fileName << endl;
		exit(-1);
		}

	if(numdimensions == 2) {

		for(int i=0;i<col;i++)
			{
			datastored[i] = new double[4];
			for (int j=0;j<4;j++)
				{
				inputFile >> datastored[i][j];
				//cout << datastored[i][j]<< " ";
				}
			//cout << endl;
			}

		} else if(numdimensions == 3) {

			for(int i=0;i<col;i++)
				{
				datastored[i] = new double[5];
				for (int j=0;j<5;j++)
					{
					inputFile >> datastored[i][j];
					//cout << datastored[i][j]<< " ";
					}
				//cout << endl;
				}
		}

    int cnt = 0;
	double* setpoints;
	if(numdimensions == 2){

		setpoints = new double[numnodes*2];
		for(int i=0;i<numnodes;i++){
			setpoints[cnt] = datastored[i][1];
			setpoints[++cnt] = datastored[i][2];
			cnt++;
			}

		}else if (numdimensions == 3){

			setpoints = new double[numnodes*3];
			for(int i=0;i<numnodes;i++){
				setpoints[cnt] = datastored[i][1];
				setpoints[++cnt] = datastored[i][2];
				setpoints[++cnt] = datastored[i][3];
				cnt++;
				}
		}

	if(mytree)
		delete mytree;

	mytree = new KDTree(setpoints,numnodes,numdimensions);

}


HTSTRCInterface::HTSTRCInterface(int numNodes, int numDim, HeatTransferDomain* theDomain)
:numnodes(numNodes), numsteps(0), numdimensions(numDim), datastored(0), 
 fileName(0), mytree(0), Temp_history(0),regionCnt2d(0),regionCnt3d(0)
{
	if(theDomain == 0) {
		cout << "FATAL:HTSTRCInterface::HTSTRCInterfacec( ) -";
		cout << " no HeatTransferDomain is associated " << endl;
		exit(-1);
		}

	HT_NodeIter& theDomainNodes = theDomain->getNodes();
	HeatTransferNode* theNode;
   
	int cnt = 0;
	double* setpoints;

	if(numdimensions == 2) {
		setpoints = new double[numnodes*2];
		while ((theNode = theDomainNodes()) != 0) {
			const Vector& crds = theNode->getCrds();
			setpoints[cnt] = crds(0);
			setpoints[++cnt] = crds(1);
			cnt++;
			}
		}else if(numdimensions == 3) {
			setpoints = new double[numnodes*3];
			while ((theNode = theDomainNodes()) != 0) {
				const Vector& crds = theNode->getCrds();
				setpoints[cnt] = crds(0);
				setpoints[++cnt] = crds(1);
				setpoints[++cnt] = crds(2);
				cnt++;
				}
		}
	

	if(mytree)
		delete mytree;

	mytree = new KDTree(setpoints,numNodes,numdimensions);	
	
	}


HTSTRCInterface::~HTSTRCInterface()
{
    if(datastored)
		delete[] datastored;

	if(mytree)
		delete mytree;
}


const vector<double>&  
HTSTRCInterface::getFiberTemperature(double x1, double x2, 
						             double y1, double y2)
{
  if(numdimensions != 2){
	  cerr << "Error:HTSTRCInterface::getFiberTemperature(x1,x2,y1,y2)--"
		   << " dimension not matching" << endl;
	  exit(-1);
	}

  Range* therange = new Range[2];
  therange[0] = new double[2];
  therange[1] = new double[2];

  therange[0][0] = x1;
  therange[0][1] = x2;
  therange[1][0] = y1;
  therange[1][1] = y2;

  vector<int> pointsInRange = mytree->get_points_in_range(therange);

  // compute the averaged temperature
  int numpts = pointsInRange.size();

  if(!Temp_history)
	  Temp_history = new vector<double>(2*numsteps);

  for(int stps=0;stps<numsteps;stps++){
	  double Tsum = 0;
	  for(int i=0;i<numpts;i++){
		  int n = pointsInRange[i] + stps * numnodes;
		  Tsum = Tsum + datastored[n][3];
		  }
      //half of the vector stores temperature
	  (*Temp_history)[stps] = Tsum / (double)numpts;
	  //another half of the vector stores time
	  (*Temp_history)[stps+numsteps] = datastored[stps*numnodes][0];
	  }

  return *Temp_history;
}


const vector<double>&  
HTSTRCInterface::getFiberTemperature(double x1, double x2, 
						             double y1, double y2,
									 double z1, double z2)
{
  if(numdimensions != 3){
	  cerr << "Error:HTSTRCInterface::getFiberTemperature(x1,x2,y1,y2,z1,z2)--"
		   << " dimension not matching" << endl;

	  exit(-1);
	}

  Range* therange = new Range[3];
  therange[0] = new double[2];
  therange[1] = new double[2];
  therange[2] = new double[2];

  therange[0][0] = x1;
  therange[0][1] = x2;
  therange[1][0] = y1;
  therange[1][1] = y2;
  therange[1][0] = z1;
  therange[1][1] = z2;

  vector<int> pointsInRange = mytree->get_points_in_range(therange);

  // compute the averaged temperature
  int numpts = pointsInRange.size();

  if(!Temp_history)
	  Temp_history = new vector<double>(2*numsteps);

  for(int stps=0;stps<numsteps;stps++){
	  double Tsum = 0;
	  for(int i=0;i<numpts;i++){
		  int n = pointsInRange[i] + stps * numnodes;
		  Tsum = Tsum + datastored[n][3];
		  }
	  (*Temp_history)[stps] = Tsum / numpts;
	  (*Temp_history)[stps+numsteps] = datastored[stps*numnodes][0];
	  }

  return *Temp_history;
}

const vector<int>& 
HTSTRCInterface::getPointsInRegion2D(double x1, double x2, 
						             double y1, double y2)
{
  Range* therange = new Range[2];
  therange[0] = new double[2];
  therange[1] = new double[2];

  therange[0][0] = x1;
  therange[0][1] = x2;
  therange[1][0] = y1;
  therange[1][1] = y2;
  vector<int> pointsInRange = mytree->get_points_in_range(therange);
  int size = pointsInRange.size();
  for(int i=0;i<size;i++){
	  pointsInRange[i]++;
	  }
  return pointsInRange;
}



const vector<int>&  
HTSTRCInterface:: getPointsInRegion3D(double x1, double x2, 
						                         double y1, double y2,
						                         double z1, double z2)
	{
	Range* therange = new Range[3];
  therange[0] = new double[2];
  therange[1] = new double[2];
  therange[2] = new double[2];

  therange[0][0] = x1;
  therange[0][1] = x2;
  therange[1][0] = y1;
  therange[1][1] = y2;
  therange[1][0] = z1;
  therange[1][1] = z2;

  vector<int> pointsInRange = mytree->get_points_in_range(therange);

  int size = pointsInRange.size();
  for(int i=0;i<size;i++){
	  pointsInRange[i]++;
	  }
  return pointsInRange;

	}

int 
HTSTRCInterface::getClosestPoint2D(double x, double y)
{
	double* point = new double[2];
	point[0] = x;
	point[1] = y;
	int pntIndex;
	int res = mytree->closest_point(point, pntIndex);
	//cout << "the index of the closest node is: " << pntIndex;  
	return pntIndex+1;
}


int 
HTSTRCInterface::getClosestPoint3D(double x, double y, double z)
{
	double* point = new double[3];
	point[0] = x;
	point[1] = y;
	point[2] = z;
	int pntIndex;
	int res = mytree->closest_point(point, pntIndex);
	//cout << "the index of the closest node is: " << pntIndex;  
	return pntIndex+1;
}


//save the temperature data to a file for a specified region
void 
HTSTRCInterface::getTemperature2DRegion(double x1, double x2, double y1, double y2)
{
   
   std::string outstr = "TemperatureFor2DRegion";
   std::stringstream f;
   ++regionCnt2d;
   f << outstr << regionCnt2d;
   std::string filename;
   filename = f.str();
   const char* fname = filename.c_str();

   std::ofstream output;
   output.open(fname, ios::out | ios::app);
   if (!output.is_open()) {
	   cerr << "HTSTRCInterface::getTemperature2DRegion(x1,x2,y1,y2) "
		   << "- could not open the file for output" << endl;
	   exit(-1);
	   }

 //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[2];

  vector<double> thetemp;
  thetemp = this->getFiberTemperature(x1,x2,y1,y2);
  for(int k=0;k<numsteps;k++){
	  fiberTemperatures[k][0] = thetemp[k+numsteps]; //store time
	  fiberTemperatures[k][1] = thetemp[k]; //store temperature
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<2;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output.close();

}




void 
HTSTRCInterface::getTemperature3DRegion(double x1, double x2, double y1, double y2,
	                                    double z1, double z2)
{
   std::string outstr = "TemperatureFor3DRegion";
   std::stringstream f;
   ++regionCnt3d;
   f << outstr << regionCnt3d;
   std::string filename;
   filename = f.str();
   const char* fname = filename.c_str();

   std::ofstream output;
   output.open(fname, ios::out | ios::app);
   if (!output.is_open()) {
	   cerr << "HTSTRCInterface::getTemperature3DRegion(x1,x2,y1,y2,z1,z2) "
		   << "- could not open the file for output" << endl;
	   exit(-1);
	   }

 //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[2];

  vector<double> thetemp;
  thetemp = this->getFiberTemperature(x1,x2,y1,y2,z1,z2);
  for(int k=0;k<numsteps;k++){
	  fiberTemperatures[k][0] = thetemp[k+numsteps]; //store time
	  fiberTemperatures[k][1] = thetemp[k]; //store temperature
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<2;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output.close();	
	
}




void 
HTSTRCInterface::getTemperatureForCompositeSection(double h, double t, double s, double b,
	                                               double slab_w, double slab_d)
{
  std::ofstream output;
  output.open("TemperatureForCompositeSection.dat");
  output << "Time" << "  "; 
  output << "LocBeam1(" << b/2.0 << "," << t/2.0 << ")" << " ";

  double dy = (h-2.0*t)/7.0;
  for(int i=0;i<7;i++){
	output << "LocBeam" << i+2 << "(" << b/2.0 << "," 
		<< t+ dy*(i+0.5) << ")" << " ";
	  }

  output << "LocBeam9(" << b/2 << "," << h-t/2.0 << ")" << " ";

  double dy1 = slab_d/20.0;
  for(int i=0;i<20;i++){
	output << "LocSlab" << i+1 << "(" << b/2.0 << "," 
		<< h + dy1*(i+0.5) << ")" << " ";
	  }
  output << endl;

  //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[30];

  vector<double> thetemp;
  for(int i=0;i<29;i++){
	  if(i==0){
		  thetemp = this->getFiberTemperature(0.0,b,0,t);
		  for(int k=0;k<numsteps;k++){
			  fiberTemperatures[k][i] = thetemp[k+numsteps]; //store time
			  fiberTemperatures[k][i+1] = thetemp[k]; //store temperature
				 }

		  } else if((i>0)&&(i<8)){
			  double xx1 = 0.5*(b-s);
			  double xx2 = 0.5*(b+s);
			  double yy1 = t+dy*i;
			  double yy2 = t+dy*(i+1);
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),t+dy*(i-1),t+dy*i);
			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = thetemp[k];
				  }

		  }else if(i==8) {
			  thetemp = this->getFiberTemperature(0.0,b,h-t,h);
			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = thetemp[k];
				  }

			  } else if(i>8){
				  thetemp = this->getFiberTemperature(0.5*(b-slab_w),0.5*(b+slab_w),h+dy1*(i-9),h+dy1*(i-8));
				  for(int k=0;k<numsteps;k++){
					  fiberTemperatures[k][i+1] = thetemp[k];
					  }
				 }
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<30;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output.close();
}


void 
HTSTRCInterface::getTemperatureForISection(double h, double t, double s, double b)
{
  std::ofstream output;
  output.open("TemperatureForISection.dat");
  output << "Time" << "  "; 
  output << "LocBeam1(" << b/2.0 << "," << t/2.0 << ")" << " ";

  double dy = (h-2.0*t)/7.0;
  for(int i=0;i<7;i++){
	output << "LocBeam" << i+2 << "(" << b/2.0 << "," 
		<< t+ dy*(i+0.5) << ")" << " ";
	  }

  output << "LocBeam9(" << b/2 << "," << h-t/2.0 << ")" << " ";
  output << endl;

  //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[10];

  vector<double> thetemp;
  for(int i=0;i<10;i++){
	  if(i==0){
		  thetemp = this->getFiberTemperature(0,b,0,t);
		  for(int k=0;k<numsteps;k++){
			  fiberTemperatures[k][i] = thetemp[k+numsteps]; //store time
			  fiberTemperatures[k][i+1] = thetemp[k]; //store temperature
				 }

		  } else if((i>0)&&(i<8)){
			  double xx1 = 0.5*(b-s);
			  double xx2 = 0.5*(b+s);
			  double yy1 = t+dy*i;
			  double yy2 = t+dy*(i+1);
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),t+dy*(i-1),t+dy*i);
			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = thetemp[k];
				  }

		  }else if(i==8) {
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),h-t,h);
			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = thetemp[k];
				  }
			  } 
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<10;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output.close();
}


void 
HTSTRCInterface::getTemperatureSeriesForCompositeSection(double h, double t, double s, double b,
	                                               double slab_w, double slab_d)
{
  std::ofstream output;
  output.open("TemperatureSeriesForCompositeSection.dat");
  output << "Time" << "  "; 
  //output << "LocBeam1(" << b/2.0 << "," << t/2.0 << ")" << " ";
  output << "LocBeam1(" << b/2.0 << "," << "0.0" << ")" << " ";

  double dy = (h-2.0*t)/7.0;
  for(int i=0;i<7;i++){
	output << "LocBeam" << i+2 << "(" << b/2.0 << "," 
		<< t+ dy*(i+0.5) << ")" << " ";
	  }

  //output << "LocBeam9(" << b/2 << "," << h-t/2.0 << ")" << " ";
  output << "LocBeam9(" << b/2 << "," << h << ")" << " ";

  double dy1 = slab_d/20.0;
  for(int i=0;i<20;i++){
	output << "LocSlab" << i+1 << "(" << b/2.0 << "," 
		<< h + dy1*(i+0.5) << ")" << " ";
	  }
  output << endl;

  //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[30];

  vector<double> thetemp;
  double dTmax[29];

  for(int i=0;i<29;i++){
	  if(i==0){
		  thetemp = this->getFiberTemperature(0.0,b,0,t);
		  // find the highest temperature in the series
		  dTmax[0] = *(max_element(thetemp.begin(),thetemp.begin()+numsteps))-20.0;

		  for(int k=0;k<numsteps;k++){
			  fiberTemperatures[k][i] = thetemp[k+numsteps]; //store time
			  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[0]; //store temperature
				 }

		  } else if((i>0)&&(i<8)){
			  double xx1 = 0.5*(b-s);
			  double xx2 = 0.5*(b+s);
			  double yy1 = t+dy*i;
			  double yy2 = t+dy*(i+1);
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),t+dy*(i-1),t+dy*i);
			  dTmax[i] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[i];
				  }

		  }else if(i==8) {
			  thetemp = this->getFiberTemperature(0.0,b,h-t,h);
			  dTmax[8] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[8];
				  }

			  } else if(i>8){
				  thetemp = this->getFiberTemperature(0.5*(b-slab_w),0.5*(b+slab_w),h+dy1*(i-9),h+dy1*(i-8));
				  dTmax[i] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

				  for(int k=0;k<numsteps;k++){
					  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[i];
					  }
				 }
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<30;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output << endl;
  output << "dTmax" << " ";
  for(int i=0;i<29;i++){
	  output << dTmax[i]<< " ";
	  }

  output.close();
}


void 
HTSTRCInterface::getTemperatureSeriesForISection(double h, double t, double s, double b)
{
  std::ofstream output;
  output.open("TemperatureSeriesForISection.dat");
  output << "Time" << "  "; 
  //output << "LocBeam1(" << b/2.0 << "," << t/2.0 << ")" << " ";
  output << "LocBeam1(" << b/2.0 << "," << "0.0" << ")" << " ";

  double dy = (h-2.0*t)/7.0;
  for(int i=0;i<7;i++){
	output << "LocBeam" << i+2 << "(" << b/2.0 << "," 
		<< t+ dy*(i+0.5) << ")" << " ";
	  }

  //output << "LocBeam9(" << b/2 << "," << h-t/2.0 << ")" << " ";
  output << "LocBeam9(" << b/2 << "," << h << ")" << " ";
  output << endl;

  //double mytemp;
  double** fiberTemperatures;
  fiberTemperatures = new double*[numsteps];
  for(int i=0;i<numsteps;i++)
	  fiberTemperatures[i] = new double[10];

  vector<double> thetemp;
  double dTmax[9];


  for(int i=0;i<10;i++){
	  if(i==0){
		  thetemp = this->getFiberTemperature(0,b,0,t);
		  dTmax[0] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

		  for(int k=0;k<numsteps;k++){
			  fiberTemperatures[k][i] = thetemp[k+numsteps]; //store time
			  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[0]; //store temperature
				 }

		  } else if((i>0)&&(i<8)){
			  double xx1 = 0.5*(b-s);
			  double xx2 = 0.5*(b+s);
			  double yy1 = t+dy*i;
			  double yy2 = t+dy*(i+1);
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),t+dy*(i-1),t+dy*i);
			  dTmax[i] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[i];
				  }

		  }else if(i==8) {
			  thetemp = this->getFiberTemperature(0.5*(b-s),0.5*(b+s),h-t,h);
			  dTmax[8] = *max_element(thetemp.begin(),thetemp.begin()+numsteps)-20.0;

			  for(int k=0;k<numsteps;k++){
				  fiberTemperatures[k][i+1] = (thetemp[k]-20.0) / dTmax[8];
				  }
			  } 
	  }

  //dump the data to a file
  for(int i=0;i<numsteps;i++){
	  for(int j=0;j<10;j++){
		  output << fiberTemperatures[i][j] << " ";
		  }
	  output << endl;
	  }

  output << endl;
  output << "dTmax" << " ";
  for(int i=0;i<10;i++){
	  output << dTmax[i]<< " ";
	  }

  output.close();
}


void 
HTSTRCInterface::getTemperatureSnapShot(int ntsp, double dtime)
{
   int idx = (ntsp-1)*numnodes;
   double time = ntsp*dtime;

   std::string outstr = "TemperatureSnapShot";
   std::stringstream f;
   f << outstr << time << "s.dat";
   std::string filename;
   filename = f.str();
   const char* fname = filename.c_str();

   std::ofstream outpt;
   outpt.open(fname, ios::out | ios::app);
   if (!outpt.is_open()) {
	   cerr << "HTSTRCInterface::getTemperatureSnapShot "
		   << "- could not open the file for output" << endl;
	   exit(-1);
	   }

   outpt << "$NodeData" << endl;

   int cntx = 0;

   for (int i=0;i<numnodes;i++) {
	   cntx++;
	   int idx1 = idx + i;
	   if  (numdimensions == 2) {
		   //outpt << datastored[idx1][1] << "  " << datastored[idx1][2] << "  " 
			  // << datastored[idx1][3] << endl;
		   outpt << cntx << "  " << datastored[idx1][3] << endl;
		   } else if (numdimensions == 3) {
			   //outpt << datastored[idx1][1] << "  " << datastored[idx1][2] << "  " 
				  // << datastored[idx1][3] << datastored[idx1][4] << endl;
			   outpt << cntx << "  " << datastored[idx1][4] << endl;
		   }
	   }
   outpt << "$EndNodeData";

}
