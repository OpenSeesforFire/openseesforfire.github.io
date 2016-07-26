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


#include <GmshBuilder.h>
#include <HeatTransferDomain.h>
#include <HeatTransferNode.h>
#include <HeatTransferElement.h>
#include <HT_ElementIter.h>
#include <QuadFour.h>
#include <QuadEight.h>
#include <BrickEight.h>
#include <HeatTransferMaterial.h>
#include <NWConcreteEC2.h>
#include <SteelASCE.h>
#include <CarbonSteelEC3.h>
#include <SimpleMaterial.h>
#include <TestMaterial.h>
#include <Radiation.h>
#include <Convection.h>
#include <PrescribedSurfFlux.h>
#include <BoundaryPattern.h>
#include <FireImposedPattern.h>
#include <LocalizedFireEC1.h>
#include <ID.h>
#include <FireModel.h>
#include <ParametricFireEC1.h>
#include <UserDefinedFire.h>
#include <AlpertCeilingJetModel.h>


#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cin;

#include <iomanip>
using std::ios;
#include <stdlib.h>
#include <string>

using std::string;
using std::endl;

//  PlaneFrameModel();
//	constructor
GmshBuilder::GmshBuilder(HeatTransferDomain& theDomain, const char* filename)
:PhysicalID(0),ConnectivityBC(0),ConnectivityDomain(0),
BoundaryEleID(0),numFluxBC(0),numPatterns(0),
TotalElements(0), NumBoundaryEle(0), NumDomainEle(0),
HTModelBuilder(theDomain), fileName(filename),
travellingStamp(false), TravelStop(false),
numNearFieldPattern(0),numFarFieldPattern(0)
{

}

// ~PlaneFrame();    
//	destructor,

GmshBuilder::~GmshBuilder()
{

}    

int
GmshBuilder::BuildModel(void) 
{
    int res = 0;

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:GmshBuilder::buildModel(void) -";
		opserr << " no associated domain \n";
		exit(-1);
		}	

//	char fileName[100];
//	cout << "Enter the name of file containing mesh information: ";
//	cin >> fileName;

	ifstream inputFile(fileName, ios::in);
	if (!inputFile) {
		cout << "FATAL:GmshBuilder::buildModel(void) -";
		cout << " could not open file: " << fileName << endl;
		exit(-1);
		}

	//
	// read in the model paramaters
	// 	  :numNodes numElements numSP_Constraints Num MP_Constraints 
	//     numLoadCases

	int NDIM;
	cout << "Enter the dimension of analysis: " << endl;
	cin  >> NDIM;
	if (NDIM != 2 && NDIM != 3) {
		cout << "Error: GmshBuilder::buildModel(void) -";
		cout <<	"The dimension should be either 2 or 3;" << endl;
		exit(-1);
		}
	cout << endl;


	//string str[30][3];
	//ifstream myfile("mesh.txt");

	//int a = 0;
	//int b = 0;
	//if(!myfile) //Always test the file open.
	//	{
	//	cout<<"Error opening output file"<<endl;
	//	system("pause");
	//	return -1;
	//	}

	inputFile >> NUMPhysicalNames;
	cout << "Number of physical entities: " << NUMPhysicalNames << endl;

	PhysicalID=new int*[NUMPhysicalNames]; //allocate NUMPhysicalNames rows
	for (int j=0;j<NUMPhysicalNames;j++)
		{
		PhysicalID[j]=new int[3]; //allocate 2 columns
		}

	string PhysicalNames;
	int counter = 0;
	string line;

	//while(!inputFile.eof())
	// read in physical entities infomation
	for (int i=0;i<NUMPhysicalNames;i++)
		{
		for (int j=0;j<2;j++)
			{
			inputFile >> PhysicalID[i][j];
			cout << PhysicalID[i][j]<< " ";
			}
		inputFile >> PhysicalNames;
		cout << PhysicalNames << endl;
		}
	cout << endl;

	// read in nodes information
	int NUMNodes;
	inputFile >> NUMNodes;
	cout << "Number of nodes: " << NUMNodes << endl;
	double** NodeCrds;
	NodeCrds = new double*[NUMNodes]; //allocate NUMPhysicalNames rows
	for (int j=0;j<NUMNodes;j++)
		{
		NodeCrds[j]=new double[4]; //allocate 2 columns
		}

	for (int i=0;i<NUMNodes;i++)
		{
		for (int j=0;j<4;j++)
			{
			inputFile >> NodeCrds[i][j];
			//cout << NodeCrds[i][j]<< " ";
			}
		//cout << endl;
		}

	cout << endl;

	// read in elements information
	int BoundaryElementType, DomainElementType;

	inputFile >> TotalElements;
	cout << "Total number of elements: " << TotalElements << endl;
	inputFile >> NumBoundaryEle;
	cout << "Number of elements on boundary: " << NumBoundaryEle << endl;

	inputFile >> BoundaryElementType;
	cout << "Type of elements on boundary: " << BoundaryElementType << endl;                 
	inputFile >> DomainElementType;
	cout << "Type of elements in the domain: " << DomainElementType << endl;

	NumDomainEle = TotalElements - NumBoundaryEle;

	ConnectivityBC = new int*[NumBoundaryEle]; //allocate NumBoundaryEle rows
	ConnectivityDomain = new int*[NumDomainEle]; //allocate NumInnerEle rows

	if (BoundaryElementType == 1){
		for (int j=0;j<NumBoundaryEle;j++)
			{
			ConnectivityBC[j]=new int[7]; //allocate 2 columns
			}
		} else if (BoundaryElementType == 3){ 
			for (int j=0;j<NumBoundaryEle;j++)
				{
				ConnectivityBC[j]=new int[9]; //allocate 2 columns
				}
		} else if (BoundaryElementType == 8){
			for (int j=0;j<NumBoundaryEle;j++)
				{
				ConnectivityBC[j]=new int[8]; //allocate 2 columns
				}
			}

		if (DomainElementType == 3){
			for (int j=0;j<NumDomainEle;j++)
				{
				ConnectivityDomain[j]=new int[9]; //allocate 9 columns
				}
			} else if (DomainElementType == 5 || DomainElementType == 16){ 
				for (int j=0;j<NumDomainEle;j++)
					{
					ConnectivityDomain[j]=new int[13]; //allocate 9 columns
					}
			}

		cout << endl;

		// read in connectivity for elements on boundary
		for (int i=0;i<NumBoundaryEle;i++)
			{
			if (BoundaryElementType == 1){
				for (int j=0;j<7;j++)
					{
					inputFile >> ConnectivityBC[i][j];
					//cout << ConnectivityBC[i][j]<< " ";
					}
				} else if (BoundaryElementType == 3) {
					for (int j=0;j<9;j++)
						{
						inputFile >> ConnectivityBC[i][j];
						//cout << ConnectivityBC[i][j]<< " ";
						}
				}  else if (BoundaryElementType == 8) {
					for (int j=0;j<8;j++)
						{
						inputFile >> ConnectivityBC[i][j];
						//cout << ConnectivityBC[i][j]<< " ";
						}
					}
				//cout << endl;
			}

		//cout << endl;

		// read in connectivity for domain elements
		for (int i=0;i<NumDomainEle;i++)
			{
			if (DomainElementType == 3){
				for (int j=0;j<9;j++)
					{
					inputFile >> ConnectivityDomain[i][j];
					//cout << ConnectivityDomain[i][j]<< " ";
					}
				} else if (DomainElementType == 5 || DomainElementType == 16) {
					for (int j=0;j<13;j++)
						{
						inputFile >> ConnectivityDomain[i][j];
						//cout << ConnectivityDomain[i][j]<< " ";
						}
				}
			//cout << endl;
			}

		// Creating node and element objects and adding them into the domain
		HeatTransferNode* NodePtr;
		bool result;
		int ndof = 1;

		// for each node read in the data, create a node & add it to the domain
		for (int i=0; i<NUMNodes; i++) {
			if (NDIM == 2) {
				NodePtr = new HeatTransferNode(NodeCrds[i][0],ndof,NodeCrds[i][1],
					NodeCrds[i][2]);}
			else if (NDIM == 3) {
				NodePtr = new HeatTransferNode(NodeCrds[i][0],ndof,NodeCrds[i][1],
					NodeCrds[i][2],NodeCrds[i][3]);}
			result = theDomain->addNode(NodePtr);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::BuildModel(void) -";
				opserr << " problems adding node " << i+1 << endln;
				}
			}

		// De-Allocate memory to prevent memory leak
		for (int i = 0; i < NUMNodes; ++i)
			delete[] NodeCrds[i];
		delete [] NodeCrds;

		double T0;
		cout << "Enter initial temperature(in Kelvin): " << endl;
		cin >> T0;
		if (T0 < 273.15) {
			cout << "Warning: GmshBuilder::BuildModel(void) - ";
			cout << "initial temperature is less than 273.15K " << endl;
			}

		theDomain->setInitial(T0); 

		// ask user to specify the material for each physical entity group
		int Mtagcounter = 0;
		HeatTransferMaterial* theMaterialONE;
		HeatTransferMaterial* theMaterialTWO;
		HeatTransferMaterial* theMaterialTHREE;
		HeatTransferMaterial* theMaterialFOUR;
		HeatTransferMaterial* theMaterialFIVE;
		bool myTag1, myTag2, myTag3, myTag4;

		for (int i=0;i<NUMPhysicalNames;i++){
			if (((NDIM == 2) && (PhysicalID[i][0]==2))
				|| ((NDIM) == 3 && (PhysicalID[i][0]==3))) {
					cout << "Enter material tag for domain elements with physical entity tag " << PhysicalID[i][1];
					cout << "( 1 for NWConcreteEC2, 2 for CarbonSteelEC3, 3 for SimpleMaterial, 4 for SteelASCE, 5 for TestMaterial):" << endl;
					cin >> PhysicalID[i][2];

					if (PhysicalID[i][2] != 1 && PhysicalID[i][2] != 2 && 
						PhysicalID[i][2] != 3 && PhysicalID[i][2] != 4 && PhysicalID[i][2] != 5)
						cout << "MaterialTag: " << PhysicalID[i][2] << "not recognized." << endl;

					int phaseTag;

					if (PhysicalID[i][2] == 1){
						//ConcreteEC2(int tag, double moisture);  
						double moist;
						cout << "Enter moisture level for the concrete(only 0, 0.015, 0.03 are supported): " << endl;
						cin >> moist;

						cout << "Phase change?(0 for no, 1 for yes) " << endl;
						cin >> phaseTag;
						if (phaseTag == 0)
							myTag1 = false;
						else
							myTag1 = true;

						theMaterialONE= new NWConcreteEC2(Mtagcounter, moist);
						Mtagcounter++;
						} else if (PhysicalID[i][2] == 2){
							theMaterialTWO = new CarbonSteelEC3(Mtagcounter);
							Mtagcounter++;

							cout << "Phase change?(0 for no, 1 for yes) " << endl;
							cin >> phaseTag;
							if (phaseTag == 0)
								myTag2 = false;
							else
								myTag2 = true;

						} else if (PhysicalID[i][2] == 3){
							// Linear2D(int tag, double rho, double cp, double kc)
							double rho, cp, kc;
							cout << "Enter material density: " << endl;
							cin >> rho;
							cout << "Enter material specific heat: " << endl;
							cin >> cp;
							cout << "Enter material conductivity: " << endl;
							cin >> kc;
							theMaterialTHREE = new SimpleMaterial(Mtagcounter, rho, cp, kc);
							Mtagcounter++;

							//myTag3 = false;
							} else if (PhysicalID[i][2] == 4){
								theMaterialFOUR = new SteelASCE(Mtagcounter);
								Mtagcounter++;
								cout << "Phase change?(0 for no, 1 for yes) " << endl;
								cin >> phaseTag;
								if (phaseTag == 0)
									myTag3 = false;
								else
									myTag3 = true;
							} else if (PhysicalID[i][2] == 5){
								theMaterialFIVE = new TestMaterial(Mtagcounter);
								Mtagcounter++;
								myTag4 = true;
								}
						cout << "Total types of materials used: " << Mtagcounter << endl;
				}
			}

		for (int i=0; i<NumDomainEle; i++) {
			HeatTransferElement* elePtr;
 			if (DomainElementType == 2){
				opserr << "Error: GmshBuilder::BuildModel(void) -";
				opserr << "3-node triangle has not been implemented.";
				} else if (DomainElementType == 3) {
					// creating 4-node quadrangles
					// QuadFour(int tag, int nd1, int nd2, int nd3, int nd4,
					// HeatTransferMaterial& m, bool phaseTransformation = false);

					int physicalEntTag = ConnectivityDomain[i][3];
					int materialchk;
					for (int y=0;y<NUMPhysicalNames;y++){
						if (PhysicalID[y][1] == physicalEntTag){
							materialchk = PhysicalID[y][2];
							break;
							}
						}

					if (materialchk == 1) {
						elePtr = new QuadFour(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
							ConnectivityDomain[i][7], ConnectivityDomain[i][8], *theMaterialONE, myTag1); 
						} else if (materialchk == 2) {
							elePtr = new QuadFour(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
								ConnectivityDomain[i][7], ConnectivityDomain[i][8], *theMaterialTWO, myTag2); 
						} else if (materialchk == 3) {
							elePtr = new QuadFour(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
								ConnectivityDomain[i][7], ConnectivityDomain[i][8], *theMaterialTHREE); 
							} else if (materialchk == 4){
								elePtr = new QuadFour(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][7], ConnectivityDomain[i][8], *theMaterialFOUR,myTag3); 
							} else if (materialchk == 5){
								elePtr = new QuadFour(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][7], ConnectivityDomain[i][8], *theMaterialFIVE,myTag4); 
							}

				} else if (DomainElementType == 5) {
					//	creating 8-node hexahedrons
					//	BrickEight(int tag, int nd1, int nd2, int nd3, int nd4,
					//	int nd5, int nd6, int nd7, int nd8, HeatTransferMaterial& m, 
					//	bool phaseTransformation = false);
					int physicalEntTag = ConnectivityDomain[i][3];
					int materialchk;
					for (int y=0;y<NUMPhysicalNames;y++){
						if (PhysicalID[y][1] == physicalEntTag){
							materialchk = PhysicalID[y][2];
							break;
							}
						}

					if (materialchk == 1) {
						elePtr = new BrickEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][9], 
							ConnectivityDomain[i][10], ConnectivityDomain[i][6], 
							ConnectivityDomain[i][8], ConnectivityDomain[i][12], 
							ConnectivityDomain[i][11], ConnectivityDomain[i][7], 
							*theMaterialONE, myTag1); 
						} else if (materialchk == 2) {
							elePtr = new BrickEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][9], 
								ConnectivityDomain[i][10], ConnectivityDomain[i][6], 
								ConnectivityDomain[i][8], ConnectivityDomain[i][12], 
								ConnectivityDomain[i][11], ConnectivityDomain[i][7], 
								*theMaterialTWO, myTag2); 

						} else if (materialchk == 3) {
							elePtr = new BrickEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][9], 
								ConnectivityDomain[i][10], ConnectivityDomain[i][6], 
								ConnectivityDomain[i][8], ConnectivityDomain[i][12], 
								ConnectivityDomain[i][11], ConnectivityDomain[i][7], 
								*theMaterialTHREE); 
							} else if(materialchk == 4){
								elePtr = new BrickEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][9], 
									ConnectivityDomain[i][10], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][8], ConnectivityDomain[i][12], 
									ConnectivityDomain[i][11], ConnectivityDomain[i][7], 
									*theMaterialFOUR,myTag3); 
							} else if(materialchk == 5){
								elePtr = new BrickEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][9], 
									ConnectivityDomain[i][10], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][8], ConnectivityDomain[i][12], 
									ConnectivityDomain[i][11], ConnectivityDomain[i][7], 
									*theMaterialFIVE,myTag4); 
							}

					} else if (DomainElementType == 16) {
						//	creating 8-node quadrangles
						//	QuadEight(int tag, int nd1, int nd2, int nd3, int nd4,
						//	int nd5, int nd6, int nd7, int nd8, HeatTransferMaterial& m,
						//	bool phaseTransformation = false);
						int physicalEntTag = ConnectivityDomain[i][3];
						int materialchk;
						for (int y=0;y<NUMPhysicalNames;y++){
							if (PhysicalID[y][1] == physicalEntTag){
								materialchk = PhysicalID[y][2];
								break;
								}
							}

						if (materialchk == 1) {
							elePtr = new QuadEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
								ConnectivityDomain[i][7], ConnectivityDomain[i][8], 
								ConnectivityDomain[i][9], ConnectivityDomain[i][10], 
								ConnectivityDomain[i][11], ConnectivityDomain[i][12], 
								*theMaterialONE, myTag1); 
							} else if (materialchk == 2) {
								elePtr = new QuadEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][7], ConnectivityDomain[i][8], 
									ConnectivityDomain[i][9], ConnectivityDomain[i][10], 
									ConnectivityDomain[i][11], ConnectivityDomain[i][12], 
									*theMaterialTWO, myTag2); 
							} else if (materialchk == 3) {
								elePtr = new QuadEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
									ConnectivityDomain[i][7], ConnectivityDomain[i][8], 
									ConnectivityDomain[i][9], ConnectivityDomain[i][10], 
									ConnectivityDomain[i][11], ConnectivityDomain[i][12], 
									*theMaterialTHREE); 
								} else if(materialchk == 4){
									elePtr = new QuadEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
										ConnectivityDomain[i][7], ConnectivityDomain[i][8], 
										ConnectivityDomain[i][9], ConnectivityDomain[i][10], 
										ConnectivityDomain[i][11], ConnectivityDomain[i][12], 
										*theMaterialFOUR,myTag3); 
								} else if(materialchk == 5){
									elePtr = new QuadEight(i+1, ConnectivityDomain[i][5], ConnectivityDomain[i][6], 
										ConnectivityDomain[i][7], ConnectivityDomain[i][8], 
										ConnectivityDomain[i][9], ConnectivityDomain[i][10], 
										ConnectivityDomain[i][11], ConnectivityDomain[i][12], 
										*theMaterialFIVE,myTag4); 
								}
					}

				result = theDomain->addElement(elePtr);
				if (result == false) {
					res =-1;
					opserr << "GmshBuilder::BuildModel(void) -";
					opserr << " problems adding element " << i+1 << endln;
					}
			}

        //************************//
	    //  change begins here!!  //
        //*********************** //

		BoundaryEleID = new int*[NumBoundaryEle]; //allocate NumBoundaryEle rows
		for (int j=0;j<NumBoundaryEle;j++)
			{
			BoundaryEleID[j]=new int[3]; //first column stores domain ele tag, 
			                             //second column stores face tag on boundary,
			                             //third column stores the physical entity tag
			}

		// now loop over boundary elements
		for (int i=0;i<NumBoundaryEle;i++)
			{
			if (BoundaryElementType == 1){ //2-noded lines
              if (DomainElementType == 3){ //4-noded quads
   				  for (int j=0;j<NumDomainEle;j++){
					  // compare the nodes on boundary and find the element with 
					  // with boundary surfaces
					  if ((  (ConnectivityBC[i][5] == ConnectivityDomain[j][5])
						  || (ConnectivityBC[i][5] == ConnectivityDomain[j][6])
						  || (ConnectivityBC[i][5] == ConnectivityDomain[j][7])
						  || (ConnectivityBC[i][5] == ConnectivityDomain[j][8]))
						  && 
						  (  (ConnectivityBC[i][6] == ConnectivityDomain[j][5])
						  || (ConnectivityBC[i][6] == ConnectivityDomain[j][6])
						  || (ConnectivityBC[i][6] == ConnectivityDomain[j][7])
						  || (ConnectivityBC[i][6] == ConnectivityDomain[j][8]))){
							  //store the element tag
							  BoundaryEleID[i][0]= j+1;
							  break;
						  }
					  }

					  // store the face number
					  ID nodesOnFace;
					  HeatTransferElement* elePtr = theDomain->getElement(BoundaryEleID[i][0]);
					  for (int k=1;k<=4;k++){
						  nodesOnFace = elePtr->getNodesOnFace(k);
						  if (( (ConnectivityBC[i][5] == nodesOnFace(0))
							  || (ConnectivityBC[i][5] == nodesOnFace(1)))
							  && ((ConnectivityBC[i][6] == nodesOnFace(0))
							  || (ConnectivityBC[i][6] == nodesOnFace(1)))){
								  BoundaryEleID[i][1]=k;
								  break;
							  }
						  }

					  //store physical entity
					  BoundaryEleID[i][2]= ConnectivityBC[i][3]; 
				  }

				} else if (BoundaryElementType == 3){
					if (DomainElementType == 5){ //8-noded bricks
						for (int j=0;j<NumDomainEle;j++){
							// compare the nodes on boundary and find the element with 
							// with boundary surfaces
							if ((  (ConnectivityBC[i][5] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][12]))
						        && 
								(  (ConnectivityBC[i][6] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][12]))
								&& 
								(  (ConnectivityBC[i][7] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][12]))
								&&
								(  (ConnectivityBC[i][8] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][8] == ConnectivityDomain[j][12]))){
									//store the element tag
									BoundaryEleID[i][0]= j+1;
									break;
								}
							}

							  // store the face number
							  ID nodesOnFace;
							  HeatTransferElement* elePtr = theDomain->getElement(BoundaryEleID[i][0]);
							  for (int k=1;k<=6;k++){
								  nodesOnFace = elePtr->getNodesOnFace(k);
								  if ((  (ConnectivityBC[i][5] == nodesOnFace(0))
									  || (ConnectivityBC[i][5] == nodesOnFace(1))
									  || (ConnectivityBC[i][5] == nodesOnFace(2))
									  || (ConnectivityBC[i][5] == nodesOnFace(3)))
									  && ( (ConnectivityBC[i][6] == nodesOnFace(0))
									  || (ConnectivityBC[i][6] == nodesOnFace(1))
									  || (ConnectivityBC[i][6] == nodesOnFace(2))
									  || (ConnectivityBC[i][6] == nodesOnFace(3)))
									  && ( (ConnectivityBC[i][7] == nodesOnFace(0))
									  || (ConnectivityBC[i][7] == nodesOnFace(1))
									  || (ConnectivityBC[i][7] == nodesOnFace(2))
									  || (ConnectivityBC[i][7] == nodesOnFace(3)))
									  && ( (ConnectivityBC[i][8] == nodesOnFace(0))
									  || (ConnectivityBC[i][8] == nodesOnFace(1))
									  || (ConnectivityBC[i][8] == nodesOnFace(2))
									  || (ConnectivityBC[i][8] == nodesOnFace(3)))){
										  BoundaryEleID[i][1]=k;
										  break;
									  }
								  }

							  //store physical entity
							  BoundaryEleID[i][2]= ConnectivityBC[i][3]; 
						}

				} else if (BoundaryElementType == 8){
					if (DomainElementType == 16){ //8-noded quads
						for (int j=0;j<NumDomainEle;j++){
							// compare the nodes on boundary and find the element with 
							// with boundary surfaces
							if ((  (ConnectivityBC[i][5] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][5] == ConnectivityDomain[j][12]))
						        && 
								(  (ConnectivityBC[i][6] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][6] == ConnectivityDomain[j][12]))
								&& 
								(  (ConnectivityBC[i][7] == ConnectivityDomain[j][5])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][6])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][7])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][8])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][9])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][10])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][11])
								|| (ConnectivityBC[i][7] == ConnectivityDomain[j][12]))){
									//store the element tag
									BoundaryEleID[i][0]= j+1;
									break;
								}
							}

							  // store the face number
							  ID nodesOnFace;
							  HeatTransferElement* elePtr = theDomain->getElement(BoundaryEleID[i][0]);
							  for (int k=1;k<=4;k++){
								  nodesOnFace = elePtr->getNodesOnFace(k);
								  if ((  (ConnectivityBC[i][5] == nodesOnFace(0))
									  || (ConnectivityBC[i][5] == nodesOnFace(1))
									  || (ConnectivityBC[i][5] == nodesOnFace(2)))
									  && ( (ConnectivityBC[i][6] == nodesOnFace(0))
									  || (ConnectivityBC[i][6] == nodesOnFace(1))
									  || (ConnectivityBC[i][6] == nodesOnFace(2)))
									  && ( (ConnectivityBC[i][7] == nodesOnFace(0))
									  || (ConnectivityBC[i][7] == nodesOnFace(1))
									  || (ConnectivityBC[i][7] == nodesOnFace(2)))){
										  BoundaryEleID[i][1]=k;
										  break;
									  }
								  }

							  //store physical entity
							  BoundaryEleID[i][2]= ConnectivityBC[i][3]; 
						}
					}
					}

    // De-Allocate memory to prevent memory leak
	for (int i = 0; i < NumBoundaryEle; ++i)
		delete[] ConnectivityBC[i];
	delete [] ConnectivityBC;

	for (int i = 0; i < NumDomainEle; ++i)
		delete[] ConnectivityDomain[i];
	delete [] ConnectivityDomain;


	return res;
	}


int 
GmshBuilder::setConvectionBC(double h, double Tf)
{
    bool result;
	int res = 0;
	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:setConvectionBC(double h, double Tf) - ";
		opserr << "no associated domain \n";
		exit(-1);
		}	

	int numPatterns = this->createBoundaryPattern();

	Convection* theConvection;
	for(int i=0;i<NumBoundaryEle;i++){
		//Convection(int tag, int eleTag, int fag, double h, double T);
		int eleTag = BoundaryEleID[i][0];
		int fTag = BoundaryEleID[i][1];
		numFluxBC++;
		theConvection = new Convection(numFluxBC, eleTag, fTag, h, Tf);
		result = theDomain->addHeatFluxBC(theConvection,numPatterns);
		if (result == false) {
			res =-1;
			opserr << "GmshBuilder::setConvectionBC(double h, double Tf) -";
			opserr << " problems adding ConvectionBC at element " << eleTag << endln;
			}
		}

	return res;
}


int 
GmshBuilder::setConvectionBC(double h, double Tf, int tag)
{
    int res = 0;
	bool result;
	int checker = 0;
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL:setConvectionBC(double h, double Tf, int tag) - ";
		opserr << " no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:GmshBuilder::setConvectionBC(double h, double Tf, int tag) -";
		opserr << " no associated domain \n";
		exit(-1);
		}	

	int numPatterns = this->createBoundaryPattern();

	Convection* theConvection;
	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {

			//Convection(int tag, int eleTag, int fag, double h, double T);
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			theConvection = new Convection(numFluxBC, eleTag, fTag, h, Tf);
			result = theDomain->addHeatFluxBC(theConvection,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setConvectionBC(double h, double Tf, int tag) -";
				opserr << " problems adding ConvectionBC at element " << eleTag << endln;
				}
			}
		}

	return res;
}
	    

int 
GmshBuilder::setRadiationBC(double epsilon, double sigma, double alpha, double qir)
{
    bool result;
	int res = 0;
	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:setRadiationBC(double epsilon, double sigma, double alpha, double qir) -";
		opserr << " no associated domain \n";
		exit(-1);
		}	

	int numPatterns = this->createBoundaryPattern();

	Radiation* theRadiation;
	for(int i=0;i<NumBoundaryEle;i++){
		//Convection(int tag, int eleTag, int fag, double h, double T);
		int eleTag = BoundaryEleID[i][0];
		int fTag = BoundaryEleID[i][1];
		numFluxBC++;
		theRadiation = new Radiation(numFluxBC, eleTag, fTag, epsilon, 
			                          sigma, alpha, qir);
		result = theDomain->addHeatFluxBC(theRadiation,numPatterns);
		if (result == false) {
			res =-1;
			opserr << "GmshBuilder::setRadiationBC(double epsilon, double sigma, double alpha, double qir) -";
			opserr << " problems adding RadiationBC at element " << eleTag << endln;
			}
		}
	return res;
}


int 
GmshBuilder::setRadiationBC(double epsilon, double sigma, double alpha, double qir, int tag)
{
	int checker = 0;
	bool result;
	int res = 0;

	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL:GmshBuilder::setRadiationBC(double epsilon, double sigma, double alpha, double qir, int tag) - ";
		opserr << " no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:setRadiationBC(double epsilon, double sigma, double alpha, double qir, int tag) -";
		opserr << " no associated domain \n";
		exit(-1);
		}	

	int numPatterns = this->createBoundaryPattern();

	Radiation* theRadiation;
	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			theRadiation = new Radiation(numFluxBC, eleTag, fTag, epsilon, 
				sigma, alpha, qir);
			result = theDomain->addHeatFluxBC(theRadiation,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setRadiationBC(double epsilon, double sigma, double alpha, double qir, int tag) -";
				opserr << " problems adding RadiationBC at element " << eleTag << endln;
				}
			}
		}

	return res;
	}


int 
GmshBuilder::setDirichletBC(double T, int tag)
{
    int checker = 0;
	bool result;
	int res = 0;
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setDirichletBC(double T, int tag) - ";
		opserr << " no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL: GmshBuilder::setDirichletBC(double T, int tag) -";
		opserr << " no associated domain \n";
		exit(-1);
		}	

	TemperatureBC* theTemperatureBC;
	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			HeatTransferElement* theEle = theDomain->getElement(eleTag);
			const ID& nodes = theEle->getNodesOnFace(fTag);
			int NumFnodes = nodes.Size();
			for(int i=0;i<NumFnodes;i++) {
				theTemperatureBC = new TemperatureBC(nodes(i), 0, T);
				result = theDomain->addTemperatureBC(theTemperatureBC);	
				if (result == false) {
					res =-1;
					opserr << "GmshBuilder::setDirichletBC(double T, int tag) -";
					opserr << " problems adding RadiationBC at node " << nodes(i) << endln;
					}
				}
			}
		}

	return res;
	}


int
GmshBuilder::setParametricFireConvcBC(ParametricFireEC1*  themodel, double h, int tag)
{
    int checker = 0;
	bool result;
	int res = 0;	
	
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setParametricFireBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setParametricFireBC() - no associated domain \n";
		exit(-1);
		}	

	ParametricFireEC1* theFireModel = themodel;
	//theFireModel = themodel;
	int numPatterns = this->createFireImposedPattern(themodel);

	Convection* theConvecBC = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			// Convection(int tag, int eleTag, int fag, double h, double T);
			double Tf = 293.15;
			theConvecBC = new Convection(numFluxBC, eleTag, fTag, h, Tf);
			result = theDomain->addHeatFluxBC(theConvecBC,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setParametricFireConvcBC() -";
				opserr << " problems adding Convection at element " << eleTag << endln;
				}
			}
		}

	return res;
}


int
GmshBuilder::setParametricFireRadBC(ParametricFireEC1*  themodel, double epsilon, 
									double sigma, double alpha,
									int tag)
{
	int checker = 0;
	bool result;
	int res = 0;	
	
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setParametricFireRadBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setParametricFireRadBC() - no associated domain \n";
		exit(-1);
		}	

	ParametricFireEC1* theFireModel = themodel;
	int numPatterns = this->createFireImposedPattern(themodel);

	Radiation* theRadiationBC = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
            //Radiation(int tag, int eleTag, int fag, double epsilon, 
			//          double sigma, double alpha, double qir);
			double qir = 0.0;
			theRadiationBC = new Radiation(numFluxBC, eleTag, fTag, epsilon, sigma,alpha,qir);
			result = theDomain->addHeatFluxBC(theRadiationBC,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setParametricFireRadBC() -";
				opserr << " problems adding Radiation at element " << eleTag << endln;
				}
			}
		}

	return res;
}

//int
//GmshBuilder::setParametricFireConvcBC(double I, double Av, double H, double At,
//		                              double Af, double Qf, double tlim, 
//									  double h, int tag)
//{
//	int checker = 0;
//	bool result;
//	int res = 0;	
//	
//	for(int i=0;i<NUMPhysicalNames;i++){
//		if (tag == PhysicalID[i][1])
//			checker++;
//		}
//
//	if (checker == 0) {
//		opserr << "FATAL: GmshBuilder::setParametricFireBC() - no tag found in the input mesh file \n";
//		exit(-1);
//		}
//
//	HeatTransferDomain* theDomain = this->getDomainPtr();    
//	if (theDomain == 0) {
//		opserr << "FATAL:FATAL: GmshBuilder::setParametricFireBC() - no associated domain \n";
//		exit(-1);
//		}	
//
//	numPatterns++;
//	FireModel* theFireModel = new  ParametricFireEC1(I, Av, H, At, Af, Qf, tlim);
//	FireImposedPattern* theFirePattern = new FireImposedPattern(numPatterns);
//    theFirePattern->setFireModel(theFireModel);
//	theDomain->addBoundaryPattern(theFirePattern);	
//
//	Convection* theConvecBC = 0;
//
//	for(int i=0;i<NumBoundaryEle;i++){
//		if (BoundaryEleID[i][2] == tag) {
//			int eleTag = BoundaryEleID[i][0];
//			int fTag = BoundaryEleID[i][1];
//			numFluxBC++;
//			// Convection(int tag, int eleTag, int fag, double h, double T);
//			double Tf = 293.15;
//			theConvecBC = new Convection(numFluxBC, eleTag, fTag, h, Tf);
//			result = theDomain->addHeatFluxBC(theConvecBC,numPatterns);
//			if (result == false) {
//				res =-1;
//				opserr << "GmshBuilder::setParametricFireConvcBC() -";
//				opserr << " problems adding Convection at element " << eleTag << endln;
//				}
//			}
//		}
//
//	return res;
//}
//
//
//int
//GmshBuilder::setParametricFireRadBC(double I, double Av, double H, double At,
//		                            double Af, double Qf, double tlim, 
//									double epsilon, double sigma, double alpha,
//									int tag)
//{
//	int checker = 0;
//	bool result;
//	int res = 0;	
//	
//	for(int i=0;i<NUMPhysicalNames;i++){
//		if (tag == PhysicalID[i][1])
//			checker++;
//		}
//
//	if (checker == 0) {
//		opserr << "FATAL: GmshBuilder::setParametricFireRadBC() - no tag found in the input mesh file \n";
//		exit(-1);
//		}
//
//	HeatTransferDomain* theDomain = this->getDomainPtr();    
//	if (theDomain == 0) {
//		opserr << "FATAL:FATAL: GmshBuilder::setParametricFireRadBC() - no associated domain \n";
//		exit(-1);
//		}	
//
//	numPatterns++;
//	FireModel* theFireModel = new  ParametricFireEC1(I, Av, H, At, Af, Qf, tlim);
//	FireImposedPattern* theFirePattern = new FireImposedPattern(numPatterns);
//    theFirePattern->setFireModel(theFireModel);
//	theDomain->addBoundaryPattern(theFirePattern);	
//
//	Radiation* theRadiationBC = 0;
//
//	for(int i=0;i<NumBoundaryEle;i++){
//		if (BoundaryEleID[i][2] == tag) {
//			int eleTag = BoundaryEleID[i][0];
//			int fTag = BoundaryEleID[i][1];
//			numFluxBC++;
//            //Radiation(int tag, int eleTag, int fag, double epsilon, 
//			//          double sigma, double alpha, double qir);
//			double qir = 0.0;
//			theRadiationBC = new Radiation(numFluxBC, eleTag, fTag, epsilon, sigma,alpha,qir);
//			result = theDomain->addHeatFluxBC(theRadiationBC,numPatterns);
//			if (result == false) {
//				res =-1;
//				opserr << "GmshBuilder::setParametricFireRadBC() -";
//				opserr << " problems adding Radiation at element " << eleTag << endln;
//				}
//			}
//		}
//
//	return res;
//}


int 
GmshBuilder::setLocalisedFireBC(double crd1, double crd2, double crd3, 
		                        double D, double Q, double H, 
						        int centerLineTag, int tag)
{
	int checker = 0;
	bool result;
	int res = 0;

	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setLocalisedFireBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setLocalisedFireBC() - no associated domain \n";
		exit(-1);
		}	

	LocalizedFireEC1* theFireModel = new LocalizedFireEC1(crd1, crd2, crd3, D, Q, H, centerLineTag);
	int numPatterns = this->createFireImposedPattern(theFireModel);

	PrescribedSurfFlux* thePflux = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			int eleType = ConnectivityBC[i][1];
			//identify number of nodes on face
			int numNPF;
			if(eleType == 1){
				numNPF = 2;
				} else if (eleType == 2) {
					numNPF = 3;
				} else if (eleType == 3) {
					numNPF = 4;
					} else if (eleType == 8) {
						numNPF = 3;
					} else if (eleType == 16) { 
						numNPF = 8;
						}
					thePflux = new PrescribedSurfFlux(numFluxBC, eleTag, fTag, numNPF);
					result = theDomain->addHeatFluxBC(thePflux,numPatterns);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setLocalisedFireBC() -";
						opserr << " problems adding PrescribedSurfFlux at element " << eleTag << endln;
						}
			}
		}

	return res;
	}


int 
GmshBuilder::setUDFFireConvcBC(UserDefinedFire* themodel, double h, int tag)
{
   int checker = 0;
   bool result;
   int res = 0;	
   for(int i=0;i<NUMPhysicalNames;i++){
	   if (tag == PhysicalID[i][1])
		   checker++;
	   }

   if (checker == 0) {
	   opserr << "FATAL: GmshBuilder::setUDFFireConvcBC() - no tag found in the input mesh file \n";
	   exit(-1);
	   }

   HeatTransferDomain* theDomain = this->getDomainPtr();    
   if (theDomain == 0) {
	   opserr << "FATAL:FATAL: GmshBuilder::setUDFFireConvcBC() - no associated domain \n";
	   exit(-1);
	   }	

   UserDefinedFire* theFireModel = new  UserDefinedFire();
   theFireModel = themodel;
   int numPatterns = this->createFireImposedPattern(themodel);

   Convection* theConvecBC = 0;

   for(int i=0;i<NumBoundaryEle;i++){
	   if (BoundaryEleID[i][2] == tag) {
		   int eleTag = BoundaryEleID[i][0];
		   int fTag = BoundaryEleID[i][1];
		   numFluxBC++;
		   // note: Tf here is a pseudo temperature
		   double Tf = 293.15;
		   theConvecBC = new Convection(numFluxBC, eleTag, fTag, h, Tf);
		   result = theDomain->addHeatFluxBC(theConvecBC,numPatterns);
		   if (result == false) {
			   res =-1;
			   opserr << "GmshBuilder::setUDFFireConvcBC() -";
			   opserr << " problems adding Convection at element " << eleTag << endln;
			   }
		   }
	   }

   return res;	
	}


int 
GmshBuilder::setUDFFireRadBC(UserDefinedFire* themodel, double epsilon, 
		                double sigma, double alpha,
						int tag)
{
   int checker = 0;
   bool result;
   int res = 0;	
   for(int i=0;i<NUMPhysicalNames;i++){
	   if (tag == PhysicalID[i][1])
		   checker++;
	   }

   if (checker == 0) {
	   opserr << "FATAL: GmshBuilder::setUDFFireRadBC() - no tag found in the input mesh file \n";
	   exit(-1);
	   }

   HeatTransferDomain* theDomain = this->getDomainPtr();    
   if (theDomain == 0) {
	   opserr << "FATAL:FATAL: GmshBuilder::setUDFFireRadBC() - no associated domain \n";
	   exit(-1);
	   }	

   UserDefinedFire* theFireModel = themodel;
   int numPatterns = this->createFireImposedPattern(themodel);


   Radiation* theRadiationBC = 0;

   for(int i=0;i<NumBoundaryEle;i++){
	   if (BoundaryEleID[i][2] == tag) {
		   int eleTag = BoundaryEleID[i][0];
		   int fTag = BoundaryEleID[i][1];
		   numFluxBC++;
			double qir = 0.0;
			theRadiationBC = new Radiation(numFluxBC, eleTag, fTag, epsilon, sigma,alpha,qir);
			result = theDomain->addHeatFluxBC(theRadiationBC,numPatterns);
		   if (result == false) {
			   res =-1;
			   opserr << "GmshBuilder::setUDFFireRadBC() -";
			   opserr << " problems adding Radiation at element " << eleTag << endln;
			   }
		   }
	   }

   return res;		
}


int 
GmshBuilder::setUDFFlux(UserDefinedFire* themodel, int tag)
{
    int checker = 0;
	bool result;
	int res = 0;

	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setUDFFlux() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setUDFFlux() - no associated domain \n";
		exit(-1);
		}	

	FireModel* theFireModel = themodel;
    int numPatterns = this->createFireImposedPattern(themodel);

	PrescribedSurfFlux* thePflux = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			int eleType = ConnectivityBC[i][1];
			//identify number of nodes on face
			int numNPF;
			if(eleType == 1){
				numNPF = 2;
				} else if (eleType == 2) {
					numNPF = 3;
				} else if (eleType == 3) {
					numNPF = 4;
					} else if (eleType == 8) {
						numNPF = 3;
					} else if (eleType == 16) { 
						numNPF = 8;
						}
					thePflux = new PrescribedSurfFlux(numFluxBC, eleTag, fTag, numNPF);
					result = theDomain->addHeatFluxBC(thePflux,numPatterns);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setUDFFlux() -";
						opserr << " problems adding PrescribedSurfFlux at element " << eleTag << endln;
						}
			}
		}

	return res;	
	
}


int
GmshBuilder::setAlpertFireConvcBC(AlpertCeilingJetModel* themodel, double h, int tag)
{
    int checker = 0;
	bool result;
	int res = 0;	
	
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setAlpertFireConvcBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setAlpertFireConvcBC() - no associated domain \n";
		exit(-1);
		}	

	AlpertCeilingJetModel* theFireModel = themodel;
    int numPatterns = this->createFireImposedPattern(theFireModel);

	Convection* theConvecBC = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
			// Convection(int tag, int eleTag, int fag, double h, double T);
			double Tf = 293.15;
			theConvecBC = new Convection(numFluxBC, eleTag, fTag, h, Tf);
			result = theDomain->addHeatFluxBC(theConvecBC,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setAlpertFireConvcBC() -";
				opserr << " problems adding Convection at element " << eleTag << endln;
				}
			}
		}

	return res;
}


int
GmshBuilder::setAlpertFireRadBC(AlpertCeilingJetModel* themodel, double epsilon, 
		                        double sigma, double alpha, int tag)
{
	int checker = 0;
	bool result;
	int res = 0;	
	
	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setAlpertFireRadBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setAlpertFireRadBC() - no associated domain \n";
		exit(-1);
		}	

	AlpertCeilingJetModel* theFireModel = themodel;
    int numPatterns = this->createFireImposedPattern(theFireModel);

	Radiation* theRadiationBC = 0;

	for(int i=0;i<NumBoundaryEle;i++){
		if (BoundaryEleID[i][2] == tag) {
			int eleTag = BoundaryEleID[i][0];
			int fTag = BoundaryEleID[i][1];
			numFluxBC++;
            //Radiation(int tag, int eleTag, int fag, double epsilon, 
			//          double sigma, double alpha, double qir);
			double qir = 0.0;
			theRadiationBC = new Radiation(numFluxBC, eleTag, fTag, epsilon, sigma,alpha,qir);
			result = theDomain->addHeatFluxBC(theRadiationBC,numPatterns);
			if (result == false) {
				res =-1;
				opserr << "GmshBuilder::setAlpertFireRadBC() -";
				opserr << " problems adding Radiation at element " << eleTag << endln;
				}
			}
		}

	return res;
}


int
GmshBuilder::setLinearTravellingFireBC(AlpertCeilingJetModel* theAlpModel, double pos1, double pos2, 
		                               double Lmax, double h, double Ta, double Tf, double epsilon, 
									   double sigma, double alpha, int DirectionFlag, int tag)
{
   int checker = 0;
	bool result;
	int res = 0;	

	for(int i=0;i<NUMPhysicalNames;i++){
		if (tag == PhysicalID[i][1])
			checker++;
		}

	if (checker == 0) {
		opserr << "FATAL: GmshBuilder::setLinearTravellingFireBC() - no tag found in the input mesh file \n";
		exit(-1);
		}

	HeatTransferDomain* theDomain = this->getDomainPtr();    
	if (theDomain == 0) {
		opserr << "FATAL:FATAL: GmshBuilder::setLinearTravellingFireBC() - no associated domain \n";
		exit(-1);
		}	


	//count the number of boundary elements exposed to travelling fires
	static int cntx;
	if(travellingStamp == false){
		for(int i=0;i<NumBoundaryEle;i++){
			if (BoundaryEleID[i][2] == tag) {
				cntx++;
				}
			}
		}

   static Convection** ConvectionBCs = new Convection*[cntx];
   static Radiation** RadiationBCs = new Radiation*[cntx];

   // create heatflux objects for surfaces elements exposed to travelling fires
   // note: just need to create them at the begning
   if(travellingStamp == false){
	   int counter = 0;
	   for(int i=0;i<NumBoundaryEle;i++){
		   if (BoundaryEleID[i][2] == tag) {
			   //initialize the pointers and associate fluxbcs with elements
			   int eleTag = BoundaryEleID[i][0];
			   int fTag = BoundaryEleID[i][1];    
			   double T0 = 293.15;
			   double qir = sigma * pow(T0,4.0);
			   numFluxBC++;
			   ConvectionBCs[counter] = new  Convection(numFluxBC, eleTag, fTag, h, T0);
			   numFluxBC++;
			   RadiationBCs[counter] = new Radiation(numFluxBC, eleTag, fTag, epsilon, sigma,alpha,qir);
			   counter++;
			   }
		   }
	   }



	if(pos1 < Lmax){
		
		//if the fire starts travelling, remove the patterns at 
		//previous time step
		if(travellingStamp == true){ 
			this->removePattern(numNearFieldPattern);
			this->removePattern(numFarFieldPattern);

            // remove heatfluxbcs from elements added in previous time steps
			for(int i=0;i<cntx;i++){
				int eleTag = ConvectionBCs[i]->getElementTag();
				HeatTransferElement* theEle = theDomain->getElement(eleTag);
				//remove the fluxbcs added in the last time step
				theEle->removeConvection(ConvectionBCs[i]);
				theEle->removeRadiation(RadiationBCs[i]);
				}

			}


		//create normal patterns for the near field
		int numPatterns1 = this->createBoundaryPattern();
		numNearFieldPattern = numPatterns1;

		//create fire patterns for the far field
		AlpertCeilingJetModel* theFireModel = theAlpModel;
		int numPatterns2 = this->createFireImposedPattern(theFireModel);
		numFarFieldPattern = numPatterns2;

		for(int i=0;i<cntx;i++){
				int eleTag = ConvectionBCs[i]->getElementTag();
				int fTag = ConvectionBCs[i]->getFaceTag();

				// check if the element faces are located within the near field
				HeatTransferElement* theEle = theDomain->getElement(eleTag);
				const ID& faceNodes = theEle->getNodesOnFace(fTag);
				int size = faceNodes.Size();

				int nearFieldChecker = 0;

				if(DirectionFlag == 1){ // if parallel to travelling direction
					
					//check if all the nodes on the face fall in the near-field
					for (int i2 = 0; i2 < size; i2++) {	
						int nodTag = faceNodes(i2);
						HeatTransferNode* theNode = theDomain->getNode(nodTag);
						if (theNode == 0) {
							opserr << "GmshBuilder::setLinearTravellingFireBC() - no node with tag " << nodTag << " exists in the domain";
							exit(-1);
							}
						const Vector& coords = theNode->getCrds();

						if((coords(0)<pos1) || (coords(0)>pos2))
							{ 
							nearFieldChecker++;
							}
						}
					} else if(DirectionFlag == 0){ // if perpendicular to travelling direction

						if(((pos1<0.0) && (pos2<0.0)) || ((pos1>0.0) && (pos2>0.0)))
							{ 
							nearFieldChecker++;
							}
					} else {
						opserr << "GmshBuilder::setLinearTravellingFireBC() - incorrect DirectionFlag " << DirectionFlag << endln;
						}

				// if the face falls in the near field
				if(nearFieldChecker == 0)
					{
					//add Radiation 		
					double qir = sigma * pow(Tf,4.0);
					RadiationBCs[i]->setIrradiation(qir);
					result = theDomain->addHeatFluxBC(RadiationBCs[i],numPatterns1);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setLinearTravellingFireBC() -";
						opserr << " problems adding Radiation at element " << eleTag << endln;
						}

					//add Convection
					ConvectionBCs[i]->setSurroundingTemp(Tf);
					result = theDomain->addHeatFluxBC(ConvectionBCs[i],numPatterns1);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setLinearTravellingFireBC() -";
						opserr << " problems adding Convection at element " << eleTag << endln;
						}
					} else { // faces are in the far field

						result = theDomain->addHeatFluxBC(ConvectionBCs[i],numPatterns2);
						if (result == false) {
							res =-1;
							opserr << "GmshBuilder::setLinearTravellingFireBC() -";
							opserr << " problems adding Convection at element " << eleTag << endln;
							}

						result = theDomain->addHeatFluxBC(RadiationBCs[i],numPatterns2);
						if (result == false) {
							res =-1;
							opserr << "GmshBuilder::setLinearTravellingFireBC() -";
							opserr << " problems adding Radiation at element " << eleTag << endln;
							}
					}
			}
		travellingStamp = true;

		} else {// if fuel burnt out, then subject to ambient temperature
			if(TravelStop==false){
				this->removePattern(numNearFieldPattern);
				this->removePattern(numFarFieldPattern);

				//remove all the previous flux bcs associated with fire exposed elements
				for(int i=0;i<NumBoundaryEle;i++){
					if (BoundaryEleID[i][2] == tag) {
						int eleTag = BoundaryEleID[i][0];

						HeatTransferElement* elePtr = theDomain->getElement(eleTag);
						elePtr->clearAllFluxBCs();  
						}
					}

				int numPatterns3 = this->createBoundaryPattern();

				for(int i=0;i<cntx;i++){
					int eleTag = ConvectionBCs[i]->getElementTag();
						
					// add Radiation 		
					double qir = sigma * pow(Ta,4.0);
					RadiationBCs[i]->setIrradiation(qir);
					result = theDomain->addHeatFluxBC(RadiationBCs[i],numPatterns3);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setLinearTravellingFireBC() -";
						opserr << " problems adding Radiation at element " << eleTag << endln;
						}

					//add Convection
					ConvectionBCs[i]->setSurroundingTemp(Ta);
					result = theDomain->addHeatFluxBC(ConvectionBCs[i],numPatterns3);
					if (result == false) {
						res =-1;
						opserr << "GmshBuilder::setLinearTravellingFireBC() -";
						opserr << " problems adding Convection at element " << eleTag << endln;
						}
					}
				TravelStop = true;
				}
			}
	return res;

	}


int
GmshBuilder::createBoundaryPattern(void)
{
    numPatterns++;
	BoundaryPattern* thePattern = new BoundaryPattern(numPatterns);
	HeatTransferDomain* theDomain = this->getDomainPtr();    
	theDomain->addBoundaryPattern(thePattern);

	return numPatterns;
}


int
GmshBuilder::createFireImposedPattern(FireModel* themodel)
{
    numPatterns++;
	FireImposedPattern* thePattern = new FireImposedPattern(numPatterns);
	thePattern->setFireModel(themodel);
	HeatTransferDomain* theDomain = this->getDomainPtr();    
	theDomain->addBoundaryPattern(thePattern);
	
	return numPatterns;
}

void 
GmshBuilder::removePattern(int tag)
{
    HeatTransferDomain* theDomain = this->getDomainPtr();    
	BoundaryPattern* removedPattern = theDomain->removeBoundaryPattern(tag);
}
