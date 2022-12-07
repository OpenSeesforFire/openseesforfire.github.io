 #include <Simple_Composite3D.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <BrickEight.h>

	/*
Simple_Composite3D::Simple_Composite3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ, 
double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw, double HTI_UBf,double HTI_UTf,double HTI_Len)
:Simple_Entity(tag,3),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_centerZ(HTI_centerZ),
HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf),HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw), HTI_UBf(HTI_UBf), HTI_UTf(HTI_UTf),
HTI_Len(HTI_Len),NumCtrlID(5),Seeds1(200),Seeds2(200),Seeds3(1000)
{
	
}
*/

Simple_Composite3D::Simple_Composite3D(int tag, double HTI_centerX, double HTI_centerY, double HTI_centerZ,
double HTI_Bf, double HTI_Tf, double HTI_Hw, double HTI_Tw,double HTI_Len, double slabW, double slabH)
:Simple_Entity(tag,7),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY),HTI_centerZ(HTI_centerZ), HTI_Bf(HTI_Bf), 
HTI_Tf(HTI_Tf),HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw),HTI_UBf(HTI_Bf), HTI_UTf(HTI_Tf),HTI_Len(HTI_Len),  SlabW(slabW),SlabH(slabH),
NumCtrlID(7),Seeds1(0),Seeds2(0),Seeds3(0),Seeds4(0), Seeds5(0)
{

}


Simple_Composite3D::~Simple_Composite3D()
{


}

int 
Simple_Composite3D::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
	double EleX = MeshCtrls(0);
	double EleY = MeshCtrls(1);
	double EleZ = MeshCtrls(2);
	double EleX_Web = MeshCtrls(3);
	double EleY_Web = MeshCtrls(4);
	double EleX_Slab = MeshCtrls(5);
	double EleY_Slab = MeshCtrls(6);
	
	
	NumCtrlID(0)= (HTI_Bf-HTI_Tw)/EleX+HTI_Tw/EleX_Web+0.5;  //Num of elements generated in flange along x-axis, web is meshed independently
	NumCtrlID(1)= HTI_Tf/EleY+0.5;                           //Num of elements generated in flange along Y-axis
	NumCtrlID(2) = HTI_Tw/EleX_Web+0.5;                      //Num of elements generated in Web along x-axis
	NumCtrlID(3) = HTI_Hw/EleY_Web+0.5;                       //Num of elements generated in Web along y-axis
	NumCtrlID(4) = HTI_Len/EleZ+0.5;                          // An integral appointed number of elments along length
	NumCtrlID(5) = (SlabW-HTI_Bf)/EleX_Slab+(HTI_Bf-HTI_Tw)/EleX+HTI_Tw/EleX_Web+0.5;    //Number of elements generated in slab along local x axis;
	NumCtrlID(6) = SlabH/EleY_Slab;							//Number of elements generated in slab along local y axis;
	
	
	if ((NumCtrlID(0)-NumCtrlID(2))%2!=0)
		{
			opserr<< "Error: Invalid dimension input for Isection " <<this->getTag()<<endln;
			return -1;
		}
	
	if(NumCtrlID(0)==0||NumCtrlID(1)==0||NumCtrlID(2)==0||NumCtrlID(3)==0||NumCtrlID(4)==0||NumCtrlID(5)==0||NumCtrlID(6)==0)
		{
			opserr<<"Error: NumCtrlID in Isection: "<< this->getTag()<<" has 0 vlaue"<<endln;
			return -2;
		}
	return 0;

}


bool 
Simple_Composite3D::InitialSeeds(void)
{
	bool result = true;
	Seeds1.resize(NumCtrlID(0)+1);
	Seeds2.resize(NumCtrlID(1)*2+NumCtrlID(3)+1);
	Seeds3.resize(NumCtrlID(4)+1);
	Seeds4.resize(NumCtrlID(5)+1);
	Seeds5.resize(NumCtrlID(6)+1);

	for(int i =0; i<=NumCtrlID(0); i++){
		
		double a0,incr;
		if (i< (NumCtrlID(0)-NumCtrlID(2))/2)
		{
		    a0 = HTI_centerX-HTI_Bf/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds1(i)= a0+incr*i;
		}
		else if(i< (NumCtrlID(0)+NumCtrlID(2))/2)
		{
			a0 = HTI_centerX-HTI_Tw/2;  //OriginCrdX from the web
		    incr = HTI_Tw/NumCtrlID(2);
			Seeds1(i)= a0+incr*(i-(NumCtrlID(0)-NumCtrlID(2))/2);
		}
		else
		{
			a0 = HTI_centerX+HTI_Tw/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds1(i)= a0+incr*(i-(NumCtrlID(0)+NumCtrlID(2))/2);
		}
	
	}

	for(int i =0; i<= NumCtrlID(1)*2+NumCtrlID(3); i++){
        double b0,incr;
		
		if (i< NumCtrlID(1))
		{
		    b0 = HTI_centerY-HTI_Hw/2-HTI_Tf;  //OriginCrdY
		    incr = HTI_Tf/NumCtrlID(1);
			Seeds2(i)= b0+incr*i;
		}
		else if(i< NumCtrlID(1)+NumCtrlID(3))
		{
			b0 = HTI_centerY-HTI_Hw/2;  //OriginCrdX from the web
		    incr = HTI_Hw/NumCtrlID(3);
			Seeds2(i)= b0+incr*(i-NumCtrlID(1));
		}
		else
		{
			b0 = HTI_centerY+HTI_Hw/2;  //OriginCrdY
		    incr = HTI_Tf/NumCtrlID(1);
			Seeds2(i)= b0+incr*(i-NumCtrlID(1)-NumCtrlID(3));
			}
	}
	
	for(int i =0; i<= NumCtrlID(4); i++){
        double c0,incr;
		Seeds3.resize(NumCtrlID(4)+1);
		 c0 = HTI_centerZ-HTI_Len/2;  //OriginCrdZ
		 incr = HTI_Len/NumCtrlID(4);
		 Seeds3(i)= c0+incr*i;
		
	}
	//Making seeds for slab along x
	for(int i =0; i<= NumCtrlID(5); i++){
		double a0,incr;
		Seeds4.resize(NumCtrlID(5)+1);
		if (i< (NumCtrlID(5)-NumCtrlID(0))/2)
		{
		    a0 = HTI_centerX-SlabW/2;  //OriginCrdX
		    incr = (SlabW-HTI_Bf)/(NumCtrlID(5)-NumCtrlID(0));
			Seeds4(i)= a0+incr*i;
		}
		
		else if (i< (NumCtrlID(5)-NumCtrlID(0))/2+(NumCtrlID(0)-NumCtrlID(2))/2)
		{
		    a0 = HTI_centerX-HTI_Bf/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds4(i)= a0+incr*(i-(NumCtrlID(5)-NumCtrlID(0))/2);
		}
		else if(i< (NumCtrlID(5)-NumCtrlID(0))/2+ (NumCtrlID(0)+NumCtrlID(2))/2)
		{
			a0 = HTI_centerX-HTI_Tw/2;  //OriginCrdX from the web
		    incr = HTI_Tw/NumCtrlID(2);
			Seeds4(i)= a0+incr*(i-(NumCtrlID(0)-NumCtrlID(2))/2-(NumCtrlID(5)-NumCtrlID(0))/2);
		}
		else if(i<(NumCtrlID(5)+NumCtrlID(0))/2)
		{
			a0 = HTI_centerX+HTI_Tw/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds4(i)= a0+incr*(i-(NumCtrlID(0)+NumCtrlID(2))/2-(NumCtrlID(5)-NumCtrlID(0))/2);
		}
		else 
		{
			a0 = HTI_centerX+HTI_Bf/2;  //OriginCrdX
		    incr = (SlabW-HTI_Bf)/(NumCtrlID(5)-NumCtrlID(0));
			Seeds4(i)= a0+incr*(i-(NumCtrlID(5)+NumCtrlID(0))/2);
			
		}
        
		
	}
	//making seeds for slab along y 
	for(int i =0; i<= NumCtrlID(6); i++){
        double b0,incr;
		 b0 = HTI_centerY+HTI_Hw/2+HTI_Tf;  //OriginCrdy
		 incr = SlabH/NumCtrlID(6);
		 Seeds5(i)= b0+incr*i;
		
	}
	
	
	return result;
}


const Vector& 
Simple_Composite3D::GetSeeds(int SeedTag)
{
	if(SeedTag==1){ 
		return Seeds1;
	}
	else if(SeedTag==2){
		return Seeds2;
	}
	else if(SeedTag==3){
		return Seeds3;
	}
	else if(SeedTag==4){
		return Seeds4;
	}
	else if(SeedTag==5){
		return Seeds5;
	}
	
		
	else 
		{
			opserr<<"Warning: invalid Seedtag For Isection"<<endln;
		}

}


const ID& 
Simple_Composite3D::GetNumCtrlID(void)
{
	return NumCtrlID;
}

int Simple_Composite3D::GetNumofNodes(void)
{
int NumNodes = ((NumCtrlID(0)+1)*(NumCtrlID(1)+1)+(NumCtrlID(0)+1)*NumCtrlID(1)+(NumCtrlID(2)+1)*(NumCtrlID(3)-1)+(NumCtrlID(6)+1)*(NumCtrlID(5)+1))*(NumCtrlID(4)+1);
	return NumNodes;
}

int Simple_Composite3D::GetNumofEles(void)
{
	int NumEles = (NumCtrlID(0)*NumCtrlID(1)*2+NumCtrlID(2)*NumCtrlID(3)+NumCtrlID(5)*NumCtrlID(6))*NumCtrlID(4);
	return NumEles;
}

int Simple_Composite3D::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs)
{
	double OriginLoc1 = 0;
	double OriginLoc2 = 0;
	if (OriginLocs.Size() == 1) {
		OriginLoc1 = OriginLocs(0);
	}
	else if (OriginLocs.Size() == 2) {
		OriginLoc1 = OriginLocs(0);
		OriginLoc2 = OriginLocs(1);
	}

	int OriginNodeTag = theHTDomain->getNumNodes() + 1;

	int NumCtrX = this->GetNumCtrlID()(0);
	int NumCtrY = this->GetNumCtrlID()(1);
	int NumCtrX_Web = this->GetNumCtrlID()(2);
	int NumCtrY_Web = this->GetNumCtrlID()(3);
	int NumCtrZ = this->GetNumCtrlID()(4);
	int NumCtrX_Slab = this->GetNumCtrlID()(5);
	int NumCtrY_Slab = this->GetNumCtrlID()(6);

	//Cordinates declared for nodes being generated here
	double NodeCrdX, NodeCrdY, NodeCrdZ;
	//for loop with k along the Beam length
	for (int k = 0; k <= NumCtrZ; k++) {
		// Nodal Cordinate y
		NodeCrdZ = (this->GetSeeds(3))(k);
		if (fabs(NodeCrdZ) < 1e-10)
			NodeCrdZ = 0;

		int LayerOriginNodeTag = OriginNodeTag + ((NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY + (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1)) * k;


		//Generating Nodes for LowerFlange
		for (int i = 0; i <= NumCtrY; i++) {

			NodeCrdY = (this->GetSeeds(2))(i);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			for (int j = 0; j <= NumCtrX; j++) {
				int NodeTag = LayerOriginNodeTag + (NumCtrX + 1) * i + j;

				NodeCrdX = (this->GetSeeds(1))(j);
				if (fabs(NodeCrdX) < 1e-10)
					NodeCrdX = 0;

				
				HeatTransferNode* TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, NodeCrdZ);
				if (theHTDomain->addNode(TempNode) < 0) {
					opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
					return -1;
				}
				

			}
		}

		//Generating Nodes for Web
		for (int i = 0; i < NumCtrY_Web - 1; i++) {
			// Nodal Cordinate y
			NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + 1);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			for (int j = 0; j <= NumCtrX_Web; j++) {

				int NodeTag = LayerOriginNodeTag + (NumCtrX_Web + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1);
				// Nodal Cordinate x
				NodeCrdX = (this->GetSeeds(1))(j + (NumCtrX - NumCtrX_Web) / 2);
				if (fabs(NodeCrdX) < 1e-10)
					NodeCrdX = 0;

				HeatTransferNode* TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, NodeCrdZ);
				if (theHTDomain->addNode(TempNode) < 0) {
					opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
					return -1;
				}
					
			}
		}

		//Generating Nodes for UpperFlange(Top surface nodes not generated)
		for (int i = 0; i < NumCtrY; i++) {
			NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + NumCtrY_Web);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			for (int j = 0; j <= NumCtrX; j++) {
				int NodeTag = LayerOriginNodeTag + (NumCtrX + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1);
				NodeCrdX = (this->GetSeeds(1))(j);
				if (fabs(NodeCrdX) < 1e-10)
					NodeCrdX = 0;

				HeatTransferNode* TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, NodeCrdZ);
				if (theHTDomain->addNode(TempNode) < 0) {
					opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
					return -1;
				}
					
			}
		}

		//Generating Nodes for Slab section
		for (int i = 0; i <= NumCtrY_Slab; i++) {

			NodeCrdY = (this->GetSeeds(5))(i);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			for (int j = 0; j <= NumCtrX_Slab; j++) {
				int NodeTag = LayerOriginNodeTag + (NumCtrX_Slab + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY;
				NodeCrdX = (this->GetSeeds(4))(j);
				if (fabs(NodeCrdX) < 1e-10)
					NodeCrdX = 0;
				
				HeatTransferNode* TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, NodeCrdZ);
				if (theHTDomain->addNode(TempNode) < 0) {
					opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
					return -1;
				}
		
			}
		}
		//End of adding nodes for Composite_Slab


	}
	//end of for(k=0;k<NumCtrZ;k++)

}

int Simple_Composite3D::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
{
	bool PhaseTransformation = false;
	bool PhaseTransformation1 = false;

	if (EleParameters != 0) {
		if (EleParameters(0) == 1) {
			PhaseTransformation = true;
		}
		else {
			if (theHTMaterial1 != 0 && EleParameters.Size() > 1) {
				if (EleParameters(1) == 1)
					PhaseTransformation1 = true;
			}
		}
	}

	int OriginNodeTag = theHTDomain->getNumNodes() - (this->GetNumofNodes()) + 1;
	int OriginEleTag = theHTDomain->getNumElements() + 1;
	HeatTransferElement* TempEle = 0;

	int NumCtrX = this->GetNumCtrlID()(0);
	int NumCtrY = this->GetNumCtrlID()(1);
	int NumCtrX_Web = this->GetNumCtrlID()(2);
	int NumCtrY_Web = this->GetNumCtrlID()(3);
	int NumCtrZ = this->GetNumCtrlID()(4);
	int NumCtrX_Slab = this->GetNumCtrlID()(5);
	int NumCtrY_Slab = this->GetNumCtrlID()(6);


	int NumNodesPerLayer = (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY + (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1);
	int NumElesPerLayer = NumCtrX * NumCtrY * 2 + NumCtrX_Web * NumCtrY_Web + NumCtrX_Slab * NumCtrY_Slab;

	int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8;
	for (int k = 0; k < NumCtrZ; k++) {
		//along the length
	//Generating Elements for LowerFlange
		for (int i = 0; i < NumCtrY; i++) {
			for (int j = 0; j < NumCtrX; j++) {
				EleTag = OriginEleTag + NumElesPerLayer * k + NumCtrX * i + j;
				NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * i + j;
				NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * i + j + 1;
				NodeTag3 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (i + 1) + j + 1;
				NodeTag4 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (i + 1) + j;
				NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * i + j;
				NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * i + j + 1;
				NodeTag7 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (i + 1) + j + 1;
				NodeTag8 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (i + 1) + j;


				TempEle = new BrickEight(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8, *theHTMaterial, PhaseTransformation);

				if (theHTDomain->addElement(TempEle) < 0) {
					opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
					return -1;
				}
			}
		}

		//Generating Elements for Web
		int JunctionLTag = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) - (NumCtrX + 1) + (NumCtrX - NumCtrX_Web) / 2;
		int JunctionUTag = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX - NumCtrX_Web) / 2;

		for (int i = 0; i < NumCtrY_Web; i++) {
			for (int j = 0; j < NumCtrX_Web; j++) {
				EleTag = OriginEleTag + NumElesPerLayer * k + NumCtrX_Web * i + j + NumCtrX * NumCtrY;

				if (i == 0)
				{
					NodeTag1 = JunctionLTag + NumNodesPerLayer * k + j;
					NodeTag2 = JunctionLTag + NumNodesPerLayer * k + j + 1;
					NodeTag3 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
					NodeTag4 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;

					NodeTag5 = JunctionLTag + NumNodesPerLayer * (k + 1) + j;
					NodeTag6 = JunctionLTag + NumNodesPerLayer * (k + 1) + j + 1;
					NodeTag7 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
					NodeTag8 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;

				}

				else if (i == NumCtrY_Web - 1)
				{
					NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
					NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
					NodeTag3 = JunctionUTag + NumNodesPerLayer * k + j + 1;
					NodeTag4 = JunctionUTag + NumNodesPerLayer * k + j;
					NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
					NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
					NodeTag7 = JunctionUTag + NumNodesPerLayer * (k + 1) + j + 1;
					NodeTag8 = JunctionUTag + NumNodesPerLayer * (k + 1) + j;

				}
				else
				{
					NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
					NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
					NodeTag3 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
					NodeTag4 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;
					NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
					NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
					NodeTag7 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
					NodeTag8 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;


				}

				TempEle = new BrickEight(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8, *theHTMaterial, PhaseTransformation);

				if (theHTDomain->addElement(TempEle) < 0) {
					opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
					return -1;
				}

			}
		}

		//Generating Elements for UpperFlange
		for (int i = 0; i < NumCtrY; i++) {
			for (int j = 0; j < NumCtrX; j++) {
				EleTag = OriginEleTag + NumElesPerLayer * k + NumCtrX * i + j + NumCtrX * NumCtrY + NumCtrX_Web * NumCtrY_Web;
				int JunctionUTagSlab = OriginNodeTag + NumNodesPerLayer - (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1) + (NumCtrX_Slab - NumCtrX) / 2;

				if (i == NumCtrY - 1)
				{
					NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
					NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
					NodeTag3 = JunctionUTagSlab + NumNodesPerLayer * k + j + 1;
					NodeTag4 = JunctionUTagSlab + NumNodesPerLayer * k + j;

					NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
					NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
					NodeTag7 = JunctionUTagSlab + NumNodesPerLayer * (k + 1) + j + 1;
					NodeTag8 = JunctionUTagSlab + NumNodesPerLayer * (k + 1) + j;


				}
				else
				{
					NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
					NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
					NodeTag3 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j + 1;
					NodeTag4 = OriginNodeTag + NumNodesPerLayer * k + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j;

					NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
					NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
					NodeTag7 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j + 1;
					NodeTag8 = OriginNodeTag + NumNodesPerLayer * (k + 1) + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j;

				}

				TempEle = new BrickEight(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8, *theHTMaterial, PhaseTransformation);

				if (theHTDomain->addElement(TempEle) < 0) {
					opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
					return -1;
				}
			}
		}
		//end of for(int i=0;i<NumCtrY; i++);


		//Generating Elements for CompositeSlab
		for (int i = 0; i < NumCtrY_Slab; i++) {
			for (int j = 0; j < NumCtrX_Slab; j++) {
				EleTag = OriginEleTag + NumElesPerLayer * (k + 1) - NumCtrX_Slab * NumCtrY_Slab + NumCtrX_Slab * i + j;
				int OriginNodeTagforSlab = NumNodesPerLayer - (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1);

				NodeTag1 = OriginNodeTag + NumNodesPerLayer * k + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j;
				NodeTag2 = OriginNodeTag + NumNodesPerLayer * k + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j + 1;
				NodeTag3 = OriginNodeTag + NumNodesPerLayer * k + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j + 1;
				NodeTag4 = OriginNodeTag + NumNodesPerLayer * k + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j;

				NodeTag5 = OriginNodeTag + NumNodesPerLayer * (k + 1) + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j;
				NodeTag6 = OriginNodeTag + NumNodesPerLayer * (k + 1) + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j + 1;
				NodeTag7 = OriginNodeTag + NumNodesPerLayer * (k + 1) + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j + 1;
				NodeTag8 = OriginNodeTag + NumNodesPerLayer * (k + 1) + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j;


				TempEle = new BrickEight(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8, *theHTMaterial, PhaseTransformation);

				if (theHTDomain->addElement(TempEle) < 0) {
					opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
					return -1;
				}
			}
		}
		//end of adding elements for composite slab

	}
	//end of for (int k = 0; k<NumCtrZ;k++);


	return 0;
}