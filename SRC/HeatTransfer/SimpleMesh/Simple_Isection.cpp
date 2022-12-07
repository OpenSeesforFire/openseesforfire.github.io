//This code is based on OpenSees Framework 
//Written by Liming Jiang(University of Edinburgh)


#include <Simple_Isection.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <QuadFour.h>


Simple_Isection::Simple_Isection(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw, double HTI_UBf,double HTI_UTf)
:Simple_Entity(tag,1),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf),
HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw), HTI_UBf(HTI_UBf), HTI_UTf(HTI_UTf),NumCtrlID(4),Seeds1(200),Seeds2(200)
{
	
}

Simple_Isection::Simple_Isection(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw)
:Simple_Entity(tag,1),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf), 
HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw),HTI_UBf(HTI_Bf), HTI_UTf(HTI_Tf),NumCtrlID(4),Seeds1(200),Seeds2(200)
{
	
}


Simple_Isection::~Simple_Isection()
{
	
	
}

int 
Simple_Isection::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
	double EleX= MeshCtrls(0);
	double EleY= MeshCtrls(1);
	double EleX_Web= MeshCtrls(2);
	double EleY_Web= MeshCtrls(3);
	NumCtrlID(0)= (HTI_Bf-HTI_Tw)/EleX+HTI_Tw/EleX_Web+0.5;  //Num of elements generated in flange along x-axis
	NumCtrlID(1)= HTI_Tf/EleY+0.5;                           //Num of elements generated in flange along Y-axis
	NumCtrlID(2) = HTI_Tw/EleX_Web+0.5;                      //Num of elements generated in Web along x-axis
	NumCtrlID(3) = HTI_Hw/EleY_Web+0.5;                       //Num of elements generated in Web along y-axis
	if ((NumCtrlID(0)-NumCtrlID(2))%2!=0)
		{
			opserr<< "WARNING: Invalid dimension input for Isection " <<this->getTag()<<endln;
			return -1;
		}
	
	if(NumCtrlID(0)==0||NumCtrlID(1)==0||NumCtrlID(2)==0||NumCtrlID(3)==0)
		{
			opserr<<"Error: NumCtrlID in Isection: "<< this->getTag()<<" has 0 vlaue"<<endln;
			return -2;
		}
	return 0;

}


bool 
Simple_Isection::InitialSeeds(void)
{
	bool result = true;

	if (NumCtrlID(0)>100||NumCtrlID(1)>100||NumCtrlID(2)>100||NumCtrlID(3)>100)
		{
		opserr<< "One of The mesh-seed vetor is overflow";
		result= false;
		}

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
	
	return result;
}


const Vector& 
Simple_Isection::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ 
			return Seeds1;
		}
	else if(SeedTag==2)
		{
			return Seeds2;
		}
	else 
		{
			opserr<<"Warning: invalid facetag For Isection"<<endln;
		}

}


const ID& 
Simple_Isection::GetNumCtrlID(void)
{
	return NumCtrlID;
}

int Simple_Isection::GetNumofNodes(void)
{
	int NumNodes = (NumCtrlID(0)+1)*(NumCtrlID(1)+1)*2+(NumCtrlID(2)+1)*(NumCtrlID(3)-1);
	return NumNodes;
}

int Simple_Isection::GetNumofEles(void)
{
	int NumEles = NumCtrlID(0)*NumCtrlID(1)*2+NumCtrlID(2)*NumCtrlID(3);
	return NumEles;
}

int Simple_Isection::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs) 
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
	double NodeCrdX, NodeCrdY;
	//Generating Nodes for LowerFlange
	for (int i = 0; i <= NumCtrY; i++) {
		for (int j = 0; j <= NumCtrX; j++) {
			int NodeTag = (NumCtrX + 1) * i + j;
			NodeCrdX = (this->GetSeeds(1))(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = (this->GetSeeds(2))(i);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			
				HeatTransferNode* TempNode = 0;
				if (OriginLocs.Size() == 1) {
					TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
				}
				else {
					TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY);
				}

				if (theHTDomain->addNode(TempNode) < 0) {
					opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
					return -1;
				}
			

		}
	}

	//Generating Nodes for Web
	for (int i = 0; i < NumCtrY_Web - 1; i++) {
		for (int j = 0; j <= NumCtrX_Web; j++) {

			int NodeTag = (NumCtrX_Web + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1);
			NodeCrdX = (this->GetSeeds(1))(j + (NumCtrX - NumCtrX_Web) / 2);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + 1);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			HeatTransferNode* TempNode = 0;
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}

			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
			
		}
	}

	//Generating Nodes for UpperFlange
	for (int i = 0; i <= NumCtrY; i++) {
		for (int j = 0; j <= NumCtrX; j++) {
			int NodeTag = (NumCtrX + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1);
			NodeCrdX = (this->GetSeeds(1))(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + NumCtrY_Web);
			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			
			HeatTransferNode* TempNode = 0;
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(OriginNodeTag + NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}
			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
		
		}
	}

}

int Simple_Isection::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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

	//Generating Elements for LowerFlange
	int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;
	for (int i = 0; i < NumCtrY; i++) {
		for (int j = 0; j < NumCtrX; j++) {
			EleTag = OriginEleTag + NumCtrX * i + j;
			NodeTag1 = OriginNodeTag + (NumCtrX + 1) * i + j;
			NodeTag2 = OriginNodeTag + (NumCtrX + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (i + 1) + j;
			

			TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);

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
			EleTag = OriginEleTag + NumCtrX_Web * i + j + NumCtrX * NumCtrY;

			if (i == 0)
			{
				NodeTag1 = JunctionLTag + j;
				NodeTag2 = JunctionLTag + j + 1;
				NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
				NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;
			}
			else if (i == NumCtrY_Web - 1)
			{
				NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
				NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
				NodeTag3 = JunctionUTag + j + 1;
				NodeTag4 = JunctionUTag + j;

			}
			else
			{
				NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j;
				NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (i - 1) + j + 1;
				NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j + 1;
				NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * i + j;

			}

			TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}

		}
	}

	//Generating Elements for UpperFlange

	for (int i = 0; i < NumCtrY; i++) {
		for (int j = 0; j < NumCtrX; j++) {
			EleTag = OriginEleTag + NumCtrX * i + j + NumCtrX * NumCtrY + NumCtrX_Web * NumCtrY_Web;
			NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
			NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j;


			TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}
		}
	}
	//end of for(int i=0;i<NumCtrY; i++);


	return 0;
}