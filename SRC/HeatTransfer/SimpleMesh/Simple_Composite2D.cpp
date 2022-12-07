#include <Simple_Composite2D.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <QuadFour.h>

/*
Simple_Composite2D::Simple_Composite2D(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw, double HTI_UBf,double HTI_UTf)
:Simple_Entity(tag,1),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf),
HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw), HTI_UBf(HTI_UBf), HTI_UTf(HTI_UTf),NumCtrlID(4),Seeds1(200),Seeds2(200)
{
	
}
*/

Simple_Composite2D::Simple_Composite2D(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf,
									   double HTI_Tf, double HTI_Tw, double HTI_Hw,double slabW, double slabH)
:Simple_Entity(tag,6),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf), 
HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw),HTI_UBf(HTI_Bf), HTI_UTf(HTI_Tf),NumCtrlID(6),SlabH(slabH), SlabW(slabW),
Seeds1(0),Seeds2(0), Seeds3(0), Seeds4(0)
{
	
}


Simple_Composite2D::~Simple_Composite2D()
{
	

}

int 
Simple_Composite2D::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
	double EleX= MeshCtrls(0);
	double EleY= MeshCtrls(1);
	double EleX_Web= MeshCtrls(2);
	double EleY_Web= MeshCtrls(3);
	double EleX_Slab = MeshCtrls(4);
	double EleY_Slab = MeshCtrls(5);
	
	NumCtrlID(0)= (HTI_Bf-HTI_Tw)/EleX+HTI_Tw/EleX_Web+0.5;  //Num of elements generated in flange along x-axis
	NumCtrlID(1)= HTI_Tf/EleY+0.5;                           //Num of elements generated in flange along Y-axis
	NumCtrlID(2) = HTI_Tw/EleX_Web+0.5;                      //Num of elements generated in Web along x-axis
	NumCtrlID(3) = HTI_Hw/EleY_Web+0.5;                       //Num of elements generated in Web along y-axis
    NumCtrlID(4) = (SlabW-HTI_Bf)/EleX_Slab+(HTI_Bf-HTI_Tw)/EleX+HTI_Tw/EleX_Web+0.5;    //Number of elements generated in slab along local x axis;
	NumCtrlID(5) = SlabH/EleY_Slab;		                     //Number of elements along y-axis

	if ((NumCtrlID(0)-NumCtrlID(2))%2!=0)
		{
			opserr<< "Error: Invalid dimension input for Isection " <<this->getTag()<<endln;
			return -1;
		}
	
	if(NumCtrlID(0)==0||NumCtrlID(1)==0||NumCtrlID(2)==0||NumCtrlID(3)==0||NumCtrlID(4)==0||NumCtrlID(5)==0)
		{
			opserr<<"Error: NumCtrlID in Isection: "<< this->getTag()<<" has 0 vlaue"<<endln;
			return -2;
		}
	return 0;

}


bool 
Simple_Composite2D::InitialSeeds(void)
{
	bool result = true;

	Seeds1.resize(NumCtrlID(0)+1);
	Seeds2.resize(NumCtrlID(1)*2+NumCtrlID(3)+1);
	Seeds3.resize(NumCtrlID(3)+1);
	Seeds4.resize(NumCtrlID(4)+1);

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

	//Making seeds for slab along x
	for(int i =0; i<= NumCtrlID(4); i++){
		double a0,incr;
		Seeds3.resize(NumCtrlID(4)+1);
		if (i< (NumCtrlID(4)-NumCtrlID(0))/2)
		{
		    a0 = HTI_centerX-SlabW/2;  //OriginCrdX
		    incr = (SlabW-HTI_Bf)/(NumCtrlID(4)-NumCtrlID(0));
			Seeds3(i)= a0+incr*i;
		}
		
		else if (i< (NumCtrlID(4)-NumCtrlID(0))/2+(NumCtrlID(0)-NumCtrlID(2))/2)
		{
		    a0 = HTI_centerX-HTI_Bf/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds3(i)= a0+incr*(i-(NumCtrlID(4)-NumCtrlID(0))/2);
		}
		else if(i< (NumCtrlID(4)-NumCtrlID(0))/2+ (NumCtrlID(0)+NumCtrlID(2))/2)
		{
			a0 = HTI_centerX-HTI_Tw/2;  //OriginCrdX from the web
		    incr = HTI_Tw/NumCtrlID(2);
			Seeds3(i)= a0+incr*(i-(NumCtrlID(0)-NumCtrlID(2))/2-(NumCtrlID(4)-NumCtrlID(0))/2);
		}
		else if(i<(NumCtrlID(4)+NumCtrlID(0))/2)
		{
			a0 = HTI_centerX+HTI_Tw/2;  //OriginCrdX
		    incr = (HTI_Bf-HTI_Tw)/(NumCtrlID(0)-NumCtrlID(2));
			Seeds3(i)= a0+incr*(i-(NumCtrlID(0)+NumCtrlID(2))/2-(NumCtrlID(4)-NumCtrlID(0))/2);
		}
		else 
		{
			a0 = HTI_centerX+HTI_Bf/2;  //OriginCrdX
		    incr = (SlabW-HTI_Bf)/(NumCtrlID(4)-NumCtrlID(0));
			Seeds3(i)= a0+incr*(i-(NumCtrlID(4)+NumCtrlID(0))/2);
			
		}
        
		
	}
	//making seeds for slab along y 
	for(int i =0; i<= NumCtrlID(5); i++){
        double b0,incr;
		Seeds4.resize(NumCtrlID(5)+1);
		 b0 = HTI_centerY+HTI_Hw/2+HTI_Tf;  //OriginCrdy
		 incr = SlabH/NumCtrlID(5);
		 Seeds4(i)= b0+incr*i;
		
	}
	
	return result;
}


const Vector& 
Simple_Composite2D::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ 
			return Seeds1;
		}
	else if(SeedTag==2)
		{
			return Seeds2;
		}
	else if(SeedTag==3)
	{
		return Seeds3;
	}
	else if(SeedTag==4)
	{
		return Seeds4;
	}
			
		
	else 
		{
			opserr<<"Warning: invalid facetag For Isection"<<endln;
		}

}


const ID& 
Simple_Composite2D::GetNumCtrlID(void)
{
	return NumCtrlID;
}

int Simple_Composite2D::GetNumofNodes(void)
{
	int NumNodes = ((NumCtrlID(0)+1)*(NumCtrlID(1)+1)+(NumCtrlID(0)+1)*NumCtrlID(1)+
		(NumCtrlID(2)+1)*(NumCtrlID(3)-1)+(NumCtrlID(4)+1)*(NumCtrlID(5)+1));
	return NumNodes;
}

int Simple_Composite2D::GetNumofEles(void)
{
	int NumEles = (NumCtrlID(0)*NumCtrlID(1)*2+NumCtrlID(2)*NumCtrlID(3)+NumCtrlID(5)*NumCtrlID(4));
	return NumEles;
}


int Simple_Composite2D::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs) 
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
	int NumCtrX_Slab = this->GetNumCtrlID()(4);
	int NumCtrY_Slab = this->GetNumCtrlID()(5);

	//Cordinates declared for nodes being generated here
	double NodeCrdX, NodeCrdY;
	//for loop with k along the Beam length
	HeatTransferNode* TempNode = 0;

	//Generating Nodes for LowerFlange
	for (int i = 0; i <= NumCtrY; i++) {

		NodeCrdY = (this->GetSeeds(2))(i);
		if (fabs(NodeCrdY) < 1e-10)
			NodeCrdY = 0;

		for (int j = 0; j <= NumCtrX; j++) {
			int NodeTag = OriginNodeTag + (NumCtrX + 1) * i + j;

			NodeCrdX = (this->GetSeeds(1))(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

		
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}

			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
				//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;

		}
	}

	//Generating Nodes for Web
	for (int i = 0; i < NumCtrY_Web - 1; i++) {
		// Nodal Cordinate y
		NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + 1);
		if (fabs(NodeCrdY) < 1e-10)
			NodeCrdY = 0;

		for (int j = 0; j <= NumCtrX_Web; j++) {

			int NodeTag = OriginNodeTag + (NumCtrX_Web + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1);
			// Nodal Cordinate x
			NodeCrdX = (this->GetSeeds(1))(j + (NumCtrX - NumCtrX_Web) / 2);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}

			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
				//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
			
		}
	}

	//Generating Nodes for UpperFlange(Top surface nodes not generated)
	for (int i = 0; i < NumCtrY; i++) {
		NodeCrdY = (this->GetSeeds(2))(i + NumCtrY + NumCtrY_Web);
		if (fabs(NodeCrdY) < 1e-10)
			NodeCrdY = 0;

		for (int j = 0; j <= NumCtrX; j++) {
			int NodeTag = OriginNodeTag + (NumCtrX + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1);
			NodeCrdX = (this->GetSeeds(1))(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}

			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
				//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
			
		}
	}

	//Generating Nodes for Slab section
	for (int i = 0; i <= NumCtrY_Slab; i++) {

		NodeCrdY = (this->GetSeeds(4))(i);
		if (fabs(NodeCrdY) < 1e-10)
			NodeCrdY = 0;

		for (int j = 0; j <= NumCtrX_Slab; j++) {
			int NodeTag = OriginNodeTag + (NumCtrX_Slab + 1) * i + j + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY;
			NodeCrdX = (this->GetSeeds(3))(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY);
			}

			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
				//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
			
		}
	}
	//End of adding nodes for Composite_Slab


}

int Simple_Composite2D::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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
	int NumCtrX_Slab = this->GetNumCtrlID()(4);
	int NumCtrY_Slab = this->GetNumCtrlID()(5);


	int NumNodesPerLayer = (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY + (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1);
	int NumElesPerLayer = NumCtrX * NumCtrY * 2 + NumCtrX_Web * NumCtrY_Web + NumCtrX_Slab * NumCtrY_Slab;

	//Generating elements for lower flange
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
			int JunctionUTagSlab = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * NumCtrY + (NumCtrX_Slab - NumCtrX) / 2;

			if (i == NumCtrY - 1)
			{
				NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
				NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
				NodeTag3 = JunctionUTagSlab + j + 1;
				NodeTag4 = JunctionUTagSlab + j;


			}
			else
			{
				NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j;
				NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * i + j + 1;
				NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j + 1;
				NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) + (NumCtrX_Web + 1) * (NumCtrY_Web - 1) + (NumCtrX + 1) * (i + 1) + j;

			}



			TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);

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
			EleTag = OriginEleTag + NumCtrX * NumCtrY * 2 + NumCtrX_Web * NumCtrY_Web + NumCtrX_Slab * i + j;
			int OriginNodeTagforSlab = NumNodesPerLayer - (NumCtrX_Slab + 1) * (NumCtrY_Slab + 1);

			NodeTag1 = OriginNodeTag + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j;
			NodeTag2 = OriginNodeTag + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + OriginNodeTagforSlab + (NumCtrX_Slab + 1) * (i + 1) + j;



			TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial1, PhaseTransformation1);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}

		}
	}
	//end of adding elements for composite slab



	return 0;
}