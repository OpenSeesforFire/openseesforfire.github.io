//This code is based on OpenSees Framework 
//Written by Liming Jiang(University of Edinburgh)

#include <Simple_IsecProtected.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>'
#include <QuadFour.h>


Simple_IsecProtected::Simple_IsecProtected(int tag, double HTI_centerX, double HTI_centerY, double HTI_Bf, double HTI_Tf, double HTI_Tw, double HTI_Hw, double HTI_coat)
:Simple_Entity(tag,11),HTI_centerX(HTI_centerX), HTI_centerY(HTI_centerY), HTI_Bf(HTI_Bf), HTI_Tf(HTI_Tf), HTI_Coat (HTI_coat),
HTI_Tw(HTI_Tw), HTI_Hw(HTI_Hw),HTI_UBf(HTI_Bf), HTI_UTf(HTI_Tf),NumCtrlID(5),Seeds1(200),Seeds2(200),
EleX(0),EleY(0),EleX_Web(0), EleY_Web(0), Ele_coat(0)
{
	
}


Simple_IsecProtected::~Simple_IsecProtected()
{
	
	
}

int 
Simple_IsecProtected::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
	 EleX= MeshCtrls(0);
	 EleY= MeshCtrls(1);
	 EleX_Web= MeshCtrls(2);
	 EleY_Web= MeshCtrls(3);
	Ele_coat = MeshCtrls(4);
  
	NumCtrlID(0)= (HTI_Bf-HTI_Tw - HTI_Coat*2)/EleX+HTI_Tw/EleX_Web+HTI_Coat/Ele_coat*2+ 0.2;  //Num of elements generated in flange along x-axis
	NumCtrlID(1)= HTI_Tf/EleY+0.2;                           //Num of elements generated in flange along Y-axis
	NumCtrlID(2) = HTI_Tw/EleX_Web+0.2;                      //Num of elements generated in Web along x-axis
	NumCtrlID(3) = (HTI_Hw-HTI_Coat*2)/EleY_Web+ HTI_Coat/Ele_coat*2 +0.2;    //Num of elements generated in Web along y-axis
	NumCtrlID(4) = HTI_Coat/Ele_coat+0.2;    
	
	EleX= (HTI_Bf-HTI_Tw - HTI_Coat*2)/(NumCtrlID(0)-NumCtrlID(2)-NumCtrlID(4)*2);
	 EleY= HTI_Tf/NumCtrlID(1);
	 EleX_Web= HTI_Tw/NumCtrlID(2);
	 EleY_Web= (HTI_Hw-HTI_Coat*2)/(NumCtrlID(3)-NumCtrlID(4)*2);
	Ele_coat = HTI_Coat/NumCtrlID(4);
	
	//Num of elements generated in outercoating along y-axis
#ifdef _DEBUG
	//opserr<<"IsecPro : NumberControlID" <<NumCtrlID<<endln;
#endif
  
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
Simple_IsecProtected::InitialSeeds(void)
{
	bool result = true;

	if (NumCtrlID(0)>200||NumCtrlID(1)>100||NumCtrlID(2)>100||NumCtrlID(3)>100)
		{
		opserr<< "One of The mesh-seed vetor is overflow";
		result= false;
		}
  
  
// along x
	for(int i =0; i<=NumCtrlID(0); i++){
    double a0,incr;
  
		if (i< (NumCtrlID(0)-NumCtrlID(2))/2-NumCtrlID(4))
		{
		    a0 = HTI_centerX-HTI_Bf/2;  //OriginCrdX for uncoated
		    incr = EleX;
        Seeds1(i)= a0+incr*i;
      
		}
    else if (i< (NumCtrlID(0)-NumCtrlID(2))/2)
    {
		    a0 = HTI_centerX-HTI_Tw/2- HTI_Coat;  //OriginCrdX for coating
		    incr =  Ele_coat;
      Seeds1(i)= a0+incr* (i-((NumCtrlID(0)-NumCtrlID(2))/2-NumCtrlID(4)));
    }
    
		else if(i< (NumCtrlID(0) + NumCtrlID(2))/2)
		{
        a0 = HTI_centerX - HTI_Tw/2;  //OriginCrdX from web
		    incr = EleX_Web;
			  Seeds1(i)= a0+incr*(i-(NumCtrlID(0)-NumCtrlID(2))/2);
		}
    
    else if (i< (NumCtrlID(0)+NumCtrlID(2))/2+ NumCtrlID(4))
    {
		    a0 = HTI_centerX+HTI_Tw/2;  //OriginCrdX for coating
		    incr = Ele_coat;
        Seeds1(i)= a0+incr* (i-(NumCtrlID(0)+NumCtrlID(2))/2 );
    }
    
		else
		{
        a0 = HTI_centerX+HTI_Tw/2 + HTI_Coat;  //OriginCrdX
		    incr = EleX;
        Seeds1(i)= a0+incr*(i-(NumCtrlID(0)+NumCtrlID(2))/2 - NumCtrlID(4));
		}
	
	}
//along y
	for(int i =0; i<= NumCtrlID(1)*2+NumCtrlID(3)+2* NumCtrlID(4); i++){
        double b0,incr;
		
		if (i< NumCtrlID(4))
		{
		    b0 = HTI_centerY-HTI_Hw/2-HTI_Tf - HTI_Coat;  //OriginCrdY for bottom coat
		    incr =  Ele_coat;
			Seeds2(i)= b0+incr*i;
		}
		else if(i< NumCtrlID(1)+NumCtrlID(4))
		{
			b0 = HTI_centerY-HTI_Hw/2 -HTI_Tf ;  //OriginCrdX from the bottom flange
			incr = EleY;
			Seeds2(i)= b0+incr*(i-NumCtrlID(4));
		}
    
		else if(i< NumCtrlID(1) + NumCtrlID(4)*2 )
		{
			b0 = HTI_centerY-HTI_Hw/2;  //OriginCrdX from the second bottom coat
		    incr =  Ele_coat;
			Seeds2(i)= b0+incr*(i-(NumCtrlID(1)+NumCtrlID(4)));
		}

		else if(i< NumCtrlID(1)+NumCtrlID(3))
		{
			b0 = HTI_centerY-HTI_Hw/2 +HTI_Coat;  //OriginCrdX from the web
		    incr = EleY_Web;
			Seeds2(i)= b0+incr*(i-(NumCtrlID(1)+NumCtrlID(4)*2));
		}

		else if(i< NumCtrlID(1) + NumCtrlID(3) + NumCtrlID(4) )
		{
			b0 = HTI_centerY+ HTI_Hw/2- HTI_Coat;  //OriginCrdX from the second bottom coat
		    incr =  Ele_coat;
			Seeds2(i)= b0+incr*(i-(NumCtrlID(1)+NumCtrlID(3)));
		}

		else if(i< NumCtrlID(1)*2 + NumCtrlID(3) + NumCtrlID(4))
		{
			b0 = HTI_centerY+ HTI_Hw/2 ;  //OriginCrdX from the bottom flange
			incr = EleY;
			Seeds2(i)= b0+incr*(i-(NumCtrlID(1) + NumCtrlID(3) + NumCtrlID(4)));
		}
		
		else if(i<= NumCtrlID(1)*2 + NumCtrlID(3) + NumCtrlID(4)*2 )
		{
			b0 = HTI_centerY+ HTI_Hw/2+ HTI_Tf ;  //OriginCrdX from the second bottom coat
		    incr =  Ele_coat;
			Seeds2(i)= b0+incr*(i-(NumCtrlID(1)*2 + NumCtrlID(3) + NumCtrlID(4)));
		}
		else
		{
			opserr<<"out of bounds"<<endln;
		}
	}
  
	//#ifdef _DEBUG
	 //opserr<<Seeds1<<endln;	
	  //opserr<<Seeds2<<endln;
//#endif
  /*we need seeds for coating along y axis
  for(int i =0; i<= NumCtrlID(4)*2; i++){
    double b0,incr;
    
    if (i< NumCtrlID(4))
    {
		    b0 = HTI_centerY-HTI_Hw/2-HTI_Tf-HTI_Coat;  //OriginCrdY
		    incr = HTI_Coat/NumCtrlID(4);
        Seeds3(i)= b0+incr*i;
    }
    else if(i< NumCtrlID(4)*2)
    {
        b0 = HTI_centerY + HTI_Hw/2 + HTI_Tf;  //OriginCrdX from the web
		    incr = HTI_Coat/NumCtrlID(4);
        Seeds3(i)= b0+incr*(i-NumCtrlID(4));
    }
    else
    {
      opserr<<"out of range"<<endln;
    }

  }
   */
  
	
	return result;
}


const Vector& 
Simple_IsecProtected::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ 
#ifdef _DEBUG
		//opserr<<Seeds1<<endln;	
#endif
		return Seeds1;
		}
	else if(SeedTag==2)
		{
#ifdef _DEBUG
		//opserr<<Seeds2<<endln;	
#endif
			return Seeds2;
		}
	else 
		{
			opserr<<"Warning: invalid facetag For Isection"<<endln;
		}

}


const ID& 
Simple_IsecProtected::GetNumCtrlID(void)
{
	return NumCtrlID;
}


int Simple_IsecProtected::GetNumofNodes(void)
{
  opserr<<NumCtrlID<<endln;
  int NumNodes = (NumCtrlID(0)+1)*((NumCtrlID(4))*2 + NumCtrlID(1)+1)*2;
	  
  NumNodes += (NumCtrlID(3)- NumCtrlID(4)*2 -1)*(NumCtrlID(2)+NumCtrlID(4)*2+1);
                                                
  return NumNodes;
}


int Simple_IsecProtected::GetNumofEles(void)
{
	int NumEles = NumCtrlID(0)*(NumCtrlID(4)*2 + NumCtrlID(3) + NumCtrlID(1)*2) - ((NumCtrlID(0)-NumCtrlID(2)-2*(NumCtrlID(4))*(NumCtrlID(3) - NumCtrlID(4)*2)));
	return NumEles;
  
  
}


int Simple_IsecProtected::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs)
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
	int NumCtr_Coat = this->GetNumCtrlID()(4);

	int NumX = NumCtrX;
	int NumY = NumCtrY + NumCtr_Coat * 2;
	int NumXweb = NumCtrX_Web + NumCtr_Coat * 2;
	int NumYweb = NumCtrY_Web - NumCtr_Coat * 2;

	double NodeCrdX, NodeCrdY;
	Vector Seeds1 = this->GetSeeds(1);
	Vector Seeds2 = this->GetSeeds(2);
	//Generating Nodes for Lower great Flange
	for (int i = 0; i <= NumY; i++) {
		for (int j = 0; j <= NumX; j++) {
			int NodeTag = (NumX + 1) * i + j;
			NodeCrdX = Seeds1(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = Seeds2(i);
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
	for (int i = 0; i < NumYweb - 1; i++) {
		for (int j = 0; j <= NumXweb; j++) {

			int NodeTag = (NumXweb + 1) * i + j + (NumX + 1) * (NumY + 1);
			NodeCrdX = Seeds1(j + (NumX - NumXweb) / 2);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = Seeds2(i + NumY + 1);
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
	for (int i = 0; i <= NumY; i++) {
		for (int j = 0; j <= NumX; j++) {
			int NodeTag = (NumX + 1) * i + j + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1);
			NodeCrdX = Seeds1(j);
			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			NodeCrdY = Seeds2(i + NumY + NumYweb);
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

	return 0;

}



int Simple_IsecProtected::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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
	int NumCtr_Coat = this->GetNumCtrlID()(4);

	int NumX = NumCtrX;
	int NumY = NumCtrY + NumCtr_Coat * 2;
	int NumXweb = NumCtrX_Web + NumCtr_Coat * 2;
	int NumYweb = NumCtrY_Web - NumCtr_Coat * 2;

	//Generating Elements for LowerFlange
	int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;

	for (int i = 0; i < NumY; i++) {
		for (int j = 0; j < NumX; j++) {
			EleTag = OriginEleTag + NumX * i + j;
			NodeTag1 = OriginNodeTag + (NumX + 1) * i + j;
			NodeTag2 = OriginNodeTag + (NumX + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + (NumX + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + (NumX + 1) * (i + 1) + j;

			

			if (i >= NumCtr_Coat && i < NumY - NumCtr_Coat)
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);
			else if ((i >= NumY - NumCtr_Coat) && ((j >= (NumCtrX - NumCtrX_Web) / 2) && (j < (NumCtrX + NumCtrX_Web) / 2)))
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);
			else
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial1, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}

		}
	}

	//Generating Elements for Web
	int JunctionLTag = OriginNodeTag + (NumX + 1) * (NumY + 1) - (NumX + 1) + (NumX - NumXweb) / 2;
	int JunctionUTag = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1) + (NumX - NumXweb) / 2;

	for (int i = 0; i < NumYweb; i++) {
		for (int j = 0; j < NumXweb; j++) {
			EleTag = OriginEleTag + NumXweb * i + j + NumX * NumY;

			if (i == 0)
			{
				NodeTag1 = JunctionLTag + j;
				NodeTag2 = JunctionLTag + j + 1;
				NodeTag3 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * i + j + 1;
				NodeTag4 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * i + j;
			}
			else if (i == NumYweb - 1)
			{
				NodeTag1 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (i - 1) + j;
				NodeTag2 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (i - 1) + j + 1;
				NodeTag3 = JunctionUTag + j + 1;
				NodeTag4 = JunctionUTag + j;

			}
			else
			{
				NodeTag1 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (i - 1) + j;
				NodeTag2 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (i - 1) + j + 1;
				NodeTag3 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * i + j + 1;
				NodeTag4 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * i + j;

			}


			
			if ((j >= NumCtr_Coat) && (j < NumXweb - NumCtr_Coat))
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);
			else
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial1, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}


		}
	}

	//Generating Elements for UpperFlange

	for (int i = 0; i < NumY; i++) {
		for (int j = 0; j < NumX; j++) {
			EleTag = OriginEleTag + NumX * i + j + NumX * NumY + NumXweb * NumYweb;
			NodeTag1 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1) + (NumX + 1) * i + j;
			NodeTag2 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1) + (NumX + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1) + (NumX + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + (NumX + 1) * (NumY + 1) + (NumXweb + 1) * (NumYweb - 1) + (NumX + 1) * (i + 1) + j;


			

			if ((i < NumCtr_Coat) && (j >= (NumCtrX - NumCtrX_Web) / 2) && (j < (NumCtrX + NumCtrX_Web) / 2))
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);
			else if (i >= NumCtr_Coat && i < NumY - NumCtr_Coat)
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);
			else
				TempEle = new QuadFour(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial1, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}

		}
	}
	//end of for(int i=0;i<NumY; i++);

	return 0;
}