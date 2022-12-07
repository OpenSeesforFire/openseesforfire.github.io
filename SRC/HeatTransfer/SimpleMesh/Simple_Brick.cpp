#include <Simple_Brick.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <BrickEight.h>

Simple_Brick::Simple_Brick(int tag, double CenterX, double CenterY, double CenterZ, double Breadth_X, double Height_Y,double Length_Z)
:Simple_Entity(tag,2),Breadth(Breadth_X), Height(Height_Y),Length(Length_Z),NumCtrlID(3),Seeds1(1),Seeds2(1),Seeds3(1)
{
	a1 = CenterX-Breadth_X/2;
	b1 = CenterY-Height_Y/2;
	c1 = CenterZ-Length_Z/2;
	a7 = CenterX+Breadth_X/2;
	b7 = CenterY+Height_Y/2;
	c7 = CenterZ+Length_Z/2;
}


Simple_Brick::~Simple_Brick()
{
}


int 
Simple_Brick::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
    double CtrX= MeshCtrls(0);
	double CtrY= MeshCtrls(1);
	double CtrZ= MeshCtrls(2);
    NumCtrlID(0)=Breadth/CtrX+0.5;
	NumCtrlID(1)=Height/CtrY+0.5;
	NumCtrlID(2)=Length/CtrZ+0.5;
	return 0;
}
bool 
Simple_Brick::InitialSeeds(void)
{
    bool result=true;
	Seeds1.resize(NumCtrlID(0)+1);
	Seeds2.resize(NumCtrlID(1)+1);
	Seeds3.resize(NumCtrlID(2)+1);

	for(int i =0; i<=NumCtrlID(0); i++){
		Seeds1(i)= a1+((a7-a1)/NumCtrlID(0))*i;
	}
	for(int i =0; i<=NumCtrlID(1); i++){
		Seeds2(i)= b1+((b7-b1)/NumCtrlID(1))*i;
	}
	for(int i =0; i<=NumCtrlID(2); i++){
		Seeds3(i)= c1+((c7-c1)/NumCtrlID(2))*i;
	}
	
	return result;
}


//Member function to adjust the seed distribution
int 
Simple_Brick::RefineSeeds(int SeedTag, const Vector& RefinedSeedsInfo)
{
	int result= 0;
	int RefinedInfoNum = (RefinedSeedsInfo.Size())/2;
	int RefinedNumCtrl=0;

	for (int i=0; i<RefinedInfoNum; i++){
			 int RecievedNum = RefinedSeedsInfo(2*i+1);
			 RefinedNumCtrl +=RecievedNum;
	}

	if(SeedTag == 1) {
		NumCtrlID(0)=RefinedNumCtrl;
		Seeds1.Zero();
		Seeds1.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds1(0)=a1;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds1(count+j+1)=a1+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd+=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds1(RefinedNumCtrl)-a7)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds1"<<endln;
			result=-1;
		}
		opserr<<" Refined Seeds1: "<<Seeds1;
	}
	else if(SeedTag ==2) {
		NumCtrlID(1)=RefinedNumCtrl;
		Seeds2.Zero();
		Seeds2.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds2(0)=b1;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds2(count+j+1)=b1+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count +=(int) RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds2(RefinedNumCtrl)-b7)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds3"<<endln;
			result=-1;
			}
	}
	else if(SeedTag ==3) {
		NumCtrlID(0)=RefinedNumCtrl;
		Seeds3.Zero();
		Seeds3.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds3(0)=c1;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds3(count+j+1)=c1+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count +=  (int)RefinedSeedsInfo(2*i+1);
		 CountCrd+=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds3(RefinedNumCtrl)-c7)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds3"<<endln;
			result=-1;
			}
	}
	else {
		opserr<<"WARNING! Simple_Block::RefineSeeds recieved an invalid Seed Tag"<<endln;
		result = -1;
	}
		
	return result;
}


const Vector& 
Simple_Brick::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ return Seeds1;}
	else if(SeedTag==2)
		{ return Seeds2;}
	else if(SeedTag ==3)
		{ return Seeds3;}
	else 
		{opserr<<"Warning: invalid Seedtag for Brick"<<endln;}

}

const ID&
Simple_Brick::GetNumCtrlID(void)
{
	return NumCtrlID;
}


int Simple_Brick::GetNumofNodes(void)
{
	return (NumCtrlID(0)+1)*(NumCtrlID(1)+1)*(NumCtrlID(2)+1);
}

int Simple_Brick::GetNumofEles(void)
{
	return NumCtrlID(0)*NumCtrlID(1)*NumCtrlID(2);
}

int Simple_Brick::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs)
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
	int NumCtrZ = this->GetNumCtrlID()(2);
	int NodeTag;

	for (int k = 0; k <= NumCtrZ; k++) {
		for (int j = 0; j <= NumCtrY; j++) {
			for (int i = 0; i <= NumCtrX; i++) {
				double NodeCrdX = 0;
				double NodeCrdY = 0;
				double NodeCrdZ = 0;


				NodeCrdX = (this->GetSeeds(1))(i);
				NodeCrdY = (this->GetSeeds(2))(j);
				NodeCrdZ = (this->GetSeeds(3))(k);

				if (fabs(NodeCrdX) < 1e-10)
					NodeCrdX = 0;

				if (fabs(NodeCrdY) < 1e-10)
					NodeCrdY = 0;

				if (fabs(NodeCrdZ) < 1e-10)
					NodeCrdZ = 0;



				NodeTag = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * k + (NumCtrX + 1) * j + i;
				HeatTransferNode* TempNode = new HeatTransferNode(NodeTag, nDoF, NodeCrdX, NodeCrdY, NodeCrdZ);
				if (theHTDomain->addNode(TempNode) < 0)
					return -1;
			}
		}


	}

#ifdef _DEBUG
	opserr << "Seed 1" << this->GetSeeds(1) << endln << "Seed 2" << this->GetSeeds(2) << endln << "Seed 3" << this->GetSeeds(3) << endln;
#endif
	return 0;
}



int Simple_Brick::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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

	int NumCtrX = this->GetNumCtrlID()(0);
	int NumCtrY = this->GetNumCtrlID()(1);
	int NumCtrZ = this->GetNumCtrlID()(2);
	int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8;

	for (int k = 0; k < NumCtrZ; k++) {
		for (int j = 0; j < NumCtrY; j++) {
			for (int i = 0; i < NumCtrX; i++) {
				EleTag = OriginEleTag + NumCtrX * NumCtrY * k + NumCtrX * j + i;
				NodeTag1 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * k + (NumCtrX + 1) * j + i;
				NodeTag2 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * k + (NumCtrX + 1) * j + i + 1;
				NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * k + (NumCtrX + 1) * (j + 1) + i + 1;
				NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * k + (NumCtrX + 1) * (j + 1) + i;
				NodeTag5 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * (k + 1) + (NumCtrX + 1) * j + i;
				NodeTag6 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * (k + 1) + (NumCtrX + 1) * j + i + 1;
				NodeTag7 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * (k + 1) + (NumCtrX + 1) * (j + 1) + i + 1;
				NodeTag8 = OriginNodeTag + (NumCtrX + 1) * (NumCtrY + 1) * (k + 1) + (NumCtrX + 1) * (j + 1) + i;

				
				HeatTransferElement* TempEle = 0;

				//TempEle = new BrickEight  (EleTag,NodeTag3,NodeTag4,NodeTag1,NodeTag2,NodeTag7,NodeTag8,NodeTag5,NodeTag6,*theHTMaterial, PhaseTransformation);
				TempEle = new BrickEight(EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8, *theHTMaterial, PhaseTransformation);

				if (theHTDomain->addElement(TempEle) < 0) {
					opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
					return -1;
				}
				
			}
		}
	}

	return 0;
}