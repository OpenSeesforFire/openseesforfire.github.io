#include <Simple_Block.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <QuadFour.h>

Simple_Block::Simple_Block(int tag, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4)
:Simple_Entity(tag,0),a1(a1), b1(b1),  a2(a2), b2(b2), a3(a3), b3(b3), a4(a4), b4(b4),NumCtrlID(2),Seeds1(1),Seeds2(1),Seeds3(1),Seeds4(1)
{
	opserr << "Define block using corner points: " << a1 << ", " << b1 << ",...,  " << a4 << ", " << b4 << endln;
}
Simple_Block::Simple_Block(int tag, double centerX, double centerY, double breadthX, double heightY)
:Simple_Entity(tag,0),CenterX(centerX), CenterY(centerY), BreadthX(breadthX), HeightY(heightY),NumCtrlID(2),Seeds1(1),Seeds2(1),Seeds3(1),Seeds4(1)
{
	a1 = a4 = CenterX-BreadthX/2;
	a2 = a3 = CenterX+BreadthX/2;
	b1 = b2 = CenterY-HeightY/2;
	b3 = b4 = CenterY+ HeightY/2;
}



Simple_Block::~Simple_Block()
{
}


int Simple_Block::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
	if (numCtrl) {
		NumCtrlID(0) = MeshCtrls(0);
		NumCtrlID(1) = MeshCtrls(1);
	}
	else {
		NumCtrlID(0) = (a2 - a1) / MeshCtrls(0) + 0.5;
		NumCtrlID(1) = (b3 - b1) / MeshCtrls(1) + 0.5;
	}
    
	return 0;
}
bool Simple_Block::InitialSeeds(void)
{
    bool result=true;
	Seeds1.resize(NumCtrlID(0)+1);
	Seeds4.resize(NumCtrlID(0)+1);
	Seeds2.resize(NumCtrlID(1)+1);
	Seeds3.resize(NumCtrlID(1)+1);

	for(int i =0; i<=NumCtrlID(0); i++){
		Seeds1(i)= a1+((a2-a1)/NumCtrlID(0))*i;
		Seeds4(i)= a4+((a3-a4)/NumCtrlID(0))*i;
	}
	for(int i =0; i<=NumCtrlID(1); i++){
		Seeds3(i)= b2+((b3-b2)/NumCtrlID(1))*i;
		Seeds2(i)= b1+((b4-b1)/NumCtrlID(1))*i;
	}

	return result;
}

//Member function to adjust the seed distribution
int Simple_Block::RefineSeeds(int SeedTag, const Vector& RefinedSeedsInfo)
{
	int result= 0;
	int RefinedInfoNum = (RefinedSeedsInfo.Size())/2;
	int RefinedNumCtrl=0;

	for (int i=0; i<RefinedInfoNum; i++){
			RefinedNumCtrl += RefinedSeedsInfo(2*i+1);
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
		if(fabs(Seeds1(RefinedNumCtrl)-a2)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds1"<<endln;
			result=-1;
		}
		opserr<<" Refined Seeds1: "<<Seeds1;
	}
	else if(SeedTag ==3) {
		NumCtrlID(1)=RefinedNumCtrl;
		Seeds3.Zero();
		Seeds3.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds3(0)=b2;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds3(count+j+1)=b2+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds3(RefinedNumCtrl)-b3)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds3"<<endln;
			result=-1;
			}
	}
	else if(SeedTag ==4) {
		NumCtrlID(0)=RefinedNumCtrl;
		Seeds4.Zero();
		Seeds4.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds4(0)=a4;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds4(count+j+1)=a4+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds4(RefinedNumCtrl)-a3)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds4"<<endln;
			result=-1;
		}
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
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds2(RefinedNumCtrl)-b4)>1e-6){
			opserr<<"WARNING! Simple_Block::RefineSeeds failed to refine Seeds2"<<endln;
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
Simple_Block::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ return Seeds1;}
	else if(SeedTag==2)
		{ return Seeds2;}
	else if(SeedTag ==3)
		{ return Seeds3;}
    else if(SeedTag ==4)
		{ return Seeds4;}
	else 
		{opserr<<"Warning: invalid facetag for Block"<<endln;}

}

const ID&
Simple_Block::GetNumCtrlID(void)
{
	return NumCtrlID;
}


int Simple_Block::GetNumofNodes(void)
{
	return (NumCtrlID(0)+1)*(NumCtrlID(1)+1);
}

int Simple_Block::GetNumofEles(void)
{
	return NumCtrlID(0)*NumCtrlID(1);
}

int Simple_Block::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs)
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

	for (int i = 0; i <= NumCtrY; i++) {
		for (int j = 0; j <= NumCtrX; j++) {
			double NodeCrdX = 0;
			double NodeCrdY = 0;
			if (this->GetSeeds(1)(j) == this->GetSeeds(4)(j))
				NodeCrdX = this->GetSeeds(1)(j);
			else {
				double crdX1i = this->GetSeeds(1)(j);
				double crdX2i = this->GetSeeds(4)(j);
				NodeCrdX = crdX1i + i / NumCtrY * (crdX2i - crdX1i);
			}

			if (fabs(NodeCrdX) < 1e-10)
				NodeCrdX = 0;

			if ((this->GetSeeds(2))(i) == (this->GetSeeds(3))(i))
				NodeCrdY = (this->GetSeeds(2))(i);
			else {
				double crdY1i = this->GetSeeds(2)(j);
				double crdY2i = this->GetSeeds(3)(j);
				NodeCrdY = crdY1i + j / NumCtrX * (crdY2i - crdY1i);
			}


			if (fabs(NodeCrdY) < 1e-10)
				NodeCrdY = 0;

			HeatTransferNode* TempNode = 0;
			if (OriginLocs.Size() == 1) {
				TempNode = new HeatTransferNode(OriginNodeTag + (NumCtrX + 1) * i + j, nDoF, NodeCrdX, NodeCrdY, OriginLoc1);
			}
			else {
				TempNode = new HeatTransferNode(OriginNodeTag + (NumCtrX + 1) * i + j, nDoF, NodeCrdX, NodeCrdY);
			}
			if (theHTDomain->addNode(TempNode) < 0) {
				opserr << "HTDomain failed to generate node with coordinates: " << NodeCrdX << ", " << NodeCrdY << endln;
				return -1;
			}
			

		}
	}


}



int Simple_Block::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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

	int OriginNodeTag = theHTDomain->getNumNodes() - (this->GetNumofNodes())+1;
	int OriginEleTag = theHTDomain->getNumElements() + 1;

	int NumCtrX = this->GetNumCtrlID()(0);
	int NumCtrY = this->GetNumCtrlID()(1);
	int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;
	for (int i = 0; i < NumCtrY; i++) {
		for (int j = 0; j < NumCtrX; j++) {
			EleTag = NumCtrX * i + j;
			NodeTag1 = OriginNodeTag + (NumCtrX + 1) * i + j;
			NodeTag2 = OriginNodeTag + (NumCtrX + 1) * i + j + 1;
			NodeTag3 = OriginNodeTag + (NumCtrX + 1) * (i + 1) + j + 1;
			NodeTag4 = OriginNodeTag + (NumCtrX + 1) * (i + 1) + j;
			HeatTransferElement* TempEle = 0;
			TempEle = new QuadFour(OriginEleTag + EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, *theHTMaterial, PhaseTransformation);

			if (theHTDomain->addElement(TempEle) < 0) {
				opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
				return -1;
			}
			
		}
	}
	return 0;
}