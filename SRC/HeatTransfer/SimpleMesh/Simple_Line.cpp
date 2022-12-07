#include <Simple_Line.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <LineTwo.h>


Simple_Line::Simple_Line(int tag, double centerX, double lengthX)
:Simple_Entity(tag,5),CenterX(centerX), LengthX(lengthX),NumCtrlID(1),Seeds1(1)
{
	a1 = CenterX-LengthX/2;
	a2 = CenterX+LengthX/2;
}



Simple_Line::~Simple_Line()
{
}


int Simple_Line::InitialMeshCtrl(Vector& MeshCtrls, bool numCtrl)
{
    NumCtrlID(0)=(a2-a1)/MeshCtrls(0)+0.5;
	return 0;
}
bool Simple_Line::InitialSeeds(void)
{
    bool result=true;
	Seeds1.resize(NumCtrlID(0)+1);
	

	for(int i =0; i<=NumCtrlID(0); i++){
		Seeds1(i)= a1+((a2-a1)/NumCtrlID(0))*i;
		
	}
	
	return result;
}

/*
//Member function to adjust the seed distribution
bool Simple_Line::RefineSeeds(int SeedTag, const Vector& RefinedSeedsInfo)
{
	bool result= true;
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
			opserr<<"WARNING! Simple_Line::RefineSeeds failed to refine Seeds1"<<endln;
			result=false;
		}
		//opserr<<" Refined Seeds1: "<<Seeds1;
	}
	else if(SeedTag ==2) {
		NumCtrlID(1)=RefinedNumCtrl;
		Seeds2.Zero();
		Seeds2.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds2(0)=b2;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds2(count+j+1)=b2+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds2(RefinedNumCtrl)-b3)>1e-6){
			opserr<<"WARNING! Simple_Line::RefineSeeds failed to refine Seeds3"<<endln;
			result=false;
			}
	}
	else if(SeedTag ==3) {
		NumCtrlID(0)=RefinedNumCtrl;
		Seeds3.Zero();
		Seeds3.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds3(0)=a4;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds3(count+j+1)=a4+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds3(RefinedNumCtrl)-a3)>1e-6){
			opserr<<"WARNING! Simple_Line::RefineSeeds failed to refine Seeds3"<<endln;
			result=false;
		}
	}
	else if(SeedTag ==4) {
		NumCtrlID(1)=RefinedNumCtrl;
		Seeds4.Zero();
		Seeds4.resize(RefinedNumCtrl+1);
		int count=0; double CountCrd=0;
		Seeds4(0)=b1;
	    for (int i=0; i<RefinedInfoNum; i++){
		 for( int j=0; j<RefinedSeedsInfo(2*i+1);j++){
			Seeds4(count+j+1)=b1+CountCrd+ (RefinedSeedsInfo(2*i))*(j+1);
		 }
		 count += (int)RefinedSeedsInfo(2*i+1);
		 CountCrd +=(RefinedSeedsInfo(2*i+1)*(RefinedSeedsInfo(2*i)));
	   }
		if(fabs(Seeds4(RefinedNumCtrl)-b4)>1e-6){
			opserr<<"WARNING! Simple_Line::RefineSeeds failed to refine Seeds3"<<endln;
			result=false;
			}
	}
	else {
		opserr<<"WARNING! Simple_Line::RefineSeeds recieved an invalid Seed Tag"<<endln;
		result = false;
	}
		
	return result;
}
*/


const Vector& 
Simple_Line::GetSeeds(int SeedTag)
{
	if(SeedTag==1)
		{ return Seeds1;}
	else 
		{opserr<<"Warning: invalid facetag for Block"<<endln;}

}

const ID&
Simple_Line::GetNumCtrlID(void)
{
	return NumCtrlID;
}


int Simple_Line::GetNumofNodes(void)
{
	return (NumCtrlID(0)+1);
}

int Simple_Line::GetNumofEles(void)
{
	return NumCtrlID(0);
}


int Simple_Line::GenerateNodes(HeatTransferDomain* theHTDomain, int nDoF, const Vector& OriginLocs)
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
	for (int i = 0; i <= NumCtrX; i++) {
		double NodeCrdX = 0;
		NodeCrdX = (this->GetSeeds(1))(i);
		if (fabs(NodeCrdX) < 1e-10)
			NodeCrdX = 0;

		HeatTransferNode* TempNode = 0;
		if (OriginLocs.Size() == 1) {
			TempNode = new HeatTransferNode(OriginNodeTag + i, nDoF, OriginLoc1, NodeCrdX);
		}
		else if (OriginLocs.Size() == 2) {
			TempNode = new HeatTransferNode(OriginNodeTag + i, nDoF, OriginLoc1, NodeCrdX, OriginLoc2);
		}
		else {
			TempNode = new HeatTransferNode(OriginNodeTag + i, nDoF, NodeCrdX);
		}

		if (theHTDomain->addNode(TempNode) < 0)
			return -1;
		

	}

	return 0;

}


int Simple_Line::GenerateEles(HeatTransferDomain* theHTDomain, const ID& EleParameters, HeatTransferMaterial* theHTMaterial, HeatTransferMaterial* theHTMaterial1)
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
	HeatTransferElement* TempEle =0;

	int NumCtrX = this->GetNumCtrlID()(0);
	int EleTag, NodeTag1, NodeTag2;
	for (int i = 0; i < NumCtrX; i++) {
		EleTag = i;
		NodeTag1 = OriginNodeTag + i;
		NodeTag2 = OriginNodeTag + i + 1;

		
		TempEle = new LineTwo(OriginEleTag + EleTag, NodeTag1, NodeTag2, *theHTMaterial, PhaseTransformation);

		if (theHTDomain->addElement(TempEle) < 0) {
			opserr << "HeatTransferDomain failed to add element" << OriginEleTag + EleTag << endln;
			return -1;
		}
	}

	return 0;
}