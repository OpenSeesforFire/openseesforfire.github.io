#include <Simple_Block.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>

Simple_Block::Simple_Block(int tag, double a1, double b1, double a3, double b3,double a2, double b2, double a4, double b4)
:Simple_Entity(tag,0),a1(a1), b1(b1), a3(a3), b3(b3), a2(a2), b2(b2), a4(a4), b4(b4),NumCtrlID(2),Seeds1(1),Seeds2(1),Seeds3(1),Seeds4(1)
{

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


int Simple_Block::InitialMeshCtrl(Vector& MeshCtrls)
{
    NumCtrlID(0)=(a2-a1)/MeshCtrls(0)+0.5;
	NumCtrlID(1)=(b3-b1)/MeshCtrls(1)+0.5;
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
		//opserr<<" Refined Seeds1: "<<Seeds1;
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


