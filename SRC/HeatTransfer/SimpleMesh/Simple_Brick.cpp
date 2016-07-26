#include <Simple_Brick.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>

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
Simple_Brick::InitialMeshCtrl(Vector& MeshCtrls)
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


