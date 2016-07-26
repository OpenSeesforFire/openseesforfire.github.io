#include <Simple_Composite3D.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>

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
Simple_Composite3D::InitialMeshCtrl(Vector& MeshCtrls)
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



