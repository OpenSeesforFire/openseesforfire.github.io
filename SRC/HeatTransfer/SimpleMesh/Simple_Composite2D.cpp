#include <Simple_Composite2D.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>

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
Simple_Composite2D::InitialMeshCtrl(Vector& MeshCtrls)
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



