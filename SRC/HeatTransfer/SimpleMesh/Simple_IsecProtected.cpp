//This code is based on OpenSees Framework 
//Written by Liming Jiang(University of Edinburgh)

#include <Simple_IsecProtected.h>
#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>


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
Simple_IsecProtected::InitialMeshCtrl(Vector& MeshCtrls)
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



