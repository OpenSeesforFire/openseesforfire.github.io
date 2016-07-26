//This code is based on OpenSees Framework 
//Written by Liming Jiang(University of Edinburgh)

#include <stdlib.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h>
#include <string.h>

#include <Simple_Mesh.h>
#include <StandardStream.h>




//StandardStream sserr;
//OPS_Stream* opserrPtr = &sserr;



Simple_Mesh::Simple_Mesh(int tag, Simple_Entity* Entity,HeatTransferDomain* theDomain,HeatTransferMaterial* theHTMaterial, Vector& MeshCtrls,
						 HeatTransferMaterial* theHTMaterial1 )
:TaggedObject(tag),theHTDomain(theDomain),isHTDomain(true), theHTMaterial(theHTMaterial),OriginLocs(0),EleParameters(0),
 theHTMaterial1(theHTMaterial1)
{
	 theEntity = Entity;
   theEntity->setMeshTag(this->getTag());
	 theEntity->InitialMeshCtrl(MeshCtrls);

	 if (theEntity->InitialSeeds()==false)
		 opserr<<"Mesh "<<this->getTag()<< "failed to Intial the mesh Seeds"<<endln;
}


//Simple_Mesh::Simple_Mesh(tag){

//}



Simple_Mesh::~Simple_Mesh()
{

}


int Simple_Mesh::GetNumOfNodes()
{
	int NumNodes= theEntity->GetNumofNodes();
	if(NumNodes==0){
		opserr<< "0 Node in Mesh"<<this->getTag()<<endln;
	}
	return NumNodes;
}

int Simple_Mesh::GetNumOfEles()
{
	int NumEles= theEntity->GetNumofEles();	
   if(NumEles==0){
		opserr<< "0 Element in Mesh"<<this->getTag()<<endln;
	}
	
	return NumEles;
}


int Simple_Mesh::SetOriginLocs(const Vector& originLocs)
{

  if (originLocs!=0||originLocs.Size()!=0) {
    OriginLocs.Zero();
    OriginLocs = originLocs;
  }
  return 0;
}


int Simple_Mesh::SetEleParameters(const ID& eleParameters)
{
  if (eleParameters!=0) {
    EleParameters.Zero();
    EleParameters = eleParameters;
  }
  return 0;
}

bool Simple_Mesh::GeneratingNodes(const Vector& originLocs )
{
  if(OriginLocs==0){
    if (originLocs!=0) {
      OriginLocs.Zero();
      OriginLocs = originLocs;
   }
  }
  else
  {
    opserr<<"SimpleMesh::GeneratingNodes encountered a redefined OriginLocs"<<endln;
  }
  
  
  bool result=true;
	
	int nDoF=1;
#ifdef _BEBUG
	opserr<<"SimpleMesh "<<this->getTag()<<" starts to generate nodes.."<<endln;
#endif
    if(isHTDomain)
		OriginNodeTag = theHTDomain->getNumNodes()+1;
	else
		OriginNodeTag = theDomain->getNumNodes()+1;

   double OriginLoc1,OriginLoc2;
	if(OriginLocs.Size()==1){
		OriginLoc1=OriginLocs(0);
		}
	else if(OriginLocs.Size()==2) {
		OriginLoc1=OriginLocs(0);
		OriginLoc2=OriginLocs(1);
		}

	//------------------Gnenerating Nodes for 1D line---------------------------//
	if((theEntity->getEntityTypeTag())==5)
	{
		int NumCtrX = theEntity->GetNumCtrlID()(0);
		for(int i=0 ;i<= NumCtrX;i++){
			double NodeCrdX =0;
			NodeCrdX = (theEntity->GetSeeds(1))(i);
			if (fabs(NodeCrdX)<1e-10)
					   NodeCrdX=0;

			if(isHTDomain){
				HeatTransferNode* TempNode=0;
				if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode(OriginNodeTag+i,nDoF,OriginLoc1,NodeCrdX);
					}
				else if(OriginLocs.Size()==2) {
					TempNode = new HeatTransferNode(OriginNodeTag+i,nDoF, OriginLoc1,OriginLoc2,NodeCrdX);
					}
				else{
					TempNode = new HeatTransferNode(OriginNodeTag+i,nDoF,NodeCrdX);
					}

				result = theHTDomain->addNode(TempNode);
			}
			else{
				opserr<<"Simple_Mesh::HeatTranferDomain doesn't support LineTwo element";
				return result =false;
			}
		}


	}
	else if((theEntity->getEntityTypeTag())==0)
	{
		//Generating Nodes for Block	
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);

		for(int i=0; i <= NumCtrY;i++) {
			for(int j=0; j<=NumCtrX; j++) {
				double NodeCrdX=0;
                double NodeCrdY=0;
				if(theEntity->GetSeeds(1)(j)== theEntity->GetSeeds(4)(j))
					NodeCrdX = theEntity->GetSeeds(1)(j);
				else
					opserr<<"Block should be a rectangular one so far"<<endln;
                 
				if (fabs(NodeCrdX)<1e-10)
					   NodeCrdX=0;
	
				if((theEntity->GetSeeds(2))(i)== (theEntity->GetSeeds(3))(i))
				    NodeCrdY =	(theEntity->GetSeeds(2))(i);
				else
					opserr<<"Block should be a rectangular one so far"<<endln;

			
				if (fabs(NodeCrdY)<1e-10)
					   NodeCrdY=0;
                
				if(isHTDomain)
				{
					
					HeatTransferNode* TempNode=0;
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode(OriginNodeTag+(NumCtrX+1)*i+j,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode(OriginNodeTag+(NumCtrX+1)*i+j,nDoF,NodeCrdX,NodeCrdY);
					}
					result = theHTDomain->addNode(TempNode);
				}
				else
				{
					Node* TempNode = new Node(OriginNodeTag+(NumCtrX+1)*i+j,nDoF,NodeCrdX,NodeCrdY);
					result = theDomain->addNode(TempNode);
				}

			}
	    } 	
		
	}
	//------------------Gnenerating Nodes for 3D Brick---------------------------//
	else if((theEntity->getEntityTypeTag())==2)
	{
		//Generating Nodes for Brick	
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrZ= theEntity->GetNumCtrlID()(2);
		int NodeTag;

		for(int k=0; k<= NumCtrZ; k++) {
			for(int j=0; j <= NumCtrY;j++) {
				for(int i=0; i<=NumCtrX; i++) {
					double NodeCrdX=0;
	                double NodeCrdY=0;
					double NodeCrdZ=0;

					
					NodeCrdX = (theEntity->GetSeeds(1))(i);
					NodeCrdY =	(theEntity->GetSeeds(2))(j);
					NodeCrdZ =	(theEntity->GetSeeds(3))(k);
				
					if (fabs(NodeCrdX)<1e-10)
						   NodeCrdX=0;
				
					if (fabs(NodeCrdY)<1e-10)
						   NodeCrdY=0;
					
					if (fabs(NodeCrdZ)<1e-10)
						   NodeCrdZ=0;
					
                
					if(isHTDomain)
					{
						NodeTag = OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*k+(NumCtrX+1)*j+i;
						HeatTransferNode* TempNode = new HeatTransferNode(NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
						result = theHTDomain->addNode(TempNode);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in brick"<<endln;
					}

				}
		    }
			
			
		} 	
		
	}
	//----------------Gnenerating Nodes for 2D I-Beam Section-------------------------//
	else if((theEntity->getEntityTypeTag())==1)
	{
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
        double NodeCrdX, NodeCrdY;
	    //Generating Nodes for LowerFlange
		for(int i=0; i <= NumCtrY;i++) {
			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = (NumCtrX+1)*i+j ;
				NodeCrdX = (theEntity->GetSeeds(1))(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = (theEntity->GetSeeds(2))(i);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				if(isHTDomain)
				{
					HeatTransferNode* TempNode=0;
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}

					result = theHTDomain->addNode(TempNode);
				}
				else
				{
					Node* TempNode = new Node (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					result = theDomain->addNode(TempNode);
				}
					
			}
		} 
			
		//Generating Nodes for Web
		for(int i=0; i < NumCtrY_Web-1;i++) {
			for(int j=0; j<=NumCtrX_Web; j++) {
				
				int NodeTag = (NumCtrX_Web+1)*i+j+(NumCtrX+1)*(NumCtrY+1) ;
				NodeCrdX = (theEntity->GetSeeds(1))(j+(NumCtrX-NumCtrX_Web)/2);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+1);	
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				
				if(isHTDomain)
				{
					HeatTransferNode* TempNode=0;
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}

					result = theHTDomain->addNode(TempNode);
				}
				else
				{
					Node* TempNode = new Node (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					result = theDomain->addNode(TempNode);
				}
			}
		} 
		
		//Generating Nodes for UpperFlange
		for(int i=0; i <= NumCtrY;i++) {
			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = (NumCtrX+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
				NodeCrdX = (theEntity->GetSeeds(1))(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+NumCtrY_Web);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				if(isHTDomain)
				{
					HeatTransferNode* TempNode=0;
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}
					result = theHTDomain->addNode(TempNode);
				}
				else
				{
					Node* TempNode = new Node (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
					result = theDomain->addNode(TempNode);			
				}
			}
		} 
		
		
		
	}
	//----------------Gnenerating Nodes for 3D I-Beam Section-------------------------//
	else if((theEntity->getEntityTypeTag())==3)
	{
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		//Cordinates declared for nodes being generated here
        double NodeCrdX, NodeCrdY,NodeCrdZ;
		//for loop with k along the Beam length
		for(int k=0; k<=NumCtrZ; k++){
			// Nodal Cordinate y
			NodeCrdZ = (theEntity->GetSeeds(3))(k);
				if (fabs(NodeCrdZ)<1e-10)
					NodeCrdZ=0;
	    //Generating Nodes for LowerFlange
		for(int i=0; i <= NumCtrY;i++) {

			NodeCrdY = (theEntity->GetSeeds(2))(i);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = OriginNodeTag+ (2*(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1))*k+ (NumCtrX+1)*i+j ;
				
				NodeCrdX = (theEntity->GetSeeds(1))(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				if(isHTDomain)
				{
				    HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
					
			}
		} 
			
		//Generating Nodes for Web
		for(int i=0; i < NumCtrY_Web-1;i++) {
			// Nodal Cordinate y
			NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+1);	
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX_Web; j++) {
				
				int NodeTag = OriginNodeTag+(2*(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1))*k+ (NumCtrX_Web+1)*i+j+(NumCtrX+1)*(NumCtrY+1) ;
				// Nodal Cordinate x
				NodeCrdX = (theEntity->GetSeeds(1))(j+(NumCtrX-NumCtrX_Web)/2);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;
								
				if(isHTDomain)
				{
					HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		} 
		
		//Generating Nodes for UpperFlange
		for(int i=0; i <= NumCtrY;i++) {
			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = OriginNodeTag+(2*(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1))*k+ (NumCtrX+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
				NodeCrdX = (theEntity->GetSeeds(1))(j);
			 	if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+NumCtrY_Web);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				if(isHTDomain)
				{
					HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		}
		
		} 
		//end of for(k=0;k<NumCtrZ;k++)
		
		
	}
	// end of if (entitytypetag ==3)
	
	//----------------Gnenerating Nodes for 3D CompositeBeam Section-------------------------//
	else if((theEntity->getEntityTypeTag())==7)
	{
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(5);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(6);
		
		//Cordinates declared for nodes being generated here
        double NodeCrdX, NodeCrdY,NodeCrdZ;
		//for loop with k along the Beam length
		for(int k=0; k<=NumCtrZ; k++){
			// Nodal Cordinate y
			NodeCrdZ = (theEntity->GetSeeds(3))(k);
			if (fabs(NodeCrdZ)<1e-10)
				NodeCrdZ=0;
			
			int LayerOriginNodeTag = OriginNodeTag+ ((NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+(NumCtrX_Slab+1)*(NumCtrY_Slab+1))*k;
				
			
	    //Generating Nodes for LowerFlange
		for(int i=0; i <= NumCtrY;i++) {

			NodeCrdY = (theEntity->GetSeeds(2))(i);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = LayerOriginNodeTag + (NumCtrX+1)*i+j ;
				
				NodeCrdX = (theEntity->GetSeeds(1))(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				if(isHTDomain)
				{
				    HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
					
			}
		} 
			
		//Generating Nodes for Web
		for(int i=0; i < NumCtrY_Web-1;i++) {
			// Nodal Cordinate y
			NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+1);	
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX_Web; j++) {
				
				int NodeTag = LayerOriginNodeTag+ (NumCtrX_Web+1)*i+j+(NumCtrX+1)*(NumCtrY+1) ;
				// Nodal Cordinate x
				NodeCrdX = (theEntity->GetSeeds(1))(j+(NumCtrX-NumCtrX_Web)/2);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;
								
				if(isHTDomain)
				{
					HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		} 
		
		//Generating Nodes for UpperFlange(Top surface nodes not generated)
		for(int i=0; i < NumCtrY;i++) {
			NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+NumCtrY_Web);
			if (fabs(NodeCrdY)<1e-10)
				NodeCrdY=0;
			
			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = LayerOriginNodeTag+ (NumCtrX+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
				NodeCrdX = (theEntity->GetSeeds(1))(j);
			 	if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;
				
				if(isHTDomain)
				{
					HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		}
		
		//Generating Nodes for Slab section
		for(int i=0; i <= NumCtrY_Slab;i++) {
			
			NodeCrdY = (theEntity->GetSeeds(5))(i);
			if (fabs(NodeCrdY)<1e-10)
				NodeCrdY=0;
			
			for(int j=0; j<=NumCtrX_Slab; j++) {
				int NodeTag = LayerOriginNodeTag+ (NumCtrX_Slab+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY;
				NodeCrdX = (theEntity->GetSeeds(4))(j);
			 	if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				if(isHTDomain)
				{
					HeatTransferNode* TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,NodeCrdZ);
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		}
		//End of adding nodes for Composite_Slab
		
		
		} 
		//end of for(k=0;k<NumCtrZ;k++)
		
		
	}
	// end of if (entitytypetag ==7)

	//
	//----------------Gnenerating Nodes for 2D CompositeBeam Section-------------------------//
	else if((theEntity->getEntityTypeTag())==6)
	{
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(4);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(5);
		
		//Cordinates declared for nodes being generated here
        double NodeCrdX, NodeCrdY;
		//for loop with k along the Beam length
		HeatTransferNode* TempNode=0;	
			
	    //Generating Nodes for LowerFlange
		for(int i=0; i <= NumCtrY;i++) {

			NodeCrdY = (theEntity->GetSeeds(2))(i);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = OriginNodeTag + (NumCtrX+1)*i+j ;
				
				NodeCrdX = (theEntity->GetSeeds(1))(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				if(isHTDomain)
				{
				    
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}
					
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
					
			}
		} 
			
		//Generating Nodes for Web
		for(int i=0; i < NumCtrY_Web-1;i++) {
			// Nodal Cordinate y
			NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+1);	
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;

			for(int j=0; j<=NumCtrX_Web; j++) {
				
				int NodeTag = OriginNodeTag+ (NumCtrX_Web+1)*i+j+(NumCtrX+1)*(NumCtrY+1) ;
				// Nodal Cordinate x
				NodeCrdX = (theEntity->GetSeeds(1))(j+(NumCtrX-NumCtrX_Web)/2);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;
								
				if(isHTDomain)
				{
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}
					
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		} 
		
		//Generating Nodes for UpperFlange(Top surface nodes not generated)
		for(int i=0; i < NumCtrY;i++) {
			NodeCrdY = (theEntity->GetSeeds(2))(i+NumCtrY+NumCtrY_Web);
			if (fabs(NodeCrdY)<1e-10)
				NodeCrdY=0;
			
			for(int j=0; j<=NumCtrX; j++) {
				int NodeTag = OriginNodeTag+ (NumCtrX+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
				NodeCrdX = (theEntity->GetSeeds(1))(j);
			 	if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;
				
				if(isHTDomain)
				{
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}
					
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		}
		
		//Generating Nodes for Slab section
		for(int i=0; i <= NumCtrY_Slab;i++) {
			
			NodeCrdY = (theEntity->GetSeeds(4))(i);
			if (fabs(NodeCrdY)<1e-10)
				NodeCrdY=0;
			
			for(int j=0; j<=NumCtrX_Slab; j++) {
				int NodeTag = OriginNodeTag+ (NumCtrX_Slab+1)*i+j+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY;
				NodeCrdX = (theEntity->GetSeeds(3))(j);
			 	if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				if(isHTDomain)
				{
					if(OriginLocs.Size()==1){
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
					}
					else{
					TempNode = new HeatTransferNode (NodeTag,nDoF,NodeCrdX,NodeCrdY);
					}
					
					result = theHTDomain->addNode(TempNode);
					//opserr<<"Adding Node "<<NodeTag <<" : "<<NodeCrdX<<" , "<<NodeCrdY<<" , " <<NodeCrdZ<<endln;
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
				}
			}
		}
		//End of adding nodes for Composite_Slab
		
		
	}
	// end of if (entitytypetag ==6)

	//----------------Gnenerating Nodes for 2D PROTECTED  Section-------------------------//
	else if((theEntity->getEntityTypeTag())==11)
	{
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtr_Coat = theEntity->GetNumCtrlID()(4);

		int NumX = NumCtrX;
		int NumY = NumCtrY+ NumCtr_Coat*2;
		int NumXweb = NumCtrX_Web + NumCtr_Coat*2;
		int NumYweb = NumCtrY_Web - NumCtr_Coat*2;

        double NodeCrdX, NodeCrdY;
		Vector Seeds1 = theEntity->GetSeeds(1);
		Vector Seeds2 = theEntity->GetSeeds(2);
	    //Generating Nodes for Lower great Flange
		for(int i=0; i <= NumY ;i++) {
			for(int j=0; j<=NumX; j++) {
				int NodeTag = (NumX+1)*i+j ;
				NodeCrdX = Seeds1(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = Seeds2(i);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				HeatTransferNode* TempNode=0;
				if(OriginLocs.Size()==1){
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
				}
				else{
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
				}

				result = theHTDomain->addNode(TempNode);
				
					
			}
		} 
			
		//Generating Nodes for Web
		for(int i=0; i < NumYweb -1;i++) {
			for(int j=0; j<=NumXweb; j++) {
				
				int NodeTag = (NumXweb+1)*i+j+(NumX+1)*(NumY+1) ;
				NodeCrdX = Seeds1(j+(NumX-NumXweb)/2);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = Seeds2(i+NumY+1);	
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
				
				HeatTransferNode* TempNode=0;
				if(OriginLocs.Size()==1){
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
				}
				else{
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
				}

				result = theHTDomain->addNode(TempNode);
				
			}
		} 
		
		//Generating Nodes for UpperFlange
		for(int i=0; i <= NumY;i++) {
			for(int j=0; j<=NumX; j++) {
				int NodeTag = (NumX+1)*i+j+ (NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1);
				NodeCrdX = Seeds1(j);
				if (fabs(NodeCrdX)<1e-10)
					NodeCrdX=0;

				NodeCrdY = Seeds2(i+NumY+NumYweb);
				if (fabs(NodeCrdY)<1e-10)
					NodeCrdY=0;
			
				HeatTransferNode* TempNode=0;
				if(OriginLocs.Size()==1){
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY,OriginLoc1);
				}
				else{
				TempNode = new HeatTransferNode (OriginNodeTag+NodeTag,nDoF,NodeCrdX,NodeCrdY);
				}
				result = theHTDomain->addNode(TempNode);
			}
		} 	
		
	}
	// end of if (entitytypetag ==11)  2d protected i section
	

		if(!result)
		opserr<<"Error:SimpleMesh::GeneratingNodes() for SimpleMesh " <<this->getTag()<<endln; 
	else {
#ifdef _DEBUG
		if(isHTDomain)
		opserr<<"SimpleMesh has successfully generated "<<(theHTDomain->getNumNodes())-OriginNodeTag+1 <<" HeatTransfer Nodes..."<<endln;
		else
		opserr<<"SimpleMesh has successfully generated "<<(theDomain->getNumNodes())-OriginNodeTag+1 <<" nodes..."<<endln;
#endif
	}
	
	return result;
}

//-----------------------Mesh tool for generating elements-----------------------------
int Simple_Mesh::GeneratingEles(const ID& eleParameters)
{
  if (EleParameters==0) {
    if (eleParameters!=0) {
      EleParameters.Zero();
      EleParameters=eleParameters;
    }
  }
  else
  {
    opserr<<"SimpleMesh::GeneratingEles encounters redefined EleParameteres"<<endln;
  }
  
  bool PhaseTransformation=false;
	bool PhaseTransformation1=false;
	
	if(EleParameters!=0){
	  if(EleParameters(0)==1){
		PhaseTransformation=true;
	   }
	   else {
		   if(theHTMaterial1!=0&&EleParameters.Size()>1){
			if(EleParameters(1)==1)
				PhaseTransformation1=true;
		   
		   }
	   }
	}
#ifdef _DEBUG
	opserr<<"SimpleMesh "<<this->getTag()<<" starts to generate elements...."<<endln;
#endif
	bool result=true;
    const char* type = "PlaneStress";    //This a temporary parameter for mesh-generating 
	
	if(isHTDomain)
		OriginEleTag = theHTDomain->getNumElements()+1;
	else
		OriginEleTag = theDomain->getNumElements()+1;
	
	if((theEntity->getEntityTypeTag())==5){
		//Generating Elements for Simple_Line
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int EleTag,NodeTag1,NodeTag2;
		for(int i=0; i<NumCtrX; i++) {
			EleTag = i;
			NodeTag1= OriginNodeTag+i;
			NodeTag2= OriginNodeTag+i+1; 

			if(isHTDomain)
			{
			  LineTwo* TempEle;	   
			  TempEle = new LineTwo(OriginEleTag+EleTag,NodeTag1,NodeTag2,*theHTMaterial, PhaseTransformation);
					
			  result=theHTDomain->addElement(TempEle);
			}
			else{		
			  opserr<<"Simple_Mesh::HeatTranferDomain doesn't support LineTwo element";
			}
		}
	}
	else if((theEntity->getEntityTypeTag())==0)
	{
		//Generating Elements for Block-
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4;
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = NumCtrX*i+j;
				NodeTag1= OriginNodeTag+(NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+(NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+(NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+(NumCtrX+1)*(i+1)+j;
			
				
				if(isHTDomain)
				{
				    QuadFour* TempEle;
				   
					TempEle = new QuadFour  (OriginEleTag+EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
					
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					FourNodeQuad* TempEle = new FourNodeQuad (OriginEleTag+EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theNDMaterial, type, EleParameters(0),EleParameters(1),EleParameters(2));
				    result= theDomain->addElement(TempEle);
					
				}
			}
		}
		
	}
	else if((theEntity->getEntityTypeTag())==2)
	{
		//Generating Elements for Block-
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrZ= theEntity->GetNumCtrlID()(2);
		int EleTag, NodeTag1,NodeTag2,NodeTag3,NodeTag4, NodeTag5, NodeTag6,NodeTag7, NodeTag8;

		for(int k=0;k<NumCtrZ; k++) {
			for(int j=0;j<NumCtrY; j++) {
				for(int i=0; i<NumCtrX; i++){
					EleTag = OriginEleTag+NumCtrX*NumCtrY*k+NumCtrX*j+i;
					NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*k+(NumCtrX+1)*j+i;
					NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*k+(NumCtrX+1)*j+i+1;
					NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*k+(NumCtrX+1)*(j+1)+i+1;
					NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*k+(NumCtrX+1)*(j+1)+i;
					NodeTag5= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*(k+1)+(NumCtrX+1)*j+i;
					NodeTag6= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*(k+1)+(NumCtrX+1)*j+i+1;
					NodeTag7= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*(k+1)+(NumCtrX+1)*(j+1)+i+1;
					NodeTag8= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*(k+1)+(NumCtrX+1)*(j+1)+i;
		
					if(isHTDomain)
					{
						BrickEight* TempEle=0;
						
						//TempEle = new BrickEight  (EleTag,NodeTag3,NodeTag4,NodeTag1,NodeTag2,NodeTag7,NodeTag8,NodeTag5,NodeTag6,*theHTMaterial, PhaseTransformation);
						TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);
						
						result=theHTDomain->addElement(TempEle);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in brick"<<endln;
						return(-1);
					}
				}
			}
		}
	}
	else if((theEntity->getEntityTypeTag())==1)
	{
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);

		//Generating Elements for LowerFlange
		int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag+NumCtrX*i+j;
				NodeTag1= OriginNodeTag+(NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+(NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+(NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+(NumCtrX+1)*(i+1)+j;
				if(isHTDomain)
				{
				    QuadFour* TempEle=0;
				    
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
					
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					FourNodeQuad* TempEle = new FourNodeQuad (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theNDMaterial, type, EleParameters(0),EleParameters(1),EleParameters(2));
				    result= theDomain->addElement(TempEle);
				}
			}
		}
		
		//Generating Elements for Web
		int JunctionLTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)-(NumCtrX+1)+(NumCtrX-NumCtrX_Web)/2;
        int JunctionUTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX-NumCtrX_Web)/2;

		for(int i=0;i<NumCtrY_Web; i++) {
			for(int j=0; j<NumCtrX_Web; j++){
				EleTag = OriginEleTag+NumCtrX_Web*i+j+NumCtrX*NumCtrY;
			
				if(i==0)
				{
					NodeTag1= JunctionLTag+j;
					NodeTag2= JunctionLTag+j+1;
					NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
				}
				else if(i==NumCtrY_Web-1)
				{
					NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
					NodeTag3= JunctionUTag+j+1;
					NodeTag4= JunctionUTag+j;
					
				}
				else
				{
					NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
					NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					
				}

				if(isHTDomain)
				{
				    QuadFour* TempEle=0;
				    
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
					
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					FourNodeQuad* TempEle = new FourNodeQuad (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theNDMaterial, type, EleParameters(0),EleParameters(1),EleParameters(2));
				    result= theDomain->addElement(TempEle);
				}
				
			}
		}
		
		//Generating Elements for UpperFlange
		
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag+NumCtrX*i+j+NumCtrX*NumCtrY+NumCtrX_Web*NumCtrY_Web;
				NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
			
				
				if(isHTDomain)
				{
					QuadFour* TempEle=0;
				    
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
		
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding elements in 2D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		//end of for(int i=0;i<NumCtrY; i++);
		
	}
	//end of else if((theEntity->getEntityTypeTag())==1)
	//---------------------------------------------------Meshing 3D I-section beam-----------------------------------------------
	else if((theEntity->getEntityTypeTag())==3)
	{
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		
		int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1)*2 +(NumCtrX_Web+1)*(NumCtrY_Web-1);

		int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8;
		for (int k = 0; k<NumCtrZ;k++){
			//along the length
		//Generating Elements for LowerFlange
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag+ (NumCtrX*NumCtrY*2 +NumCtrX_Web*NumCtrY_Web)*k + NumCtrX*i+j;
				NodeTag1= OriginNodeTag+ NumNodesPerLayer*k + (NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*(i+1)+j;
				NodeTag5= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*i+j;
				NodeTag6= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*i+j+1;
				NodeTag7= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*(i+1)+j+1;
				NodeTag8= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*(i+1)+j;
				
				
				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);
					
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding elements in 3D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		
		//Generating Elements for Web
		int JunctionLTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)-(NumCtrX+1)+(NumCtrX-NumCtrX_Web)/2;
        int JunctionUTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX-NumCtrX_Web)/2;

		for(int i=0;i<NumCtrY_Web; i++) {
			for(int j=0; j<NumCtrX_Web; j++){
				EleTag = OriginEleTag+(NumCtrX*NumCtrY*2+NumCtrX_Web*NumCtrY_Web)*k + NumCtrX_Web*i+j+NumCtrX*NumCtrY;
			
				if(i==0)
				{
					NodeTag1= JunctionLTag+ NumNodesPerLayer *k +j;
					NodeTag2= JunctionLTag+ NumNodesPerLayer *k +j+1;
					NodeTag3= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
					
					NodeTag5= JunctionLTag+NumNodesPerLayer *(k+1) +j;
					NodeTag6= JunctionLTag+NumNodesPerLayer *(k+1) +j+1;
					NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag8= OriginNodeTag+NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
					
				}
				
				else if(i==NumCtrY_Web-1)
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
					NodeTag3= JunctionUTag + NumNodesPerLayer *k +j+1;
					NodeTag4= JunctionUTag + NumNodesPerLayer *k +j;
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
					NodeTag7= JunctionUTag+ NumNodesPerLayer *(k+1) +j+1;
					NodeTag8= JunctionUTag+ NumNodesPerLayer *(k+1) +j;
					
				}
				else
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
					NodeTag3= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
					NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag8= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					
					
				}

				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);
				
					result=theHTDomain->addElement(TempEle);
				    
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
					return(-1);
				}
				
			}
		}
		
		//Generating Elements for UpperFlange
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag + (NumCtrX*NumCtrY*2+NumCtrX_Web*NumCtrY_Web)*k +NumCtrX*i+j+NumCtrX*NumCtrY+NumCtrX_Web*NumCtrY_Web;
				NodeTag1= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
				
				NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
				NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
				NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
				NodeTag8= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
				
			
				
				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);

 					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		//end of for(int i=0;i<NumCtrY; i++);
	}
	//end of for (int k = 0; k<NumCtrZ;k++);
		
	}
	//end of else if((theEntity->getEntityTypeTag())==3), I-section 3D
	
	//---------------------------------------------------Meshing 3D Composite beam-----------------------------------------------
	else if((theEntity->getEntityTypeTag())==7)
	{
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(5);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(6);
		
		
		int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1) +(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+(NumCtrX_Slab+1)*(NumCtrY_Slab+1) ;
		int NumElesPerLayer = NumCtrX*NumCtrY*2 +NumCtrX_Web*NumCtrY_Web+ NumCtrX_Slab*NumCtrY_Slab;
		
		int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4, NodeTag5, NodeTag6, NodeTag7, NodeTag8;
		for (int k = 0; k<NumCtrZ;k++){
			//along the length
		//Generating Elements for LowerFlange
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag+ NumElesPerLayer*k + NumCtrX*i+j;
				NodeTag1= OriginNodeTag+ NumNodesPerLayer*k + (NumCtrX+1)*i+j;
				NodeTag2= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*i+j+1;
				NodeTag3= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+ NumNodesPerLayer*k+ (NumCtrX+1)*(i+1)+j;
				NodeTag5= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*i+j;
				NodeTag6= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*i+j+1;
				NodeTag7= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*(i+1)+j+1;
				NodeTag8= OriginNodeTag+ NumNodesPerLayer*(k+1) + (NumCtrX+1)*(i+1)+j;
				
				
				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);
					
					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding elements in 3D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		
		//Generating Elements for Web
		int JunctionLTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)-(NumCtrX+1)+(NumCtrX-NumCtrX_Web)/2;
        int JunctionUTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX-NumCtrX_Web)/2;

		for(int i=0;i<NumCtrY_Web; i++) {
			for(int j=0; j<NumCtrX_Web; j++){
				EleTag = OriginEleTag+NumElesPerLayer*k + NumCtrX_Web*i+j+NumCtrX*NumCtrY;
			
				if(i==0)
				{
					NodeTag1= JunctionLTag+ NumNodesPerLayer *k +j;
					NodeTag2= JunctionLTag+ NumNodesPerLayer *k +j+1;
					NodeTag3= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
					
					NodeTag5= JunctionLTag+NumNodesPerLayer *(k+1) +j;
					NodeTag6= JunctionLTag+NumNodesPerLayer *(k+1) +j+1;
					NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag8= OriginNodeTag+NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
					
				}
				
				else if(i==NumCtrY_Web-1)
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
					NodeTag3= JunctionUTag + NumNodesPerLayer *k +j+1;
					NodeTag4= JunctionUTag + NumNodesPerLayer *k +j;
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
					NodeTag7= JunctionUTag+ NumNodesPerLayer *(k+1) +j+1;
					NodeTag8= JunctionUTag+ NumNodesPerLayer *(k+1) +j;
					
				}
				else
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
					NodeTag3= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag4= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
					NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
					NodeTag8= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					
					
				}

				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);
				
					result=theHTDomain->addElement(TempEle);
				    
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
					return(-1);
				}
				
			}
		}
		
		//Generating Elements for UpperFlange
		for(int i=0;i<NumCtrY; i++) {
			for(int j=0; j<NumCtrX; j++){
				EleTag = OriginEleTag + NumElesPerLayer*k +NumCtrX*i+j+NumCtrX*NumCtrY+NumCtrX_Web*NumCtrY_Web;
				int JunctionUTagSlab= OriginNodeTag+NumNodesPerLayer- (NumCtrX_Slab+1)*(NumCtrY_Slab+1)+(NumCtrX_Slab-NumCtrX)/2;
				
				if(i==NumCtrY-1)
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1; 
					NodeTag3= JunctionUTagSlab + NumNodesPerLayer *k +j+1;
					NodeTag4= JunctionUTagSlab + NumNodesPerLayer *k +j;
					
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1; 
					NodeTag7= JunctionUTagSlab+ NumNodesPerLayer *(k+1) +j+1;
					NodeTag8= JunctionUTagSlab+ NumNodesPerLayer *(k+1) +j;
					
					
				}
				else
				{
					NodeTag1= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
					NodeTag2= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
					NodeTag3= OriginNodeTag+ NumNodesPerLayer *k + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
					NodeTag4= OriginNodeTag+ NumNodesPerLayer *k +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
				
					NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
					NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
					NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
					NodeTag8= OriginNodeTag+ NumNodesPerLayer *(k+1) +(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
					
				}
				
				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial, PhaseTransformation);

 					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		//end of for(int i=0;i<NumCtrY; i++);
		
		
		//Generating Elements for CompositeSlab
		for(int i=0;i<NumCtrY_Slab; i++) {
			for(int j=0; j<NumCtrX_Slab; j++){
				EleTag = OriginEleTag + NumElesPerLayer*(k+1) -NumCtrX_Slab*NumCtrY_Slab+NumCtrX_Slab*i+j;
				int OriginNodeTagforSlab =  NumNodesPerLayer-(NumCtrX_Slab+1)*(NumCtrY_Slab+1);
				
				NodeTag1= OriginNodeTag+ NumNodesPerLayer *k + OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j;
				NodeTag2= OriginNodeTag+ NumNodesPerLayer *k + OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j+1;
				NodeTag3= OriginNodeTag+ NumNodesPerLayer *k + OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+ NumNodesPerLayer *k + OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j;
				
				NodeTag5= OriginNodeTag+ NumNodesPerLayer *(k+1) + OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j;
				NodeTag6= OriginNodeTag+ NumNodesPerLayer *(k+1) + OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j+1;
				NodeTag7= OriginNodeTag+ NumNodesPerLayer *(k+1) + OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j+1;
				NodeTag8= OriginNodeTag+ NumNodesPerLayer *(k+1) +OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j;
					
				
				if(isHTDomain)
				{
					BrickEight* TempEle=0;
					
					TempEle = new BrickEight  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,NodeTag5,NodeTag6,NodeTag7,NodeTag8,*theHTMaterial1, PhaseTransformation1);

 					result=theHTDomain->addElement(TempEle);
				}
				else
				{
					opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 3D I-section Beam"<<endln;
					return(-1);
				}
			}
		}
		//end of adding elements for composite slab
		
	}
	//end of for (int k = 0; k<NumCtrZ;k++);
		
	}
	//end of else if((theEntity->getEntityTypeTag())==7), Composite_Section 3D
	
	//------Generating elements for 2D PROTECTED ISECTION 
	//---------------------------------------------------Meshing 2D PROTECTED beam-----------------------------------------------
	else if((theEntity->getEntityTypeTag())==6)
	{
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(4);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(5);
		
		
		int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1) +(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+(NumCtrX_Slab+1)*(NumCtrY_Slab+1) ;
		int NumElesPerLayer = NumCtrX*NumCtrY*2 +NumCtrX_Web*NumCtrY_Web+ NumCtrX_Slab*NumCtrY_Slab;
		
		//Generating elements for lower flange
			int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;
			for(int i=0;i<NumCtrY; i++) {
				for(int j=0; j<NumCtrX; j++){
					EleTag = OriginEleTag+NumCtrX*i+j;
					NodeTag1= OriginNodeTag+(NumCtrX+1)*i+j;
					NodeTag2= OriginNodeTag+(NumCtrX+1)*i+j+1;
					NodeTag3= OriginNodeTag+(NumCtrX+1)*(i+1)+j+1;
					NodeTag4= OriginNodeTag+(NumCtrX+1)*(i+1)+j;
					if(isHTDomain)
					{
					    QuadFour* TempEle=0;
				    
						TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
					
						result=theHTDomain->addElement(TempEle);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 2D composite Beam"<<endln;
						return(-1);
						
					}
				}
			}
		
			//Generating Elements for Web
			int JunctionLTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)-(NumCtrX+1)+(NumCtrX-NumCtrX_Web)/2;
	        int JunctionUTag=OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX-NumCtrX_Web)/2;

			for(int i=0;i<NumCtrY_Web; i++) {
				for(int j=0; j<NumCtrX_Web; j++){
					EleTag = OriginEleTag+NumCtrX_Web*i+j+NumCtrX*NumCtrY;
			
					if(i==0)
					{
						NodeTag1= JunctionLTag+j;
						NodeTag2= JunctionLTag+j+1;
						NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
						NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j; 
					}
					else if(i==NumCtrY_Web-1)
					{
						NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
						NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1; 
						NodeTag3= JunctionUTag+j+1;
						NodeTag4= JunctionUTag+j;
					
					}
					else
					{
						NodeTag1= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j;
						NodeTag2= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(i-1)+j+1;
						NodeTag3= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j+1;
						NodeTag4= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*i+j;
					
					}

					if(isHTDomain)
					{
					    QuadFour* TempEle=0;
				    
						TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
					
						result=theHTDomain->addElement(TempEle);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding nodes in 2D composite Beam"<<endln;
						return(-1);
						
					}
				
				}
			}
		
			//Generating Elements for UpperFlange
		
			for(int i=0;i<NumCtrY; i++) {
				for(int j=0; j<NumCtrX; j++){
					EleTag = OriginEleTag+NumCtrX*i+j+NumCtrX*NumCtrY+NumCtrX_Web*NumCtrY_Web;
					int JunctionUTagSlab= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+ (NumCtrX+1)*NumCtrY+(NumCtrX_Slab-NumCtrX)/2;
				
					if(i==NumCtrY-1)
					{
						NodeTag1= OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
						NodeTag2= OriginNodeTag + (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1; 
						NodeTag3= JunctionUTagSlab +j+1;
						NodeTag4= JunctionUTagSlab  +j;
					
					
					}
					else
					{
						NodeTag1= OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j;
						NodeTag2= OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*i+j+1;
						NodeTag3= OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j+1;
						NodeTag4= OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(i+1)+j;
					
					}
					
			
				
					if(isHTDomain)
					{
						QuadFour* TempEle=0;
				    
						TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
		
						result=theHTDomain->addElement(TempEle);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding elements in 2D composite Beam"<<endln;
						return(-1);
					}
				}
			}
			//end of for(int i=0;i<NumCtrY; i++);
			
			//Generating Elements for CompositeSlab
			for(int i=0;i<NumCtrY_Slab; i++) {
				for(int j=0; j<NumCtrX_Slab; j++){
					EleTag = OriginEleTag+NumCtrX*NumCtrY*2+NumCtrX_Web*NumCtrY_Web+NumCtrX_Slab*i+j;
					int OriginNodeTagforSlab =  NumNodesPerLayer-(NumCtrX_Slab+1)*(NumCtrY_Slab+1);
				
					NodeTag1= OriginNodeTag+  OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j;
					NodeTag2= OriginNodeTag+  OriginNodeTagforSlab +(NumCtrX_Slab+1)*i+j+1;
					NodeTag3= OriginNodeTag+ OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j+1;
					NodeTag4= OriginNodeTag+  OriginNodeTagforSlab +(NumCtrX_Slab+1)*(i+1)+j;
				
					
				
					if(isHTDomain)
					{
						QuadFour* TempEle=0;
				    
						TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial1, PhaseTransformation1);
		
						result=theHTDomain->addElement(TempEle);
					}
					else
					{
						opserr<<"WARNING: The Simple_Mesh asks for HT_Domain when adding elements in 2D Composite Beam"<<endln;
						return(-1);
					}
					
				}
			}
			//end of adding elements for composite slab
			
		
		}
		//end of else if((theEntity->getEntityTypeTag())==6) 

	    else if((theEntity->getEntityTypeTag())==11)
		{
		
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtr_Coat = theEntity->GetNumCtrlID()(4);

		int NumX = NumCtrX;
		int NumY = NumCtrY+ NumCtr_Coat*2;
		int NumXweb = NumCtrX_Web + NumCtr_Coat*2;
		int NumYweb = NumCtrY_Web - NumCtr_Coat*2;

		//Generating Elements for LowerFlange
		int EleTag, NodeTag1, NodeTag2, NodeTag3, NodeTag4;
		
		for(int i=0;i<NumY; i++) {
			for(int j=0; j<NumX; j++){
				EleTag = OriginEleTag+NumX*i+j;
				NodeTag1= OriginNodeTag+(NumX+1)*i+j;
				NodeTag2= OriginNodeTag+(NumX+1)*i+j+1;
				NodeTag3= OriginNodeTag+(NumX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+(NumX+1)*(i+1)+j;
				
			    QuadFour* TempEle=0;

				if(i>=NumCtr_Coat &&i < NumY-NumCtr_Coat)
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
				else if((i>= NumY-NumCtr_Coat) && ((j >= (NumCtrX - NumCtrX_Web)/2) && (j < (NumCtrX + NumCtrX_Web)/2)))
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
				else 
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial1, PhaseTransformation);
				
				result=theHTDomain->addElement(TempEle);
				
			}
		}
		
		//Generating Elements for Web
		int JunctionLTag=OriginNodeTag+(NumX+1)*(NumY+1)-(NumX+1)+(NumX-NumXweb)/2;
        int JunctionUTag=OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1)+(NumX-NumXweb)/2;

		for(int i=0;i<NumYweb; i++) {
			for(int j=0; j<NumXweb; j++){
				EleTag = OriginEleTag+NumXweb*i+j+NumX*NumY;
			
				if(i==0)
				{
					NodeTag1= JunctionLTag+j;
					NodeTag2= JunctionLTag+j+1;
					NodeTag3= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*i+j+1;
					NodeTag4= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*i+j; 
				}
				else if(i==NumYweb-1)
				{
					NodeTag1= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(i-1)+j+1; 
					NodeTag3= JunctionUTag+j+1;
					NodeTag4= JunctionUTag+j;
					
				}
				else
				{
					NodeTag1= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(i-1)+j;
					NodeTag2= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(i-1)+j+1;
					NodeTag3= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*i+j+1;
					NodeTag4= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*i+j;
					
				}

				
			    QuadFour* TempEle=0;
			    if( (j >= NumCtr_Coat )&& (j < NumXweb - NumCtr_Coat))
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
				else 
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial1, PhaseTransformation);
				
				result=theHTDomain->addElement(TempEle);
				
				
			}
		}
		
		//Generating Elements for UpperFlange
		
		for(int i=0;i<NumY; i++) {
			for(int j=0; j<NumX; j++){
				EleTag = OriginEleTag+NumX*i+j+NumX*NumY+NumXweb*NumYweb;
				NodeTag1= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1)+(NumX+1)*i+j;
				NodeTag2= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1)+(NumX+1)*i+j+1;
				NodeTag3= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1)+(NumX+1)*(i+1)+j+1;
				NodeTag4= OriginNodeTag+(NumX+1)*(NumY+1)+(NumXweb+1)*(NumYweb-1)+(NumX+1)*(i+1)+j;
			
				
				QuadFour* TempEle=0;

				if((i< NumCtr_Coat )&& (j >= (NumCtrX - NumCtrX_Web)/2) && (j < (NumCtrX + NumCtrX_Web)/2))
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
				else if(i>=NumCtr_Coat &&i < NumY-NumCtr_Coat)
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial, PhaseTransformation);
				else
					TempEle = new QuadFour  (EleTag,NodeTag1,NodeTag2,NodeTag3,NodeTag4,*theHTMaterial1, PhaseTransformation);
	
				result=theHTDomain->addElement(TempEle);
		
			}
		}
		//end of for(int i=0;i<NumY; i++);
		
	}
		//end of else if((theEntity->getEntityTypeTag())==11) 
		
#ifdef _DEBUG
		if(isHTDomain)
		opserr<<"SimpleMesh has successfully generated "<<(theHTDomain->getNumElements())-OriginEleTag+1 <<" HeatTransfer Elements..."<<endln;
		else
		opserr<<"SimpleMesh has successfully generated "<<(theDomain->getNumElements())-OriginNodeTag+1 <<" elements..."<<endln;
#endif
}







//The following code is written for selecting Nodes and elements

int
Simple_Mesh::SelectingNodes(ID& NodesRange,int crdTag,double MinValue, double MaxValue,double Tolerance)
{
	//const ID& TempNodesRange = NodesRange; 
	int iniNumOfNodes = 0;
	bool iniSelecting = true;
	//check if it is selecting nodes in the whole domain
	if (NodesRange.Size() !=0) 
	{
		iniNumOfNodes = NodesRange.Size();
		iniSelecting = false;
	}
	else 
	{
		iniNumOfNodes = this->GetNumOfNodes();
		iniSelecting = true;
	}

	vector<int> SelectedNodes;
	int NodeTag =0;
	for(int i=0;i<iniNumOfNodes;i++){
		if(isHTDomain)
		{
			if(iniSelecting)
				NodeTag = OriginNodeTag+i;
			else
				NodeTag	= NodesRange(i);
			//opserr<<(theHTDomain->getNode(NodeTag)->getCrds())<<"     ";
			double NodalCrd = (theHTDomain->getNode(NodeTag)->getCrds())(crdTag);
			if((NodalCrd<=MaxValue+Tolerance)&&(NodalCrd>=MinValue-Tolerance)){
				SelectedNodes.push_back(NodeTag);
			}
		 }
		else 
		{
			if(iniSelecting)
				NodeTag = OriginNodeTag+i;
			else
				NodeTag	= NodesRange(i);

			Node* theNode = theDomain->getNode(NodeTag);
			if((theNode->getCrds()(crdTag)<=MaxValue+Tolerance)&&(theNode->getCrds()(crdTag)>=MinValue-Tolerance))
				SelectedNodes.push_back(NodeTag);
		}
	  
	}
	int NewIDsize = SelectedNodes.size();
	opserr<<"SimpleMesh::SelectingNodes has selected "<<NewIDsize<<" Nodes.."<<endln;

	NodesRange.resize(NewIDsize);
	for (int i = 0; i< NewIDsize; i++) {
		NodesRange(i)=SelectedNodes[i];
	}
#ifdef _DEBUG
	opserr<<NodesRange<<endln;
#endif
  return 0;
}

int Simple_Mesh::SelectingNodesbyFace(ID& NodesRange, int FaceTag) {
	//1D Line
	if((theEntity->getEntityTypeTag())==5){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		if (FaceTag==1) {
			NodesRange.resize(1);
			NodesRange(0) = OriginNodeTag;
		}else if (FaceTag==2) {
			NodesRange.resize(1);
			NodesRange(0) = OriginNodeTag+NumCtrX;
        }
		else {
		opserr<<"WARNING: SelectingNodesbyFace is only available for face 1 and 2 currently" <<endln;
		}
	}

   //2D block
	else if((theEntity->getEntityTypeTag())==0){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int NumCtrY= theEntity->GetNumCtrlID()(1);
		if (FaceTag==1) {
			NodesRange.resize(NumCtrX+1);
			for(int i =0;i<(NumCtrX+1);i++) {
			 NodesRange(i) = OriginNodeTag+i;
			}
		}else if (FaceTag==4) {
			NodesRange.resize(NumCtrX+1);
			for(int i =0;i<NumCtrX+1;i++) {
				NodesRange(i) = OriginNodeTag+(NumCtrX+1)*NumCtrY+i;
			}
		}else if (FaceTag==3) {
			NodesRange.resize(NumCtrY+1);
			for(int i =0;i<NumCtrX+1;i++) {
			 NodesRange(i) = OriginNodeTag+(NumCtrX+1)*i+NumCtrX;
			}
		}else if (FaceTag==2) {
			NodesRange.resize(NumCtrY+1);
			for(int i =0;i<NumCtrY+1;i++) {
			 NodesRange(i) = OriginNodeTag+(NumCtrX+1)*i;
        }
      }
	 else {
		opserr<<"WARNING: SelectingNodesbyFace is only available for face 1 and 4 currently" <<endln;
	  }
	}
	//2D I section Beam
	else if((theEntity->getEntityTypeTag())==1 || theEntity->getEntityTypeTag()==11 ){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);

		if(theEntity->getEntityTypeTag()==11){
			
			int NumCtr_Coat = theEntity->GetNumCtrlID()(4);

			int NumX = NumCtrX;
			int NumY = NumCtrY+ NumCtr_Coat*2;
			int NumXweb = NumCtrX_Web + NumCtr_Coat*2;
			int NumYweb = NumCtrY_Web - NumCtr_Coat*2;

			NumCtrX = NumX ; NumCtrY = NumY ;
			NumCtrX_Web = NumXweb ; NumCtrY_Web = NumYweb;

		}


		if (FaceTag==1) {
			NodesRange.resize(NumCtrX+1);
			for(int i =0;i<NumCtrX+1;i++) {
				NodesRange(i) = OriginNodeTag+i;
			}
		}
		else if (FaceTag==12) {
			NodesRange.resize(NumCtrX+1);
			//Temporarily starts from nodetag = The existing nodes in other mesh+ Upperflange other than top line+lower flange+web;
			int TempOriginNode = OriginNodeTag+(NumCtrX+1)*NumCtrY+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
			for(int i =0;i<(NumCtrX+1);i++) {
				NodesRange(i) = TempOriginNode+i;
			}
		}
		else {
			opserr<<"WARNING: SelectingNodesbyFace is only available for facetag 1 and 2 currently" <<endln;
		}
	}
	//3D Brick
	else if((theEntity->getEntityTypeTag())==2){
			int NumCtrX= theEntity->GetNumCtrlID()(0);
			int NumCtrY= theEntity->GetNumCtrlID()(1);
			int NumCtrZ= theEntity->GetNumCtrlID()(2);
			int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1);

			if (FaceTag==1) {
				NodesRange.resize((NumCtrX+1)*(NumCtrZ+1));
				for(int j =0;j<(NumCtrZ+1);j++){
					for(int i =0;i<(NumCtrX+1);i++) {
						NodesRange((NumCtrX+1)*j+i) = OriginNodeTag+NumNodesPerLayer*j+i;					
					}
				}
			}
			else if(FaceTag==2) {
				NodesRange.resize((NumCtrX+1)*(NumCtrZ+1));
				int TempNodeOrigin = OriginNodeTag+ (NumCtrX+1)*NumCtrY;
				for(int j =0;j<(NumCtrZ+1);j++){
					for(int i =0;i<(NumCtrX+1);i++) {
						NodesRange((NumCtrX+1)*j+i) = TempNodeOrigin+NumNodesPerLayer*j+i;					
					}
				}
			}
			else if (FaceTag==3) {
			    NodesRange.resize((NumCtrX+1)*(NumCtrY+1));
				for(int i =0;i<(NumCtrX+1)*(NumCtrY+1);i++) {
					NodesRange(i) = OriginNodeTag+i;
				}
			}
			else if (FaceTag==5) {
				NodesRange.resize((NumCtrX+1)*(NumCtrY+1));
				for(int i =0;i<(NumCtrX+1)*(NumCtrY+1);i++) {
					NodesRange(i) = OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)*NumCtrZ+i;
				}
			}	
			else {
				opserr<<"WARNING: SelectingNodesbyFace is only available for facetag 1 and 2 currently" <<endln;
			}
	}
	//3D I Section Beam
	else if((theEntity->getEntityTypeTag())==3){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ= theEntity->GetNumCtrlID()(4);
		int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1)*2 +(NumCtrX_Web+1)*(NumCtrY_Web-1);

		if (FaceTag==1) {
			NodesRange.resize((NumCtrX+1)*(NumCtrZ+1));
			for(int j= 0;j<NumCtrZ+1;j++) {
				for(int i =0;i<NumCtrX+1;i++) {
				NodesRange((NumCtrX+1)*j+ i) = OriginNodeTag+NumNodesPerLayer*j+i;
				}
			}
		}
		else if (FaceTag==12) {
			NodesRange.resize(NumCtrX+1);
			//Temporarily starts from nodetag = The existing nodes in other mesh+ Upperflange other than top line+lower flange+web;
			int TempOriginNode = OriginNodeTag+(NumCtrX+1)*NumCtrY+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1);
			
			NodesRange.resize((NumCtrX+1)*(NumCtrZ+1));
			for(int j= 0;j<NumCtrZ+1;j++) {
				for(int i =0;i<NumCtrX+1;i++) {
				NodesRange((NumCtrX+1)*j+ i) =	TempOriginNode + NumNodesPerLayer*j+i;
				}
			}
		}
		else {
			opserr<<"WARNING: SelectingNodesbyFace is only available for facetag 1 and 2 currently" <<endln;
		}
	}
	
	else
	{
		opserr<<"WARNING: SelectingNodesbyFace failed to identify the entity face tag" <<endln;
	}
	return 0;
}

int 
Simple_Mesh::SelectingEles(ID& ElesRange, const ID& NodesRange, int eleFaceTag)
{
	int iniNumOfEles = 0;
	bool iniSelecting = true;
	//check if it is selecting nodes in the whole domain
	if (ElesRange.Size() !=0) 
	{
		iniNumOfEles = ElesRange.Size();
		iniSelecting = false;
	}
	else 
	{
		iniNumOfEles = this->GetNumOfEles();
		iniSelecting = true;
	}

	vector<int> SelectedEles;
	int DetectedNodes =0;
	for(int i=0; i< iniNumOfEles; i++){
		if(isHTDomain)
		{
			int EleTag;
			if(iniSelecting)
				EleTag = OriginEleTag+i;
			else
				EleTag	= ElesRange(i);

			HeatTransferElement* theEle = theHTDomain->getElement(EleTag);
			
			ID NodesOnElementFace = theEle->getNodesOnFace(eleFaceTag);
			int NumPF = NodesOnElementFace.Size();
			for(int j=0; j<NodesRange.Size();j++){
				//For each selected node, the NodesonElementFace would at most have one in the nodes range
				for(int m =0; m< NumPF; m++){
					if (NodesOnElementFace(m)==NodesRange(j)){
						DetectedNodes++;
						break;
					}
				}
				if(DetectedNodes == NumPF)
					break;

			}
			//Seclecting eles by face or edge(the whole face or edge should lie on the selected plane
			if(DetectedNodes ==NumPF){
					SelectedEles.push_back(EleTag);
			}
			DetectedNodes=0;
			
		}else
		{
			opserr<<"WARNING: INVALID DOMAIN Recieved in SelectingEles"<<endln;
		}
			
	}
	int NewIDsize = SelectedEles.size();
	opserr<<"SimpleMesh::SelectingEles has selected "<<NewIDsize<<" elements.."<<endln;
    ElesRange.Zero();
	ElesRange.resize(NewIDsize);
	for (int i = 0; i< NewIDsize; i++) {
		ElesRange(i)=SelectedEles[i];
	}
  return 0;
}


int Simple_Mesh::SelectingElesbyFace(ID& ElesRange, int FaceTag, int& eleFaceID) {
	if((theEntity->getEntityTypeTag())==5){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		if (FaceTag==1) {
			eleFaceID=1;
			ElesRange.resize(1);
			ElesRange(0) = OriginEleTag;
		}else if (FaceTag==2) {
			eleFaceID=2;
			ElesRange.resize(1);
			ElesRange(0) = OriginEleTag+NumCtrX-1;
			
		}
		else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of block" <<endln;
		}
	}
	else if((theEntity->getEntityTypeTag())==1||(theEntity->getEntityTypeTag())==11){
    //2D Isection Beam;
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);

		if(theEntity->getEntityTypeTag()==11){
			
			int NumCtr_Coat = theEntity->GetNumCtrlID()(4);

			int NumX = NumCtrX;
			int NumY = NumCtrY+ NumCtr_Coat*2;
			int NumXweb = NumCtrX_Web + NumCtr_Coat*2;
			int NumYweb = NumCtrY_Web - NumCtr_Coat*2;

			NumCtrX = NumX ; NumCtrY = NumY ;
			NumCtrX_Web = NumXweb ; NumCtrY_Web = NumYweb;

		}
		if (FaceTag==1) {
			ElesRange.resize(NumCtrX);
      eleFaceID=1;
			for(int i =0;i<NumCtrX;i++) {
				ElesRange(i) = OriginEleTag+i;
			}
		}
    else if (FaceTag==4) {
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
       eleFaceID=3;
      for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
        ElesRange(i) = OriginEleTag+NumCtrX*(NumCtrY-1)+i;
      }
    }
    else if (FaceTag==5) {
       eleFaceID=3;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
      for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
        ElesRange(i) = OriginEleTag+NumCtrX*(NumCtrY-1)+NumCtrX/2+NumCtrX_Web/2+i;
      }
    }		else if (FaceTag==6) {
      eleFaceID=4;
			ElesRange.resize(NumCtrY_Web);
			for(int i =0;i<NumCtrY_Web;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*i;
			}
		}
		else if (FaceTag==7) {
      eleFaceID=2;
			ElesRange.resize(NumCtrY_Web);
			for(int i =0;i<NumCtrY_Web;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*i+NumCtrX_Web-1;
			}
		}
    else if (FaceTag==8) {
      eleFaceID=1;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
      for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
        ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+i;
      }
    }
    else if (FaceTag==9) {
      eleFaceID=1;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
      for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
        ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX/2+NumCtrX_Web/2+i;
      }
    }
		else if (FaceTag==12) {
      eleFaceID=3;
			ElesRange.resize(NumCtrX);
			for(int i =0;i<NumCtrX;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX*(NumCtrY-1)+ i;
			}
		}
		else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of I section beam" <<endln;
		}
	}
	//Block2D
	else if((theEntity->getEntityTypeTag())==0){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		if (FaceTag==1) {
      eleFaceID=1;
			ElesRange.resize(NumCtrX);
			for(int i =0;i<NumCtrX;i++) {
				ElesRange(i) = OriginEleTag+i;
				}
		}else if (FaceTag==3) {
      eleFaceID=2;
			ElesRange.resize(NumCtrY);
			for(int i =0;i<NumCtrY;i++) {
				ElesRange(i) = OriginEleTag+ NumCtrX*i+NumCtrX-1;
			}
		}
		else if (FaceTag==2) {
      eleFaceID=4;
			ElesRange.resize(NumCtrY);
			for(int i =0;i<NumCtrY;i++) {
				ElesRange(i) = OriginEleTag+ NumCtrX*i;
			}
		}
		else if (FaceTag==4) {
      eleFaceID=3;
			ElesRange.resize(NumCtrX);
			for(int i =0;i<NumCtrX;i++) {
				ElesRange(i) = OriginEleTag+ NumCtrX*(NumCtrY-1)+ i;
			}
		}
		else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of block" <<endln;
		}
	}
	//Brick3D
	else if((theEntity->getEntityTypeTag())==2){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrZ= theEntity->GetNumCtrlID()(2);
		if (FaceTag==1) {
      eleFaceID=3;
			ElesRange.resize(NumCtrX*NumCtrZ);
			for(int j= 0; j<NumCtrZ;j++) {
				for(int i =0;i<NumCtrX;i++){
					ElesRange(NumCtrX*j+i) = OriginEleTag+NumCtrX*NumCtrY*j+ i;
				}
			}
		}else if(FaceTag ==2) {
      eleFaceID=5;
			ElesRange.resize(NumCtrX*NumCtrZ);
			for(int j= 0; j<NumCtrZ;j++) {
				for(int i =0;i<NumCtrX;i++){
					ElesRange(NumCtrX*j+i) = OriginEleTag+NumCtrX*NumCtrY*j+ NumCtrX*(NumCtrY-1)+ i;
				}
			}
		}else if(FaceTag==3) {
      eleFaceID=1;
			ElesRange.resize(NumCtrX*NumCtrY);
			for(int i =0;i<NumCtrX*NumCtrY;i++) {
				ElesRange(i) = OriginEleTag+i;
			}
		}else if (FaceTag==5) {
      eleFaceID=2;
			ElesRange.resize(NumCtrX*NumCtrY);
			for(int i =0;i<NumCtrX*NumCtrY;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY*(NumCtrZ-1)+i;
			}
		}else {
			opserr<<"WARNING: SelectingElesbyFace is only available for facetag 1 and 2 currently" <<endln;
		}
	}
	//I section Beam3D
	else if((theEntity->getEntityTypeTag())==3){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ= theEntity->GetNumCtrlID()(4);
    int NumElesPerLayer = (NumCtrX*NumCtrY*2+NumCtrX_Web*NumCtrY_Web);
		
    if (FaceTag==1) {
      eleFaceID=3;
			ElesRange.resize(NumCtrX*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX;i++) {
			   ElesRange(NumCtrX*j+i) = OriginEleTag+NumElesPerLayer*j+i;
			 }
			}
    }
    else if (FaceTag==4) {
      eleFaceID=5;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
      for(int j=0; j<NumCtrZ;j++){
        for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
          ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*(NumCtrY-1)+i;
        }
      }
    }else if (FaceTag==5) {
      eleFaceID=5;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
      for(int j=0; j<NumCtrZ;j++){
        for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
          ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*(NumCtrY-1)+(NumCtrX+NumCtrX_Web)/2+i;
        }
      }
		}
    else if (FaceTag==6) {
      eleFaceID=6;
			ElesRange.resize(NumCtrY_Web*NumCtrZ);
			for(int j=0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrY_Web;i++) {
			  ElesRange(NumCtrY_Web*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*i;
			 }
			}
		}
    else if (FaceTag==7) {
      eleFaceID=4;
			ElesRange.resize(NumCtrY_Web*NumCtrZ);
			for(int j=0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrY_Web;i++) {
			  ElesRange(NumCtrY_Web*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*i+NumCtrX_Web-1;
			 }
			}
		}
    else if (FaceTag==8) {
      eleFaceID=3;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
      for(int j=0; j<NumCtrZ;j++) {
        for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
          ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+i;
        }
      }
		}
    else if (FaceTag==9) {
      eleFaceID=3;
      ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
      for(int j=0; j<NumCtrZ;j++) {
        for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
          ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ (NumCtrX+NumCtrX_Web)/2+i;
        }
      }
		}
    else if (FaceTag==12) {
      eleFaceID=5;
			ElesRange.resize(NumCtrX*NumCtrZ);
			for(int j =0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX;i++) {
			  ElesRange(NumCtrX*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX*(NumCtrY-1)+ i;
			 }
			}
		}else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of I section beam" <<endln;
		}
	}
	//Composite Beam3d
	else if((theEntity->getEntityTypeTag())==7){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ= theEntity->GetNumCtrlID()(4);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(5);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(6);
		
        int NumElesPerLayer = (NumCtrX*NumCtrY*2+NumCtrX_Web*NumCtrY_Web+NumCtrX_Slab*NumCtrY_Slab);
		if (FaceTag==1) {
      eleFaceID=3;
			ElesRange.resize(NumCtrX*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX;i++) {
			   ElesRange(NumCtrX*j+i) = OriginEleTag+NumElesPerLayer*j+i;
			 }
			}
		}else if (FaceTag==4) {
      eleFaceID=5;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
			  ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*(NumCtrY-1)+i;
			 }
			}
		}else if (FaceTag==5) {
      eleFaceID=5;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
			  ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*(NumCtrY-1)+(NumCtrX+NumCtrX_Web)/2+i;
			 }
			}
			
		}else if (FaceTag==6) {
      eleFaceID=6;
			ElesRange.resize(NumCtrY_Web*NumCtrZ);
			for(int j=0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrY_Web;i++) {
			  ElesRange(NumCtrY_Web*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*i;
			 }
			}
		}else if (FaceTag==7) {
      eleFaceID=4;
			ElesRange.resize(NumCtrY_Web*NumCtrZ);
			for(int j=0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrY_Web;i++) {
			  ElesRange(NumCtrY_Web*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*i+NumCtrX_Web-1;
			 }
			}
		}else if (FaceTag==8) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++) {
			 for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
			  ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+i;
			 }
			}
		}else if (FaceTag==9) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2*NumCtrZ);
			for(int j=0; j<NumCtrZ;j++) {
			 for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
			  ElesRange((NumCtrX-NumCtrX_Web)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ (NumCtrX+NumCtrX_Web)/2+i;
			 }
			}
			
		}else if (FaceTag==12) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX_Slab-NumCtrX)/2*NumCtrZ);
			for(int j =0;j<NumCtrZ;j++){
			 for(int i =0;i<(NumCtrX_Slab-NumCtrX)/2;i++) {
			  ElesRange((NumCtrX_Slab-NumCtrX)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX*NumCtrY+ i;
			 }
			}
		}else if (FaceTag==13) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX_Slab-NumCtrX)/2*NumCtrZ);
			for(int j =0;j<NumCtrZ;j++){
			 for(int i =0;i<(NumCtrX_Slab-NumCtrX)/2;i++) {
			  ElesRange((NumCtrX_Slab-NumCtrX)/2*j+i) = OriginEleTag+NumElesPerLayer*j+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX*NumCtrY+ (NumCtrX_Slab+NumCtrX)/2+i;
			 }
			}
		}else if (FaceTag==16) {
      eleFaceID=5;
			ElesRange.resize(NumCtrX_Slab*NumCtrZ);
			for(int j =0;j<NumCtrZ;j++){
			 for(int i =0;i<NumCtrX_Slab;i++) {
			  ElesRange(NumCtrX_Slab*j+i) = OriginEleTag+NumElesPerLayer*(j+1)-NumCtrX_Slab+i;
			 }
			}	
		}else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of composite beam" <<endln;
		}
	}
	//Composite Beam2D
	else if((theEntity->getEntityTypeTag())==6){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(4);
		int NumCtrY_Slab = theEntity->GetNumCtrlID()(5);
		
		if (FaceTag==1) {
      eleFaceID=1;
			ElesRange.resize(NumCtrX);
			for(int i =0;i<NumCtrX;i++) {
				ElesRange(i) = OriginEleTag+i;
			}
		}
		else if (FaceTag==4) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
			for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*(NumCtrY-1)+i;
			}
		}
		else if (FaceTag==5) {
      eleFaceID=3;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
			for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*(NumCtrY-1)+NumCtrX/2+NumCtrX_Web/2+i;
			}
		}
		
		else if (FaceTag==6) {
      eleFaceID=4;
			ElesRange.resize(NumCtrY_Web);
			for(int i =0;i<NumCtrY_Web;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*i;
			}
		}
		else if (FaceTag==7) {
      eleFaceID=2;
			ElesRange.resize(NumCtrY_Web);
			for(int i =0;i<NumCtrY_Web;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*i+NumCtrX_Web-1;
			}
		}
		else if (FaceTag==8) {
      eleFaceID=1;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
			for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+i;

			}
			
		}
		else if (FaceTag==9) {
      eleFaceID=1;
			ElesRange.resize((NumCtrX-NumCtrX_Web)/2);
			for(int i =0;i<NumCtrX/2-NumCtrX_Web/2;i++) {
				ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX/2+NumCtrX_Web/2+i;
			}
			
		}
	
		else if (FaceTag==12) {
      eleFaceID=1;
      ElesRange.resize((NumCtrX_Slab-NumCtrX)/2);
		
		 for(int i =0;i<(NumCtrX_Slab-NumCtrX)/2;i++) {
		  ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY*2+ NumCtrX_Web*NumCtrY_Web+ i;
		 }
		 }
		else if (FaceTag==13) {
      eleFaceID=1;
      ElesRange.resize((NumCtrX_Slab-NumCtrX)/2);
		 for(int i =0;i<(NumCtrX_Slab-NumCtrX)/2;i++) {
		  ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY+ NumCtrX_Web*NumCtrY_Web+ NumCtrX*NumCtrY+ (NumCtrX_Slab+NumCtrX)/2+i;
		 }
		}
		else if (FaceTag==16) {
      eleFaceID=3;
      ElesRange.resize(NumCtrX_Slab);
		 for(int i =0;i<NumCtrX_Slab;i++) {
		  ElesRange(i) = OriginEleTag+NumCtrX*NumCtrY*2+ NumCtrX_Web*NumCtrY_Web+ NumCtrX_Slab*(NumCtrY_Slab-1)+i;
		 }

		}
		else {
				opserr<<"WARNING: invalid face tag for selecting eles on face of I section beam" <<endln;
		}
	}
	//end of selecting eles by face for 2d Composite beam


	
	else
	{
		opserr<<"WARNING: SelectingNodesbyFace failed to identify the entity face tag" <<endln;
	}
	return 0;
}



int Simple_Mesh::GetNodesForRecorder(ID& NodesRange, int dimTag, double PlaneLoc) {
	if((theEntity->getEntityTypeTag())==1||(theEntity->getEntityTypeTag())==6){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		if(dimTag==2){
			int MonitorNodes[9];
			MonitorNodes[0]= OriginNodeTag+NumCtrX/2;
			MonitorNodes[1]= OriginNodeTag+(NumCtrX+1)*(NumCtrY)+NumCtrX/2;
			MonitorNodes[2]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(1*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[3]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(2*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[4]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(3*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[5]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(4*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[6]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(5*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[7]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX/2;
			MonitorNodes[8]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+NumCtrX/2;
			NodesRange.resize(9);
			for(int i =0;i<9;i++){
				NodesRange(i)=MonitorNodes[i];
			}
		}
		else if (dimTag==3){
			int MonitorNodes[15];
			MonitorNodes[0]= OriginNodeTag+(NumCtrX+1)*NumCtrY+NumCtrX/2;
			MonitorNodes[1]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[2]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/2-1)+NumCtrX_Web/2;
			MonitorNodes[3]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(3*NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[4]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX/2;
			MonitorNodes[5]= OriginNodeTag+(NumCtrX+1)*(NumCtrY/2)+NumCtrX/10;
			MonitorNodes[6]= OriginNodeTag+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[7]= OriginNodeTag+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[8]= OriginNodeTag+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[9]= OriginNodeTag+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			MonitorNodes[10]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+1*NumCtrX/10;
			MonitorNodes[11]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[12]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[13]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[14]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			NodesRange.resize(15);
			for(int i =0;i<15;i++){
				NodesRange(i)=MonitorNodes[i];
			}

		}
	}
	//end of ISection
	else if((theEntity->getEntityTypeTag())==11){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtr_Coat = theEntity->GetNumCtrlID()(4);
		if(dimTag==2){
			opserr<<"Not Available"<<endln;
		}
		else if (dimTag==3){
			int MonitorNodes[15];  
			int OriginTag = OriginNodeTag+(NumCtrX+1)*(NumCtrY+NumCtr_Coat*2+1);
			MonitorNodes[0]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+NumCtr_Coat)+NumCtrX/2;
			MonitorNodes[1]= OriginTag +(NumCtrX_Web+ 2*NumCtr_Coat +1)*(NumCtrY_Web/4-NumCtr_Coat/2-1)+  NumCtrX_Web/2+NumCtr_Coat;
			MonitorNodes[2]= OriginTag +(NumCtrX_Web+ 2*NumCtr_Coat +1)*(NumCtrY_Web/2-NumCtr_Coat-1)+  NumCtrX_Web/2+NumCtr_Coat;
			MonitorNodes[3]= OriginTag +(NumCtrX_Web+ 2*NumCtr_Coat +1)*(3*NumCtrY_Web/4-3*NumCtr_Coat/2-1)+ NumCtrX_Web/2+NumCtr_Coat;
			MonitorNodes[4]= OriginNodeTag+(NumCtrX+1)*(NumCtrY+NumCtr_Coat*2+1)+(NumCtrX_Web+ 2*NumCtr_Coat +1)*(NumCtrY_Web-NumCtr_Coat*2-1)
								+(NumCtrX+1)* NumCtr_Coat + NumCtrX/2 ;
			
			OriginTag = OriginNodeTag+(NumCtrX+1)*(NumCtrY/2 + NumCtr_Coat);
			MonitorNodes[5]= OriginTag+1;
			MonitorNodes[6]= OriginTag+((NumCtrX-NumCtrX_Web)/2- NumCtr_Coat)/2;
			MonitorNodes[7]= OriginTag + NumCtrX/2+1;
			MonitorNodes[8]= OriginTag +NumCtrX-((NumCtrX-NumCtrX_Web)/2- NumCtr_Coat)/2+1;
			MonitorNodes[9]= OriginTag+NumCtrX+1;
			
			OriginTag = OriginNodeTag+(NumCtrX+1)*(NumCtrY+2*NumCtr_Coat+1)+(NumCtrX_Web+2*NumCtr_Coat+1)*(NumCtrY_Web-2*NumCtr_Coat-1)+ 
				        (NumCtrX+1)*(NumCtrY/2 + NumCtr_Coat);
			MonitorNodes[10]= OriginTag;
			MonitorNodes[11]= OriginTag+((NumCtrX-NumCtrX_Web)/2- NumCtr_Coat)/2;
			MonitorNodes[12]= OriginTag + NumCtrX/2;
			MonitorNodes[13]= OriginTag +NumCtrX-((NumCtrX-NumCtrX_Web)/2- NumCtr_Coat)/2;
			MonitorNodes[14]= OriginTag+NumCtrX;
			NodesRange.resize(15);
			for(int i =0;i<15;i++){
				NodesRange(i)=MonitorNodes[i];
			}

		}
	}
	else if((theEntity->getEntityTypeTag())==3){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		//int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1)*2 +(NumCtrX_Web+1)*(NumCtrY_Web-1);
		Vector Seeds= theEntity->GetSeeds(3);
		int PlaneTag;
		for(int i=0; i<=NumCtrZ; i++){
			if(fabs(PlaneLoc-Seeds(i))<1e-5){
				PlaneTag=i;
				break;
			}
		}
		int PlaneOriginTag = OriginNodeTag+	((NumCtrX+1)*(NumCtrY+1)*2 +(NumCtrX_Web+1)*(NumCtrY_Web-1))*PlaneTag;
		if(dimTag==2){
			int MonitorNodes[9];
			MonitorNodes[0]= PlaneOriginTag+NumCtrX/2;
			MonitorNodes[1]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY)+NumCtrX/2;
			MonitorNodes[2]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(1*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[3]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(2*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[4]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(3*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[5]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(4*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[6]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(5*NumCtrY_Web/6-1)+NumCtrX_Web/2;
			MonitorNodes[7]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX/2;
			MonitorNodes[8]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+NumCtrX/2;
			NodesRange.resize(9);
			for(int i =0;i<9;i++){
				NodesRange(i)=MonitorNodes[i];
			}
		}
		else if (dimTag==3){
			int MonitorNodes[15];
			
			//First 5 points are monitoring temperature in beam web;
			MonitorNodes[0]= PlaneOriginTag+(NumCtrX+1)*NumCtrY+NumCtrX/2;
			MonitorNodes[1]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[2]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/2-1)+NumCtrX_Web/2;
			MonitorNodes[3]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(3*NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[4]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX/2;
			//Second 5 points are monitoring temperature in beam lower flange;
			MonitorNodes[5]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+NumCtrX/10;
			MonitorNodes[6]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[7]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[8]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[9]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			//Last 5 points are monitoring temperature in beam upper flange;
			MonitorNodes[10]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+1*NumCtrX/10;
			MonitorNodes[11]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[12]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[13]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[14]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			NodesRange.resize(15);
			for(int i =0;i<15;i++){
				NodesRange(i)=MonitorNodes[i];
			}

		}
	}
	else if((theEntity->getEntityTypeTag())==7){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		int NumCtrX_Slab = theEntity->GetNumCtrlID()(5);
		int NumCtrY_Slab = theEntity->GetNumCtrlID()(6);
		//int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1)*2 +(NumCtrX_Web+1)*(NumCtrY_Web-1);
		Vector Seeds= theEntity->GetSeeds(3);
		int PlaneTag;
		for(int i=0; i<=NumCtrZ; i++){
			if(fabs(PlaneLoc-Seeds(i))<1e-5){
				PlaneTag=i;
				break;
			}
		}
		int PlaneOriginTag = OriginNodeTag+	((NumCtrX+1)*(NumCtrY+1) +(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+ (NumCtrX_Slab+1)*(NumCtrY_Slab+1))*PlaneTag;
		if(dimTag==2){
			int MonitorNodes[9];
			MonitorNodes[0]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+NumCtrX/2;
			MonitorNodes[1]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/10-1)+NumCtrX_Web/2;
			MonitorNodes[2]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(2*NumCtrY_Web/10-1)+NumCtrX_Web/2;
			MonitorNodes[3]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(4*NumCtrY_Web/10-1)+NumCtrX_Web/2;
			MonitorNodes[4]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(5*NumCtrY_Web/10)+NumCtrX/2;
			MonitorNodes[5]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(7*NumCtrY_Web/10-1)+NumCtrX/2;
			MonitorNodes[6]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(9*NumCtrY_Web/10-1)+NumCtrX_Web/2;
			MonitorNodes[7]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX_Web/2;
			MonitorNodes[8]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+NumCtrX/2;
			NodesRange.resize(9);
			for(int i =0;i<9;i++){
				NodesRange(i)=MonitorNodes[i];
			}
		}
		else if (dimTag==3){
			int MonitorNodes[15];
			
			//First 5 points are monitoring temperature in beam web;
			MonitorNodes[0]= PlaneOriginTag+(NumCtrX+1)*NumCtrY+NumCtrX/2;
			MonitorNodes[1]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[2]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web/2-1)+NumCtrX_Web/2;
			MonitorNodes[3]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(3*NumCtrY_Web/4-1)+NumCtrX_Web/2;
			MonitorNodes[4]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+NumCtrX/2;
			//Second 5 points are monitoring temperature in beam lower flange;
			MonitorNodes[5]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+NumCtrX/10;
			MonitorNodes[6]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[7]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[8]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[9]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			//Last 5 points are monitoring temperature in beam upper flange;
			MonitorNodes[10]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+1*NumCtrX/10;
			MonitorNodes[11]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+3*NumCtrX/10;
			MonitorNodes[12]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+5*NumCtrX/10;
			MonitorNodes[13]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+7*NumCtrX/10;
			MonitorNodes[14]= PlaneOriginTag+(NumCtrX+1)*(NumCtrY+1)+(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*(NumCtrY/2)+9*NumCtrX/10;
			NodesRange.resize(15);
			for(int i =0;i<15;i++){
				NodesRange(i)=MonitorNodes[i];
			}

		}
	}
	//end of get node for record compsosite beam
	else
		opserr<<"invalid input for recording nodes"<<endln;
	return 0;
}


int Simple_Mesh::GetNodesForRecorder(ID& NodesRange, const Vector& MonitorLocX, double PlaneLoc) {
	if((theEntity->getEntityTypeTag())==0){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int PlaneOriginTag=0;
		Vector yLocTags = Vector(9);
	
		PlaneOriginTag = OriginNodeTag;
		
		Vector Seeds1= theEntity->GetSeeds(1);
		Vector Seeds2= theEntity->GetSeeds(2);

		int MonitorLocX_Size=MonitorLocX.Size();
		ID TagLocx(MonitorLocX_Size);

		for(int k=0; k<MonitorLocX_Size;k++){
		 for(int i=0; i<=NumCtrX; i++){
			if(fabs(MonitorLocX(k)-Seeds1(i))<1e-5){
				TagLocx(k)=i;
			}
		 }
		}
		double Entity_thickness= fabs(Seeds2(NumCtrY)-Seeds2(0));
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrY; i++){
				if(fabs((Entity_thickness/8*j)+Seeds2(0)-Seeds2(i))<(Entity_thickness/NumCtrY/2)){
					yLocTags(j)=i;					
				}
			}
		}

#ifdef _DEBUG
	opserr<<"yLocTags: "<<yLocTags;
	opserr<<"SimpleMesh:TagLocx: "<<TagLocx;
#endif
		ID MonitorNodes(MonitorLocX_Size*9);
		for(int i=0;i<MonitorLocX_Size;i++){
			for(int j=0; j<9 ;j++){
			MonitorNodes(9*i+j)= PlaneOriginTag+(NumCtrX+1)*yLocTags(j)+TagLocx(i);
			}
		}
		
		NodesRange=MonitorNodes;
	
	}
	else if((theEntity->getEntityTypeTag())==2){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrZ= theEntity->GetNumCtrlID()(2);
		int PlaneOriginTag=0;
		int PlaneTag=0;
		Vector yLocTags = Vector(9);

		Vector Seeds1= theEntity->GetSeeds(1);
		Vector Seeds2= theEntity->GetSeeds(2);
		Vector Seeds3= theEntity->GetSeeds(3);
		
		double Entity_y_thickness= fabs(Seeds2(NumCtrY)-Seeds2(0));
		
		for(int i=0; i<=NumCtrZ; i++){
			if(fabs(PlaneLoc-Seeds3(i))<1e-5){
				PlaneTag=i;
				break;
			}
		}
		PlaneOriginTag = OriginNodeTag+	(NumCtrX+1)*(NumCtrY+1)*PlaneTag;
		
		int MonitorLocX_Size=MonitorLocX.Size();
		ID TagLocx(MonitorLocX_Size);

		for(int k=0; k<MonitorLocX_Size;k++){
		 for(int i=0; i<=NumCtrX; i++){
			if(fabs(MonitorLocX(k)-Seeds1(i))<1e-5){
				TagLocx(k)=i;
			}
		 }
		}
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrY; i++){
				if(fabs((Entity_y_thickness/8*j)+Seeds2(0)-Seeds2(i))<(Entity_y_thickness/NumCtrY/2)){
					yLocTags(j)=i;					
				}
			}
		}
		#ifdef _DEBUG
		opserr<<"yLocTags: "<<yLocTags;
		opserr<<"SimpleMesh:TagLocx: "<<TagLocx;
		#endif

		ID MonitorNodes(MonitorLocX_Size*9);
		for(int i=0;i<MonitorLocX_Size;i++){
			for(int j=0; j<9;j++){
			MonitorNodes(9*i+j)= PlaneOriginTag+(NumCtrX+1)*yLocTags(j)+TagLocx(i);
			}
		}
		
		NodesRange=MonitorNodes;
	
	}
	else if((theEntity->getEntityTypeTag())==7){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrZ = theEntity->GetNumCtrlID()(4);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(5);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(6);
		
		int NumNodesPerLayer = (NumCtrX+1)*(NumCtrY+1) +(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY+(NumCtrX_Slab+1)*(NumCtrY_Slab+1) ;
		
		int PlaneOriginTag=0;
		int PlaneTag=0;
		Vector yLocTags = Vector(9);

		Vector Seeds4= theEntity->GetSeeds(4); //composite slab along x
		Vector Seeds5= theEntity->GetSeeds(5);  //compsoite slab along y
		Vector Seeds3= theEntity->GetSeeds(3);   //composite slab along z
		
		double Entity_y_thickness= fabs(Seeds5(NumCtrY_Slab)-Seeds5(0));
		
		for(int i=0; i<=NumCtrZ; i++){
			if(fabs(PlaneLoc-Seeds3(i))<1e-5){
				PlaneTag=i;
				break;
			}
		}
		
		PlaneOriginTag = OriginNodeTag+	NumNodesPerLayer *(PlaneTag+1)-(NumCtrX_Slab+1)*(NumCtrY_Slab+1);
		
		int MonitorLocX_Size=MonitorLocX.Size();
		ID TagLocx(MonitorLocX_Size);
		//Looking for the corrsponding x loc tag;
		for(int k=0; k<MonitorLocX_Size;k++){
		 for(int i=0; i<=NumCtrX_Slab; i++){
			if(fabs(MonitorLocX(k)-Seeds4(i))<1e-5){
				TagLocx(k)=i;
			}
		 }
		}
		//Looking for the corrsponding y loc tag;
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrY_Slab; i++){
				if(fabs((Entity_y_thickness/8*j)+Seeds5(0)-Seeds5(i))<(Entity_y_thickness/NumCtrY_Slab/2)){
					yLocTags(j)=i;					
				}
			}
		}
		#ifdef _DEBUG
		opserr<<"yLocTags: "<<yLocTags;
		opserr<<"SimpleMesh:TagLocx: "<<TagLocx;
		#endif
			
		ID MonitorNodes(MonitorLocX_Size*9);
		for(int i=0;i<MonitorLocX_Size;i++){
			for(int j=0; j<9;j++){
			MonitorNodes(9*i+j)= PlaneOriginTag+(NumCtrX_Slab+1)*yLocTags(j)+TagLocx(i);
			}
		}
		
		NodesRange=MonitorNodes;
	
	}
	 else if((theEntity->getEntityTypeTag())==6){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrX_Web = theEntity->GetNumCtrlID()(2);
		int NumCtrY_Web = theEntity->GetNumCtrlID()(3);
		int NumCtrX_Slab= theEntity->GetNumCtrlID()(4);
        int NumCtrY_Slab= theEntity->GetNumCtrlID()(5);
		
		int PlaneOriginTag=0;
		Vector yLocTags = Vector(9);
	
		PlaneOriginTag = OriginNodeTag+ (NumCtrX+1)*(NumCtrY+1) +(NumCtrX_Web+1)*(NumCtrY_Web-1)+(NumCtrX+1)*NumCtrY;
		
		Vector Seeds3= theEntity->GetSeeds(3);
		Vector Seeds4= theEntity->GetSeeds(4);

		int MonitorLocX_Size=MonitorLocX.Size();
		ID TagLocx(MonitorLocX_Size);

		for(int k=0; k<MonitorLocX_Size;k++){
		 for(int i=0; i<=NumCtrX_Slab; i++){
			if(fabs(MonitorLocX(k)-Seeds3(i))<1e-5){
				TagLocx(k)=i;
			}
		 }
		}
		double Entity_thickness= fabs(Seeds4(NumCtrY_Slab)-Seeds4(0));
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrY_Slab; i++){
				if(fabs((Entity_thickness/8*j)+Seeds4(0)-Seeds4(i))<(Entity_thickness/NumCtrY_Slab/2)){
					yLocTags(j)=i;					
				}
			}
		}

#ifdef _DEBUG
	opserr<<"yLocTags: "<<yLocTags;
	opserr<<"SimpleMesh:TagLocx: "<<TagLocx;
#endif
		ID MonitorNodes(MonitorLocX_Size*9);
		for(int i=0;i<MonitorLocX_Size;i++){
			for(int j=0; j<9 ;j++){
			MonitorNodes(9*i+j)= PlaneOriginTag+(NumCtrX_Slab+1)*yLocTags(j)+TagLocx(i);
			}
		}
		
		NodesRange=MonitorNodes;
	
	}
	//end of if entitytype==6;
	
	else
		opserr<<"invalid input for recording nodes"<<endln;
	return 0;
}

int Simple_Mesh::GetNodesForRecorder(ID& NodesRange, double xLoc, double yLoc) {
	if((theEntity->getEntityTypeTag())==5){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
		Vector Seeds1= theEntity->GetSeeds(1);
		double Entity_x_thickness= fabs(Seeds1(NumCtrX)-Seeds1(0));
		//LocOriginTag = OriginNodeTag;
		Vector xLocTags = Vector(9);
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrX; i++){
				if(fabs((Entity_x_thickness/8*j)+Seeds1(0)-Seeds1(i))<(Entity_x_thickness/NumCtrX/2)){
					xLocTags(j)=i;					
				}
			}
		}
#ifdef _DEBUG
		opserr<<"xLocTags: "<<xLocTags;
#endif
		ID MonitorNodes(9);
		for(int i=0;i<9;i++){
			MonitorNodes(i)= OriginNodeTag+xLocTags(i);
		}
		
		NodesRange=MonitorNodes;
	}
	 else if((theEntity->getEntityTypeTag())==2){
		int NumCtrX= theEntity->GetNumCtrlID()(0);
        int NumCtrY= theEntity->GetNumCtrlID()(1);
		int NumCtrZ= theEntity->GetNumCtrlID()(2);
		int xLocTag=-1;
		int yLocTag=-1;
		Vector zLocTags = Vector(9);
		int LocOriginTag =0;
		Vector Seeds1= theEntity->GetSeeds(1);
		Vector Seeds2= theEntity->GetSeeds(2);
		Vector Seeds3= theEntity->GetSeeds(3);
		double Entity_z_thickness= fabs(Seeds3(NumCtrZ)-Seeds3(0));
		for(int i=0; i<=NumCtrX; i++){
			if(fabs(xLoc-Seeds1(i))<1e-5){
				xLocTag=i;
				break;
			}
		}
		for(int i=0; i<=NumCtrY; i++){
			if(fabs(yLoc-Seeds2(i))<1e-5){
				yLocTag=i;
				break;
			}
		}
		for(int j=0; j<9; j++){
			for(int i=0; i<=NumCtrZ; i++){
				if(fabs((Entity_z_thickness/8*j)+Seeds3(0)-Seeds3(i))<(Entity_z_thickness/NumCtrZ/2)){
					zLocTags(j)=i;					
				}
			}
		}
#ifdef _DEBUG
		opserr<<"zLocTags: "<<zLocTags;
#endif
		if(xLocTag<0||yLocTag<0){
			opserr<<"Simple_Mesh::GetNodesForRecorder fail to find xLocTag or yLocTag"<<endln;
		}
		LocOriginTag = OriginNodeTag+ (NumCtrX+1)*yLocTag+ xLocTag;

		ID MonitorNodes(9);
		for(int i=0;i<9;i++){
			MonitorNodes(i)= LocOriginTag+(NumCtrX+1)*(NumCtrY+1)*zLocTags(i);
		}
		
		NodesRange=MonitorNodes;
	
	}
	else
		opserr<<"invalid input for recording nodes"<<endln;
	return 0;
}
/*
int Simple_Mesh::getMeshTag()
{
	return MeshTag;
}

*/

//The following function is for getting number of elements (Mesh) along X-axis
const ID&
Simple_Mesh::getNumCtrlID()
	{return theEntity->GetNumCtrlID();}



