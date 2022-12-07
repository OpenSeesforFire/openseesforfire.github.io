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
//Entity Type:
//Line:5; Block:0;Brick:2;Isection:1;Isection3D:3; IsecProtected:11;composite 2D:6;composite3D:7


Simple_Mesh::Simple_Mesh(int tag, Simple_Entity* Entity,HeatTransferDomain* theDomain,HeatTransferMaterial* theHTMaterial, Vector& MeshCtrls,
						 HeatTransferMaterial* theHTMaterial1, bool numCtrl)
:TaggedObject(tag),theHTDomain(theDomain),isHTDomain(true), theHTMaterial(theHTMaterial),OriginLocs(0),EleParameters(0),
 theHTMaterial1(theHTMaterial1)
{
	theEntity = Entity;
	theEntity->setMeshTag(this->getTag());
	 theEntity->InitialMeshCtrl(MeshCtrls,numCtrl);

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

int Simple_Mesh::GeneratingNodes(const Vector& originLocs )
{
  if(OriginLocs==0){
    if (originLocs!=0) {
      OriginLocs.Zero();
      OriginLocs = originLocs;
   }
  }
  else
  {
#ifdef _DEBUG
    opserr<<"SimpleMesh::GeneratingNodes encountered redefined OriginLocs"<<endln;
#endif
  }
  OriginNodeTag = theHTDomain->getNumNodes()+1;
  int nDoF=1;

	if (theEntity->GenerateNodes(theHTDomain, nDoF, OriginLocs) < 0) {
		opserr << "Simple Mesh: " << this->getTag() << "failed to generate nodes" << endln;
		return -1;
	}
	else {

#ifdef _DEBUG
	 opserr << "SimpleMesh has successfully generated " << (theHTDomain->getNumNodes()) - OriginNodeTag+1 << " HeatTransfer Nodes..." << endln;
#endif

	}

	return 0;
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
    opserr<< "SimpleMesh " << this->getTag()<< "encounters redefined EleParameteres"<<endln;
  }
  

#ifdef _DEBUG
	opserr<<"SimpleMesh "<<this->getTag()<<" starts to generate elements...."<< theEntity->getEntityTypeTag() <<endln;
#endif

  
	if(isHTDomain)
		OriginEleTag = theHTDomain->getNumElements()+1;

		//Generating Elements for Simple_Line
	int result = theEntity->GenerateEles(theHTDomain, EleParameters, theHTMaterial,theHTMaterial1);

	if (result < 0) {
		opserr << "HeatTransferDomain Failed to generate elements";
		return -1;
	}

	
		
#ifdef _DEBUG
		if(isHTDomain)
		opserr<<"SimpleMesh has successfully generated "<<(theHTDomain->getNumElements())-OriginEleTag+1 <<" HeatTransfer Elements..."<<endln;
#endif

return 0;
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
			/*
			if(iniSelecting)
				NodeTag = OriginNodeTag+i;
			else
				NodeTag	= NodesRange(i);

			Node* theNode = theDomain->getNode(NodeTag);
			if((theNode->getCrds()(crdTag)<=MaxValue+Tolerance)&&(theNode->getCrds()(crdTag)>=MinValue-Tolerance))
				SelectedNodes.push_back(NodeTag);*/
		}
	  
	}
	int NewIDsize = SelectedNodes.size();

	//opserr<<"SimpleMesh::SelectingNodes has selected "<<NewIDsize<<" Nodes for crd"<<crdTag<<endln;

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
			for(int i =0;i<NumCtrY+1;i++) {
			 NodesRange(i) = OriginNodeTag+(NumCtrX+1)*i+NumCtrX;
			}
		}else if (FaceTag==2) {
			NodesRange.resize(NumCtrY+1);
			for(int i =0;i<NumCtrY+1;i++) {
			 NodesRange(i) = OriginNodeTag+(NumCtrX+1)*i;
        }
      }
	 else {
		opserr<<"WARNING: SelectingNodesbyFace is unable to identify the face tag" <<endln;
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
	//opserr<<"SimpleMesh::SelectingEles has selected "<<NewIDsize<<" elements.."<<endln;
    ElesRange.Zero();
	ElesRange.resize(NewIDsize);
	for (int i = 0; i< NewIDsize; i++) {
		ElesRange(i)=SelectedEles[i];
	}
  return 0;
}


int Simple_Mesh::SelectingElesbyFace(ID& ElesRange, int FaceTag, int& eleFaceID) {
	if((theEntity->getEntityTypeTag())==5){
		//1D Line; Face Tag 1&2
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
    //2D Isection Beam or IsectionBeamProtected;
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



