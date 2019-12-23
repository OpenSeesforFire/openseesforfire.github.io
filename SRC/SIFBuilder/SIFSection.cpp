/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                   
/**********************************************************************
** This project is aiming to provide a Tcl interface to define models  **
** for simulating structural behaviours under fire action.           **
** Developed by:                  `                                   **
**   Liming Jiang (liming.jiang@ed.ac.uk)                            **
**   Praven Kamath(Praveen.Kamath@ed.ac.uk)                          **
**   Xu Dai(X.Dai@ed.ac.uk)                                          **
**   Asif Usmani(asif.usmani@ed.ac.uk)                               **
**********************************************************************/
// $Revision: 2.4.0.1 $
// This file constructs the class SIFSection for storing and generating section defintion
// Created by Liming Jiang @UoE


#include <SIFSection.h>
#include <ElasticSection3d.h>
#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h> 
#include <FiberSectionGJThermal.h> 
#include <NDMaterial.h>

#include <QuadPatch.h> 
#include <QuadCell.h> 
#include <UniaxialFiber3d.h> 
#include <UniaxialMaterial.h>
//#include <ConcreteSThermal.h>
#include <PlasticDamageConcretePlaneStressThermal.h>
#include <SteelECThermal.h>
#include <PlateRebarMaterialThermal.h>
#include <PlateFromPlaneStressMaterialThermal.h>
#include <LayeredShellFiberSectionThermal.h>
#include <PlateRebarMaterialThermal.h>

#include <ElasticPlateSection.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>
#include <MembranePlateFiberSectionThermal.h> 
#include <ElasticIsotropic3DThermal.h>
#include <ElasticIsotropicMaterial.h>
#include <PlateFiberMaterialThermal.h>



SIFSection::SIFSection(int tag, int SectionTypeTag, SIFMaterial *theSifMaterial):TaggedObject(tag), 
SectionTypeTag(SectionTypeTag) 
{
	//Section type:  Rectangular section (1), I section(2),Composite section (3), protected Isection(22);
  //Section type can also be section for slab, because slab is also defined with sifmember,then we call it a slab section.
  //Section type: Slab section (10)
  

	if (theSifMaterial!=0) {
    SIFMaterialPtr = theSifMaterial;
  } else {
    opserr<<"WARNING:: the SIFMaterial can't be assigned to SIFSection "
    <<this->getTag()<<endln;
  }

}

SIFSection::~SIFSection()
{
	//
}

int 
SIFSection::getSectionTypeTag()
{
	return SectionTypeTag;
}



int 
SIFSection::AssignSectionPars(const Vector& sectionPars)
{
  if (sectionPars==0) {
	  opserr<<"SIFSection::AssignSectionPars recived an empty sectionPars"<<endln;
	  return -1;
  }
    //if the section is rectangular;
    //2 SectionPars:  Breadth, Height;
   if (SectionTypeTag==1) {
     if (sectionPars.Size()!=2) {
        opserr<<"WARNING:: SIFSection "<<this->getTag()
        <<" wants 2 arguments for definition "
        <<endln;
        return -1;

     } else {
       theSectionPars=sectionPars;
     }
   }
    //if the section is I-section
    //4 SectionPars: Beam Height, Flange Breadth, web thickness, flange thickness
   else if(SectionTypeTag==2) {
      if (sectionPars.Size()!=4) {
         opserr<<"WARNING:: SIFSection "<<this->getTag()
         <<" wants 4 arguments for definition"
         <<endln;
         return -1;
          
       } else {
         theSectionPars=sectionPars;
       }
    
  }
     else if(SectionTypeTag==22) {
      if (sectionPars.Size()!=5) {
         opserr<<"WARNING:: SIFSection "<<this->getTag()
         <<" wants 5 arguments for definition"
         <<endln;
         return -1;
          
       } else {
         theSectionPars=sectionPars;
       }
    
  }
   //for composite beam
    else if(SectionTypeTag==3) {
      if (sectionPars.Size()!=6) {
         opserr<<"WARNING:: SIFSection "<<this->getTag()
         <<" wants 4 arguments for definition"
         <<endln;
         return -1;
          
       } else {
         theSectionPars=sectionPars;
       }
    
  }
  //if the section is Slab-section
  //1 SectionPar: slab thickness
   else if(SectionTypeTag==10) {
     if (sectionPars.Size()!=2) {
       opserr<<"WARNING:: SIFSection "<<this->getTag()
       <<" wants 1 arguments for definition"
       <<endln;
       return -1;
       
     } else {
       theSectionPars=sectionPars;
     }
     
   }
   else {
    opserr<<"WARNING:: the SIFMaterial can't be assigned to SIFSection "
    <<this->getTag()<<endln;
    return -1;
  }
  
  return 0;
	
}

const Vector& 
SIFSection::getSectionPars(void)
{
	//This function will return a Vector which describes the dimension of the section
  //This Vector varies with the section type
  
  //for rectangular section
  return theSectionPars;
	
}


//The following function is to assign SIFMaterial to the section
int
SIFSection::assignSecondMaterial(SIFMaterial* theSifMaterial)
{
  SecSIFMaterialPtr = theSifMaterial;
  
  return 0;
}

//The following function is to return the pointer to the SIFMaterial assigned to this section
SIFMaterial*
SIFSection::getSIFMaterialPtr(int MaterialTag )
{
	if (MaterialTag ==1)
		return SIFMaterialPtr;
	else if(MaterialTag ==2)
		return SecSIFMaterialPtr;
	else{
		opserr<<"WARNING::Unknown material tag when requesting SIFMaterial from SIFSection"<<endln;
		return 0;
	}
  
}


//ElasticSection3d(int tag, double E, double A, double Iz,  double Iy, double G, double J);
SectionForceDeformation* 
SIFSection::DefineBeamSection(bool isElastic)
{
	SectionForceDeformation* theSection =0;
	int sectionTag = this->getTag();

	//This should be extended to Fiber Section

	//theSection = new FiberSection3d( )

    UniaxialMaterial* theUniMaterial = SIFMaterialPtr->getUniaxialMaterial();
	UniaxialMaterial* theUniMaterial1 = SIFMaterialPtr->getUniaxialMaterial(isElastic);
	int matTag = theUniMaterial->getTag();
	//-------------------------------------Rectangular Section--------------------------------
	if(SectionTypeTag==1||SectionTypeTag==10){		// which means this is a rectangular section;
		if(isElastic){
			double E0= SIFMaterialPtr->getInitialModulus();
			double secb = theSectionPars(0);
			double secd = theSectionPars(1);
			double Area = secb* secd;
			double Modulus_z = secb*secd*secd*secd/12.0;
			double Modulus_y = secd*secb*secb*secb/12.0;

			theSection =  new ElasticSection3d(sectionTag, E0, Area, Modulus_z,  Modulus_y, 1e-5, 1e-5);
			return theSection;
		}

		
		double Breadth = theSectionPars(0);
		double Height = theSectionPars(1);
		int numSubdivY = 8;				// number of subdivisions (fibers) in the local y direction;
		int numSubdivZ = 4;				// number of subdivisions (fibers) in the local z direction;
		if(SectionTypeTag==10){
			Breadth = theSectionPars(1);
			Height = theSectionPars(0);
			numSubdivZ = 24;
		}

		static Matrix vertexCoords(4,2);
		int j;

		vertexCoords(0,0) = -Height/2;
		vertexCoords(0,1) = -Breadth/2;
		vertexCoords(2,0) = Height/2;
	    vertexCoords(2,1) = Breadth/2;
		vertexCoords(1,0) = vertexCoords(2,0);
		vertexCoords(1,1) = vertexCoords(0,1);	
		vertexCoords(3,0) = vertexCoords(0,0);
		vertexCoords(3,1) = vertexCoords(2,1);		

		// create patch
		QuadPatch *patch = new QuadPatch(matTag, numSubdivY, numSubdivZ, vertexCoords);

		int i, k = 1;
		int numFibers = patch->getNumCells();
        //opserr << "\nnumFibers: " << numFibers;
      
		static Vector fiberPosition(2);
      
		ID     fibersMaterial(1);
		Matrix fibersPosition(2,numFibers);
		double fibersArea;

		Cell **cell = patch->getCells();			// get all the cells from the rectangular patch
    	if (cell == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            return 0;
        }    

		Fiber **fiber = new Fiber *[numFibers];
		if (fiber == 0) {
		opserr <<  "WARNING unable to allocate fibers \n";
		return 0;
		}
      

        for (i = 0; i < numFibers; i++)
        { 
            fibersArea     = cell[i]->getArea();
            fiberPosition     = cell[i]->getCentroidPosition();
#ifdef _DEBUG
			//opserr<<"SIFSection::position: "<<fiberPosition;
#endif
			fiber[i] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
	        k++;
         }

		 // create 3d section      
		double GJ =2e11;
		 theSection =  new FiberSectionGJThermal(sectionTag, numFibers, fiber,GJ);
		 // Delete fibers
		for (i = 0; i < numFibers; i++)
			delete fiber[i];

	    // Delete fiber array
		delete [] fiber;
	 }
	//-----------------------------Isection---------------------------------------------------
	 else if(SectionTypeTag==2||SectionTypeTag==22){		// which means this is an I section;
		
		double d = theSectionPars(0);
		double tw = theSectionPars(1);
		double bf = theSectionPars(2);
		double tf = theSectionPars(3);
		double wd = d - 2*tf;				//web depth;

		int nfdw = 8;			//number of fibers along web depth;
		int nftw = 4;			//number of fibers along web thickness;
		int nfbf = 8;			//number of fibers along flange width;
		int nftf = 4;			//number of fibers along flange thickness;

		int numFibers = nfdw*nftw + nfbf*nftf*2;		//number of fibers for the whole I section;

		Fiber **fiber = new Fiber *[numFibers];
		if (fiber == 0) {
			opserr <<  "WARNING unable to allocate fibers \n";
			exit (-1);
		}

		// patch for web;
		static Matrix vertexCoords1(4,2);
		
		vertexCoords1(0,0) = -wd/2;
		vertexCoords1(0,1) = -tw/2;
		vertexCoords1(2,0) = wd/2;
	    vertexCoords1(2,1) = tw/2;
		vertexCoords1(1,0) = vertexCoords1(2,0);
		vertexCoords1(1,1) = vertexCoords1(0,1);	
		vertexCoords1(3,0) = vertexCoords1(0,0);
		vertexCoords1(3,1) = vertexCoords1(2,1);		

		QuadPatch *patch1 = new QuadPatch(matTag, nfdw, nftw, vertexCoords1);

		int i, k = 0;
		int numFibersWeb = patch1->getNumCells();
       // opserr << "\nnumFibers for web: " << numFibersWeb;
      
		static Vector fiberPosition(2);
      
		double fibersArea;

		Cell **cell1 = patch1->getCells();			// get all the cells from the web
    	if (cell1 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        }    

		// Delete patches
		delete patch1;

		
      
        for (i = 0; i < numFibersWeb; i++)
        { 
            fibersArea = cell1[i]->getArea();
            fiberPosition = cell1[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in web "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersWeb; i++)
			delete cell1[i];
		delete [] cell1;

		// patch for top flange;
		static Matrix vertexCoords2(4,2);
		
		vertexCoords2(0,0) = wd/2;
		vertexCoords2(0,1) = -bf/2;
		vertexCoords2(1,0) = d/2;
		vertexCoords2(1,1) = -bf/2;
		vertexCoords2(2,0) = d/2;
	    vertexCoords2(2,1) = bf/2;
		vertexCoords2(3,0) = wd/2;
		vertexCoords2(3,1) = bf/2;		

		QuadPatch *patch2 = new QuadPatch(matTag, nftf, nfbf, vertexCoords2);

		int numFibersTF = patch2->getNumCells();
        //opserr << "\nnumFibers for top flange: " << numFibersTF;
      
		Cell **cell2 = patch2->getCells();			// get all the cells from the web
    	if (cell2 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        } 

		// Delete patches
		delete patch2;
      
        for (i = 0; i < numFibersTF; i++)
        { 
            fibersArea = cell2[i]->getArea();
            fiberPosition = cell2[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in top flange "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersTF; i++)
			delete cell2[i];
		delete [] cell2;

		// patch for bottom flange;
		static Matrix vertexCoords3(4,2);
		
		vertexCoords3(0,0) = -d/2;
		vertexCoords3(0,1) = -bf/2;
		vertexCoords3(1,0) = -wd/2;
		vertexCoords3(1,1) = -bf/2;
		vertexCoords3(2,0) = -wd/2;
	    vertexCoords3(2,1) = bf/2;
		vertexCoords3(3,0) = -d/2;
		vertexCoords3(3,1) = bf/2;		

		QuadPatch *patch3 = new QuadPatch(matTag, nftf, nfbf, vertexCoords3);

		int numFibersBF = patch3->getNumCells();
       // opserr << "\nnumFibers for bottom flange: " << numFibersBF;
      
		Cell **cell3 = patch3->getCells();			// get all the cells from the web
    	if (cell3 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        }

		// Delete patches
		delete patch3;
      
        for (i = 0; i < numFibersBF; i++)
        { 
            fibersArea = cell3[i]->getArea();
            fiberPosition = cell3[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in bottom flange "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersBF; i++)
			delete cell3[i];
		delete [] cell3;

		 // create 3d section      
		 double GJ =2e11;
		 theSection =  new FiberSectionGJThermal(sectionTag, numFibers, fiber,GJ);
		 
		 // Delete fibers
		for (i = 0; i < numFibers; i++)
			delete fiber[i];

	    // Delete fiber array
		delete []fiber;

	 }
	//---------------------------------end of defining Isection-----------------------------------
	  else if(SectionTypeTag==3){		// which means this is an composite section;
		
		double d = theSectionPars(0);
		double tw = theSectionPars(1);
		double bf = theSectionPars(2);
		double tf = theSectionPars(3);
		double wd = d - 2*tf;				//web depth;
		double slabW = theSectionPars(5);
		double slabT = theSectionPars(6);

		int nfdw = 8;			//number of fibers along web depth;
		int nftw = 4;			//number of fibers along web thickness;
		int nfbf = 8;			//number of fibers along flange width;
		int nftf = 4;			//number of fibers along flange thickness;

		int numFibers = nfdw*nftw + nfbf*nftf*2 +10*20;		//number of fibers for the whole I section;

		Fiber **fiber = new Fiber *[numFibers];
		if (fiber == 0) {
			opserr <<  "WARNING unable to allocate fibers \n";
			exit (-1);
		}

		// patch for web;
		static Matrix vertexCoords1(4,2);
		
		vertexCoords1(0,0) = -wd/2;
		vertexCoords1(0,1) = -tw/2;
		vertexCoords1(2,0) = wd/2;
	    vertexCoords1(2,1) = tw/2;
		vertexCoords1(1,0) = vertexCoords1(2,0);
		vertexCoords1(1,1) = vertexCoords1(0,1);	
		vertexCoords1(3,0) = vertexCoords1(0,0);
		vertexCoords1(3,1) = vertexCoords1(2,1);		

		QuadPatch *patch1 = new QuadPatch(matTag, nfdw, nftw, vertexCoords1);

		int i, k = 0;
		int numFibersWeb = patch1->getNumCells();
       // opserr << "\nnumFibers for web: " << numFibersWeb;
      
		static Vector fiberPosition(2);
      
		double fibersArea;

		Cell **cell1 = patch1->getCells();			// get all the cells from the web
    	if (cell1 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        }    

		// Delete patches
		delete patch1;

		
      
        for (i = 0; i < numFibersWeb; i++)
        { 
            fibersArea = cell1[i]->getArea();
            fiberPosition = cell1[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in web "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersWeb; i++)
			delete cell1[i];
		delete [] cell1;

		// patch for top flange;
		static Matrix vertexCoords2(4,2);
		
		vertexCoords2(0,0) = wd/2;
		vertexCoords2(0,1) = -bf/2;
		vertexCoords2(1,0) = d/2;
		vertexCoords2(1,1) = -bf/2;
		vertexCoords2(2,0) = d/2;
	    vertexCoords2(2,1) = bf/2;
		vertexCoords2(3,0) = wd/2;
		vertexCoords2(3,1) = bf/2;		

		QuadPatch *patch2 = new QuadPatch(matTag, nftf, nfbf, vertexCoords2);

		int numFibersTF = patch2->getNumCells();
        //opserr << "\nnumFibers for top flange: " << numFibersTF;
      
		Cell **cell2 = patch2->getCells();			// get all the cells from the web
    	if (cell2 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        } 

		// Delete patches
		delete patch2;
      
        for (i = 0; i < numFibersTF; i++)
        { 
            fibersArea = cell2[i]->getArea();
            fiberPosition = cell2[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in top flange "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersTF; i++)
			delete cell2[i];
		delete [] cell2;

		// patch for bottom flange;
		static Matrix vertexCoords3(4,2);
		
		vertexCoords3(0,0) = -d/2;
		vertexCoords3(0,1) = -bf/2;
		vertexCoords3(1,0) = -wd/2;
		vertexCoords3(1,1) = -bf/2;
		vertexCoords3(2,0) = -wd/2;
	    vertexCoords3(2,1) = bf/2;
		vertexCoords3(3,0) = -d/2;
		vertexCoords3(3,1) = bf/2;		

		QuadPatch *patch3 = new QuadPatch(matTag, nftf, nfbf, vertexCoords3);

		int numFibersBF = patch3->getNumCells();
       // opserr << "\nnumFibers for bottom flange: " << numFibersBF;
      
		Cell **cell3 = patch3->getCells();			// get all the cells from the web
    	if (cell3 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        }

		// Delete patches
		delete patch3;
      
        for (i = 0; i < numFibersBF; i++)
        { 
            fibersArea = cell3[i]->getArea();
            fiberPosition = cell3[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in bottom flange "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersBF; i++)
			delete cell3[i];
		delete [] cell3;

		//FOR CONCRETE SLAB
		static Matrix vertexCoords4(4,2);
		
		vertexCoords4(0,0) = d/2;
		vertexCoords4(0,1) = -slabW/2;
		vertexCoords4(1,0) = d/2+slabT;
		vertexCoords4(1,1) = -slabW/2;
		vertexCoords4(2,0) = d/2+slabT;
	    vertexCoords4(2,1) = slabW/2;
		vertexCoords4(3,0) = d/2;
		vertexCoords4(3,1) = slabW/2;		

		QuadPatch *patch4 = new QuadPatch(matTag, 10, 20, vertexCoords2);

	//	int numFibersTF = patch4->getNumCells();
        //opserr << "\nnumFibers for top flange: " << numFibersTF;
      
		Cell **cell4 = patch4->getCells();			// get all the cells from the web
    	if (cell4 == 0){
			opserr <<  "WARNING out of run to create fibers\n";
            exit (-1);
        } 

		// Delete patches
		delete patch4;
      
        for (i = 0; i < numFibersTF; i++)
        { 
            fibersArea = cell4[i]->getArea();
            fiberPosition = cell4[i]->getCentroidPosition();
			//opserr<<"SIFSection::position: - fibers in top flange "<<fiberPosition;
			fiber[k] = new UniaxialFiber3d(k, *theUniMaterial, fibersArea, fiberPosition);
			k++;
         }

		// Delete cells
		for (i = 0; i < numFibersTF; i++)
			delete cell4[i];
		delete [] cell4;

		 // create 3d section      
		 double GJ =20000;
		 theSection =  new FiberSectionGJThermal(sectionTag, numFibers, fiber,GJ);
		 
		 // Delete fibers
		for (i = 0; i < numFibers; i++)
			delete fiber[i];

	    // Delete fiber array
		delete []fiber;

	 }
	 //end of defining composite section
	 return  theSection;
}


SectionForceDeformation* 
SIFSection::DefineShellSection(bool isElastic)
{
	SectionForceDeformation* theShellSection=0;
	 double thick = theSectionPars(0);
	int tag, nLayers, matTag;
	double h;
	double *thickness;
	NDMaterial **theMats;

	if(isElastic){
		double thick = theSectionPars(0);
		ElasticIsotropic3DThermal* theMaterial3D = new ElasticIsotropic3DThermal(11, 1.92e10, 0.3, 0, 1.45e-5, true);
		//PlateFiberMaterialThermal(   int    tag, NDMaterial &the3DMaterial ) ;
		PlateFiberMaterialThermal* thePlaterMaterial = new PlateFiberMaterialThermal(12, *theMaterial3D);
		
		nLayers =10;
        
	    theMats   = new NDMaterial*[nLayers];
        thickness = new double[nLayers];
		
		for(int i=0; i<nLayers; i++){
			theMats[i] = thePlaterMaterial ;
			thickness[i]= thick/nLayers;		

		}

		theShellSection = new LayeredShellFiberSectionThermal(this->getTag(), 10, thickness, theMats);
      if (thickness != 0) delete thickness;
      if (theMats != 0) delete [] theMats;

  }
  else{
		double thick = theSectionPars(0);
		
		//ConcreteSThermal::ConcreteSThermal(int tag, double rE, double rnu, double rfc, double rft, double rEs) :
		//NDMaterial* theNDConcrete = new ConcreteSThermal (11, 1.92e10,0.3,30e6, 3e6, 1.92e9);
		double ft = 3e6;
		double fc = 30e6;
		double E = 1.92e10;
		double gt = ft*ft / E*4.0;
		double gc = fc*fc / E * 6;

		NDMaterial* theNDConcrete = new  PlasticDamageConcretePlaneStressThermal(11, E, 0.2, ft, fc,gt,gc,0.4, 4.0, 0.4,0.1);
		//PlateFromPlaneStressMaterialThermal::PlateFromPlaneStressMaterialThermal(    
				   //int tag, NDMaterial &ndMat, double g )
		PlateFromPlaneStressMaterialThermal* thePlateMatCon = new PlateFromPlaneStressMaterialThermal (12, *theNDConcrete, 1.2e10);
		//ElasticIsotropic3DThermal* theMaterial3D = new ElasticIsotropic3DThermal(11, 1.92e10, 0.3, 0, 1.45e-5, 2);
		//PlateFiberMaterialThermal(   int    tag, NDMaterial &the3DMaterial ) ;
		//PlateFiberMaterialThermal* thePlateMatCon = new PlateFiberMaterialThermal(12, *theMaterial3D);
		//SteelECThermal(int tag, int TypeTag , double FY, double E)
		UniaxialMaterial* theRebarMat = new SteelECThermal(13, 22, 450e6, 2.06e11);
		
		//PlateRebarMaterialThermal(int tag, UniaxialMaterial &uniMat, double ang)
		NDMaterial* thePlateRebar0 = new PlateRebarMaterialThermal(20, *theRebarMat,0);
		NDMaterial* thePlateRebar90 = new PlateRebarMaterialThermal(21, *theRebarMat,90);

		//LayeredShellFiberSectionThermal::LayeredShellFiberSectionThermal(    
          //                         int tag, 
            //                       int iLayers, 
              //                     double *thickness, 
                //                   NDMaterial **fibers )
			
			nLayers =14;
        
	    theMats   = new NDMaterial*[nLayers];
        thickness = new double[nLayers];
		double rebarT0 = 1.98e-4;
		double rebarT90 = 1.98e-4*2;
		for(int i=0; i<nLayers; i++){
			if(i==2||i==nLayers-3){
				theMats[i] = thePlateRebar0;
				thickness[i]= rebarT0;
			}
			else if(i==3||i==nLayers-4){
				theMats[i] = thePlateRebar90;
				thickness[i]= rebarT90;
			}
			else if(i==4||i==nLayers-5){
				theMats[i] = thePlateMatCon;
				thickness[i]= thick/(nLayers-4) - rebarT0- rebarT90;
			}
			else {
				theMats[i] = thePlateMatCon;
				thickness[i]= thick/(nLayers-4);		
			}

		}

		theShellSection = new LayeredShellFiberSectionThermal(this->getTag(), 14, thickness, theMats);
      if (thickness != 0) delete thickness;
      if (theMats != 0) delete [] theMats;

		//end of layer section

		//theShellSection = new MembranePlateFiberSectionThermal(10, thickness, *thePlaterMaterial);
		
    
  }


return theShellSection;

}