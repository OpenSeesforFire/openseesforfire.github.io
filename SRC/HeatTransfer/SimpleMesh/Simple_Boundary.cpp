#include <Simple_Boundary.h>
#include <StandardStream.h>
#include <math.h>

//StandardStream sserr;
//OPS_Stream* opserrPtr = &sserr;




Simple_Boundary::Simple_Boundary(int tag,HeatTransferDomain* theDomain)

:BoundaryTag(tag),theHTDomain(theDomain),isHTDomain(true)

{

}

Simple_Boundary::Simple_Boundary(int tag,Domain* theDomain)

:BoundaryTag(tag),theDomain(theDomain),isHTDomain(false)

{

}


Simple_Boundary::~Simple_Boundary()
{

}









void 

Simple_Boundary::GeneratingSP_BC(const ID& NodesRange,int nDoF, double value, bool ISconstant, int LoadPatternTag)

{

	int NumSP_BC = NodesRange.Size();
	if (NumSP_BC==0)
	{
	    opserr<<"Warning:No SP_BC will  be generated in Simple_Boundary "<<this->getTag()<<endln;
	}
	else
	{
		for(int i= 0;i<NumSP_BC; i++){
			if(isHTDomain)
			{
				TemperatureBC* Temp_BC = new TemperatureBC(NodesRange(i),nDoF,value,ISconstant);
				theHTDomain->addTemperatureBC(Temp_BC);
			}
			else
			{
				SP_Constraint* SP_BC = new SP_Constraint(NodesRange(i),nDoF,value,ISconstant);	

				if(LoadPatternTag == -1)
					theDomain->addSP_Constraint(SP_BC);
				else
					theDomain->addSP_Constraint(SP_BC,LoadPatternTag);
			}
		}
	}
}



void 
Simple_Boundary::GeneratingMP_BC(const ID& NodesRange, const ID& NodesRange1, int DoFTag0,int DoFTag1)
{
    ID EqualDoF;
	if(DoFTag1==-1)
	{
		ID TempID(1);
		TempID(0)=DoFTag0;
        EqualDoF = TempID;
	}
	else
	{
		ID TempID(2);
		TempID(0)=DoFTag0;
		TempID(1)=DoFTag1;
		EqualDoF = TempID;
	}
	Matrix Ccr(EqualDoF.Size(),EqualDoF.Size()); 
	Ccr.Zero();
	for (int j = 0; j < EqualDoF.Size(); j++){    
		Ccr (j,j) = 1.0;
	}
	
	if (NodesRange.Size()==0||NodesRange1.Size()==0)
		opserr<<"Waring:No MP_BC will  be generated in Simple_Boundary "<<this->getTag()<<endln;
    
	else if (NodesRange.Size()!=NodesRange1.Size())
	{
		opserr<<"Waning:The size of Node vectors for generatingMP_BC in Simple_Boundary "<<this->getTag()
			  <<" is not compatible"<<endln;
		exit(-9);
	}
	else
	{
		for (int i =0; i<NodesRange.Size();i++) {
			if(isHTDomain)
			{
		        //Heat Transfer has only one DOF
				MP_TemperatureBC* MP_TempBC = new MP_TemperatureBC(NodesRange(i),NodesRange1(i),Ccr,EqualDoF,EqualDoF);
				if(MP_TempBC==0)
					opserr << "Warning: ran out of memory for MP_TemperatureBC"<<endln;
				if(theHTDomain->addMP_TemperatureBC(MP_TempBC) ==false) 
				{
					opserr<< "Warning: could not add MP_TemperatureBC to domain" <<endln;
					delete MP_TempBC;
				}
			}
			else
			{
				MP_Constraint *theMP = new MP_Constraint (NodesRange(i),NodesRange1(i), Ccr, EqualDoF, EqualDoF);
				if (theMP == 0) 
					opserr << "WARNING ran out of memory for equalDOF MP_Constraint ";
				if (theDomain->addMP_Constraint (theMP) == false) 
				{
					opserr << "WARNING could not add equalDOF MP_Constraint to domain ";
					delete theMP;
				}

			}

		}

	}

}





void Simple_Boundary::GeneratingHeatFluxBC(const ID& ElesRange, int eleFaceTag, int HeatFluxTypeTag,int PatternTag,const Vector& HeatFluxConstants, int FireType)

{
    BoundaryPattern* thePattern = theHTDomain->getBoundaryPattern(PatternTag);
	int NumHeatFluxBCs = ElesRange.Size();
	if (NumHeatFluxBCs==0)
	{
	    opserr<<"Error: No HeatFluxBc will be defined in Simple_Boundary"<<this->getTag()<<endln;
	}

	int ExistingHeatFluxBCs = thePattern->getNumHeatFluxBCs();

   // HeatFluxConstants is a Vector transfering the cnstants:  h, Ta, sigma, epsilon, alpha, qir

	
	//1:Convection, 2:Radiation
	
	if(HeatFluxTypeTag==1)
	{
        for(int i= 0;i<NumHeatFluxBCs; i++){
			Convection* Convec_BC = new Convection(ExistingHeatFluxBCs+i,ElesRange(i),eleFaceTag,HeatFluxConstants(0),HeatFluxConstants(1));
			theHTDomain->addHeatFluxBC(Convec_BC,PatternTag);
		}
	}
	else if(HeatFluxTypeTag==2)
	{

		for(int i= 0;i<NumHeatFluxBCs; i++){

			Radiation* Radiat_BC = new Radiation(ExistingHeatFluxBCs+i,ElesRange(i),eleFaceTag,HeatFluxConstants(2),HeatFluxConstants(3),HeatFluxConstants(4),HeatFluxConstants(5));
			theHTDomain->addHeatFluxBC(Radiat_BC,PatternTag);

		}
	}
	else if(HeatFluxTypeTag ==3)
	{
		int NumNPF = (theHTDomain->getElement(ElesRange(0)))->getNumNodesperFace();

		for(int i= 0;i<NumHeatFluxBCs; i++){
			PrescribedSurfFlux* PreSurfFlux_BC = new PrescribedSurfFlux(ExistingHeatFluxBCs+i,ElesRange(i),eleFaceTag,NumNPF,FireType);
			theHTDomain->addHeatFluxBC(PreSurfFlux_BC,PatternTag);
		}	
	}
	else 
	{
		opserr<<"WARNING: Invalid HeatFluxTypeTag in simple_Boundary "<<this->getTag() <<endln;
	}



}







int Simple_Boundary::getTag()

{
     return BoundaryTag;
}
