#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees


wipe;

HeatTransfer 2D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1.
HTMaterial ConcreteEC2 2 0.0; 

#HTEntity Isection $tag $centreX $centreY $Bf $Hb  $Tw  $Tf;
HTEntity  Isection   1    0.0   0.15    0.3 0.30 0.02 0.02;
#HTEntity Block $tag $centreX $centreY $slabB $slabH;
HTEntity  Block   2    0.0      0.35     0.7    0.1;


#HTMesh $MeshTag $EntityTag  $MaterialTag -SecondMat 2
HTMesh 1 1 1 -phaseChange 0 -MeshCtrls 0.01 0.005 0.005 0.014 
HTMesh 2 2 1 -phaseChange 1 -MeshCtrls 0.02 0.02 

HTRefineMesh -Entity 2 -SeedTag 1 4 -space 0.02 10 0.01 9 0.005 4 0.01 9 0.02 10;
HTMeshAll;

#RefineMesh 1 -seed 1 ;
SetInitialT 293.15;

HTNodeSet 1 -Entity 1 -face 12;
HTNodeSet 2 -Entity 2 -face 1 -locx -0.1 0.1;
HTCoupleT -NodeSet 1 2;

HTConstants 1 4.0 293.15 0.7 5.67e-8 0.7 ;
HTConstants 2 25.0 293.15 0.7 5.67e-8 0.7;

HTPattern AmbientBC 1 {
	HeatFluxBC -HTEntity 2 -faceTag 4 -type ConvecAndRad -HTConstants 1;	
}

FireModel standard 1; 

HTNodeSet 3 -Entity 2 -Locx -0.3 -0.1;
HTEleSet 1 -Entity 2 -NodeSet 3 -face 1;
HTNodeSet 4 -Entity 2 -Locx 0.1 0.3;
HTEleSet 2 -Entity 2 -NodeSet 4 -face 1;

HTPattern fire 2 model 1 {
	HeatFluxBC -HTEntity 1 -face 1 4 5 6 7 8 9 -type ConvecAndRad -HTConstants 2;	
	HeatFluxBC -HTEleSet 1 -face 1 -type ConvecAndRad -HTConstants 2;
	HeatFluxBC -HTEleSet 2 -face 1 -type ConvecAndRad -HTConstants 2;	
}

#Recorder
HTRecorder -file temp0.out -NodeSet 1;
HTRecorder -file temp1.out -NodeSet 2;

HTAnalysis HeatTransfer
HTAnalyze 20 30;	

wipeHT;


