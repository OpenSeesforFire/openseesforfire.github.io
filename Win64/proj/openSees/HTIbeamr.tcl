#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees

#Dimension of beam
set BF 0.20;
set HB 0.40;
set Tw 0.013;
set Tf 0.015;
set Len 6.0;
set Centrex 0.0
set Centrey [expr $HB/2.0]
set Centrez [expr $Len/2.0]

set elex [expr ($BF-$Tw)/8.0]
set eleweby [expr ($HB-$Tf*2.0)/10.0]
set eley [expr $Tf/2.0]
set elewebx [expr $Tw/2.0]
set elez [expr $Len/120.0]


wipe;

HeatTransfer 3D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1.
#HTMaterial ConcreteEC2 2 0.0; 

#HTEntity Isection3d $tag $centreX $centreY $centreZ $Bf    $Hb   $Tw  $Tf $Len;
HTEntity  Isection3D   1    $Centrex   $Centrey  $Centrez   $BF   $HB  $Tw    $Tf   $Len;
#HTEntity Block $tag $centreX $centreY $slabB $slabH;


#HTMesh $MeshTag $EntityTag  $MaterialTag -SecondMat 2
HTMesh 1 1 1 -phaseChange 0 -MeshCtrls $elex $eley $elez $elewebx $eleweby
#HTMesh 2 2 1 -phaseChange 1 -MeshCtrls 0.02 0.02 

HTMeshAll;
puts "mesh done";

#RefineMesh 1 -seed 1 ;
SetInitialT 289.15;


HTConstants 1 0.1 293.15 0.85 0.85 ;
HTConstants 2 15.0 293.15 0.85 0.85;

HTPattern AmbientBC 1 {
	#HeatFluxBC -HTEntity 1 -face  1 4 5 6 7 8 9  -type -ConvecAndRad -HTConstants 2;
}

#FireModel standard 1; 
FireModel Localised 2 -origin 0.0 -2.0 $Centrez -firePars 1.0 1127.0E3 2.4 2 ;
FireModel LocalisedSFPE 3 -origin 0.0 -2.0 $Centrez -firePars 1.0 1127.0E3 2.40 2.0 2 ;
FireModel LocalisedSFPE 4 -origin 0.0 -1.0 $Centrez -firePars 0.5 200E3 1.15 1.00 2 ;
FireModel NaturalFire 5 -firePars 1.0 569.0E3 2.40 2 231.3;
SetFirePars firemodel 5 0.0 -2.0 $Centrez 569.0E3 1.0 231.3 6E3;


#FireModel idealised 5 -origin 0.0 -0.6 $Centrez -Q 5000 -Linear 2.0 0.5;
puts "Fire done";
HTPattern fire 2 model 5 {
	HeatFluxBC -HTEntity 1 -face 1  -type -Prescribed -HTConstants 2 -par 4; #Fire action using input heat flux
	HeatFluxBC -HTEntity 1 -face 4 5  -type -Prescribed -HTConstants 2 -par 3;
	HeatFluxBC -HTEntity 1 -face 6 7  -type -Prescribed -HTConstants 2 -par 2;
	HeatFluxBC -HTEntity 1 -face 8 9  -type -Prescribed -HTConstants 2 -par 1;
	HeatFluxBC -HTEntity 1 -face  1 4 5 6 7 8 9  -type -ConvecAndRad -HTConstants 2;

}

#Recorder
set flloc [expr ($BF-$Tw)/4.0+$Tw/2]
set webloc [expr  $HB/2]
HTNodeSet 10 -Locx $flloc -Locy 0.0075;
HTNodeSet 11 -Locx 0.0 -Locy $webloc;
HTNodeSet 12 -Locx $flloc -Locy 0.3925;

HTRecorder -file tempBLF.out -NodeSet 10;    #Nodes of beam on the beam-slab interface 
HTRecorder -file tempWeb.out -NodeSet 11;
HTRecorder -file tempBUF.out -NodeSet 12;
#HTRecorder -file temp1.out -NodeSet 2;    #Nodes of slab on the beam-slab interface 
#HTRecorder -file temp2.out -NodeSet 5;    #All the nodes along x=0.0

HTAnalysis HeatTransfer
HTAnalyze 72 15;	
#HTAnalyze 80 15;
wipeHT;


