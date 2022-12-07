#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees

HeatTransfer 2D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1.
HTMaterial GenericMaterial 2 250 1000 0.09;
#HTMaterial ConcreteEC2 2 0.015; 
#HTMaterial GenericMaterial 3 7850 600 54;  #Defining HeatTransfer Material with Material tag 1.
#HTMaterial GenericMaterial 4 2300 900 1.4;

# beam sections: (I section here)  units: m is required in HT model
set d 355e-3;	# depth
set bf 171.5e-3;	# flange width
set tf 11.5e-3;	# flange thickness
set tw 7.4e-3;	# web thickness
set coat 0.025;
#set slabh 0.3;

set halfD [expr $d/2];
#set halfB [expr $bf/2]

#HTEntity Block 1 0.2 0.2 0.4 0.4;  # Tag, $ x-coord of section centre, $ y-coord of section centre, $width of block, $ height of block  
                                    # surface definition of block: 1,2,3 and 4 is lower, left, right and upper surface respectively. 
#HTEntity Composite2D $centreX $centreY $HB $BF $TF $TW SlabW SlabT;
HTEntity ProtectedIsection 2 0.0 [expr $d/2] $bf $d $tw $tf $coat ;    # $halfB and $halfD is the position of centre point of steel section, so we can know the local coordinate for CB section locates at the bottom flange.

set elex [expr ($bf-$tw-2*$coat)/10];
set eley [expr ($tf)/4];
set elewx [expr ($tw)/2];
set elewy [expr ($d-2*$tf-2*$coat)/20];
#set elewy [expr ($d-2*$tf)/40];
set eleCoat [expr $coat/5];
#set eleSlaby [expr $slabh/32];
#set eleSlaby [expr $slabh/20];

#HTMesh $MeshTag $EntityTag  $MaterialTag(first material tag) -SecondMat $secondmaterialTag, -phaseChange, 4FisrtMaterial (no need for steel), 4SecondMaterial
#HTMesh 1 1 1 -phaseChange 0 -MeshCtrls 0.02 0.02;   #-origin 0.02
#HTMesh 2 2 1 -phaseChange 1 -MeshCtrls $elex $eley $elewx $elewy $eleSlabx $eleSlaby;   #-origin 0.02
HTMesh 2 2 1 -SecondMat 2 -phaseChange 0 0 -MeshCtrls $elex $eley $elewx $elewy $eleCoat;   #-origin 0.02 

HTMeshAll;

puts "Mesh"

set SetInitialT 295.15;

# HTconstants tag, convection coefficient,ambient temp.,emissity,boltzman,absorptivity, irradiation
#HTConstants 1 4.0 293.15 0.7 5.67e-8 0.7 418.738; # concrete
#HTConstants 2 25.0 293.15 0.9 5.67e-8 0.9 418.738; # steel
#HTConstants 3 4.0  $SetInitialT 0.7 5.67e-8 0.7 418.738; # concrete
#HTConstants 4 25.0 $SetInitialT 0.7 5.67e-8 0.7 418.738; # steel hydrocarbon
HTConstants 3 4.0  $SetInitialT 0.7 0.7; # concrete
HTConstants 4 25.0 $SetInitialT 0.7 0.7; # steel hydrocarbon
#HTConstants 5 50.0 293.15 0.7 5.67e-8 0.7 418.738; # concrete hydrocarbon


HTPattern AmbientBC 1 {
	#HeatFluxBC eleset 1 -faceTag 3 -type -ConvecAndRad -HTConstants 1;
	HeatFluxBC -HTEntity 2 -face 12 -type -ConvecAndRad -HTConstants 3; # 16 is the top surface of concrete
}


FireModel Standard 1; 



#HTEleSet 2 -HTEntity 1 -face 1;
HTPattern fire 2 model 1 {
	#HeatFluxBC eleset 2 -faceTag 1 -type -ConvecAndRad -HTConstants 2;
	HeatFluxBC -HTEntity 2  -face 1 4 5 6 7 8 9 -type -ConvecAndRad -HTConstants 4;
}


  HTNodeSet 10 -Locx 0.0  -Locy $tf  [expr $d-$tf] 
  HTNodeSet 11 -Locx [expr -0.5*$bf]  [expr 0.5*$bf]  -Locy [expr 0.5*$tf]
  HTNodeSet 12 -Locx [expr -0.5*$bf]  [expr 0.5*$bf]  -Locy [expr $d-0.5*$tf]
  file mkdir CardingtonStandardBeam1F
  HTRecorder -file CardingtonStandardBeam1F/Web.out -NodeSet 10
  HTRecorder -file CardingtonStandardBeam1F/Bottom.out -NodeSet 11
  HTRecorder -file CardingtonStandardBeam1F/Top.out -NodeSet 12

  HTAnalysis HeatTransfer
  HTAnalyze 172 30;	

wipeHT;


