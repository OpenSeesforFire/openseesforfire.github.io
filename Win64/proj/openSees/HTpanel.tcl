
#This Tcl file is written for testing tcl commands for Heat Transfer module in OpenSees


wipe;
HeatTransfer 3D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1;
HTMaterial ConcreteEC2 2 0.0; 
HTMaterial GenericMaterial 3 500 800 1.2；  #$density  $specific heat $ conductivity

#Secondary Beam section with slab
# beam sections: (I section here)
set slabw 3.15;  #slab width
set slabl 4.15;	#slab length
set slabt 0.1;	#slab thickness



HTEntity Brick 1 0.0 0.0 0.0 $slabw $slabt $slabl ; #To create a brick whose y-axis is through depth with its centroid at 0,0,0
#HTEntity Block 2 0.0 0.0 $slabl $slabt;

#secondary beam mesh components
set elex [expr $slabw/20];
set eley [expr $slabt/16];
set elez [expr $slabl/20]

#Secondary beam and mesh
#HTMesh $MeshTag $EntityTag  $MaterialTag(first material tag) -SecondMat $secondmaterialTag, -phaseChange, 4FisrtMaterial (no need for steel), 4SecondMaterial

HTMesh 1 1 2 -phaseChange 1 -MeshCtrls $elex $eley $elez;  

HTMeshAll;

puts "Mesh complete"

SetInitialT 293.15;

# HTconstants tag, convection coefficient,ambient temp.,emissity,absorptivity
HTConstants 1 4.0 293.15 0.7 0.7 ; # concrete
HTConstants 2 35.0 293.15 0.7 0.7;

HTPattern AmbientBC 1 {
	HeatFluxBC -HTEntity 1 -face 2 -type ConvecAndRad -HTConstants 1; 
	#HeatFluxBC -HTEntity 1 -face 1 -type -ConvecAndRad -HTConstants 2; #secondary concrete beam
	# 16 is the top surface of concrete slab
}

#puts "first entity done"

# FireModel standard 1; 
 FireModel Localised 2 -origin  0 -2.0  0 -firePars 1.0 1127E3 2.4 2 ;
# FireModel LocalisedSFPE 3 -origin 0.0 -2.0 $Centrez -firePars 1.0 1127.0E3 2.40 2.0 2 ;
# FireModel LocalisedSFPE 4 -origin 0.0 -1.0 $Centrez -firePars 0.5 200E3 1.15 1.00 2 ;
# Fire model with smoke modification
FireModel NaturalFire 5 -firePars 1.0 1127.0E3 2.40 2 230; # Diameter, HRR, H, symetrical axis, smokeT
#Set(overwrite fire parameters)
SetFirePars firemodel 5 0.0 -2.0 0.0 1127.0E3 1.0 327 12E3;   #locx, locy,locz, HRr,D, SmokeT, additional Heat flux (smoke induced)
FireModel	UserDefined 6	-file hf.dat -type 4; 


# fire for secondary concrete beam and slab
HTPattern fire 2 model 6 {
	#HeatFluxBC -HTEntity 1 -face 1 -type -Prescribed -HTConstants 2; #secondary concrete beam
	HeatFluxBC -HTEntity 1 -face 1 -type -Prescribed -HTConstants 2;
	HeatFluxBC -HTEntity 1 -face 1 -type ConvecAndRad -HTConstants 2; 
}


#Record concrete beam temperatures
#HTRecorder -file temp2.out -xloc 0.0 -yloc [expr $cd/18] [expr 3*$cd/18] [expr 5*$cd/18] [expr 7*$cd/18] [expr 9*$cd/18] [expr 11*$cd/18] [expr 13*$cd/18] [expr 15*$cd/18] [expr 17*$cd/18];



# End of building the heat transfer models;
HTNodeSet 1 -Locx 0.0 -Locz 0.0;
HTNodeSet 2 -Locx [expr $slabw/4] -Locz [expr $slabl/4];
HTNodeSet 3 -Locx [expr $slabw/2] -Locz [expr $slabl/2];
# HTNodeSet 4 -Entity 1 -face 1;
# HTNodeSet 5 -Locy -0.025;
# HTNodeSet 6 -Entity 1 -face 2;

HTRecorder -file slab01.out -NodeSet 1;
HTRecorder -file slab02.out -NodeSet 2;
HTRecorder -file slab03.out -NodeSet 3;
# HTRecorder -file bottom.out -NodeSet 4;
# HTRecorder -file bar.out -NodeSet 5;
# HTRecorder -file upper.out -NodeSet 6;


#HTNodeSet 1 -Entity 1 -Face 1 -Locx 0.0 -Locy -0.05 -Locz 0.0;

	
HTAnalysis HeatTransfer
#HTAnalyze $numSteps $timeStep;
HTAnalyze 120 15;	

wipe;

	
	