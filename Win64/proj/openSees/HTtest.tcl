
#This Tcl file is written for testing tcl commands for Heat Transfer module in OpenSees


wipe;
HeatTransfer 3D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1;
HTMaterial ConcreteEC2 2 0.015; 

#Secondary Beam section with slab
# beam sections: (I section here)
set slabw 0.02;  #slab width
set slabl 0.02;	#slab length
set slabt 0.1;	#slab thickness

HTEntity Line 1 0.0 $slabt; #To create a brick whose y-axis is through depth whith its centroid at 0,0,0
HTEntity Block 2 0.0 0.0 $slabw $slabt;
HTEntity Brick 3 0.0 0.0 0.0 $slabw $slabt $slabl;

#secondary beam mesh components
set elex [expr $slabw/1];
set eley [expr $slabt/100];
set elez [expr $slabl/1]

#Secondary beam and mesh
#HTMesh $MeshTag $EntityTag  $MaterialTag(first material tag) -SecondMat $secondmaterialTag, -phaseChange, 4FisrtMaterial (no need for steel), 4SecondMaterial

#HTMesh 1 1 2 -phaseChange 1 -MeshCtrls $eley;  

#HTMesh 1 2 2 -phaseChange 1 -MeshCtrls $elex $eley;  

HTMesh 1 3 2 -phaseChange 1 -MeshCtrls $elex $eley $elez;  

HTMeshAll;

puts "Mesh complete"

SetInitialT 293.15;

# HTconstants tag, convection coefficient,ambient temp.,emissity,boltzman,absorptivity, irradiation
HTConstants 1 4.0 293.15 0.8 0.8; # concrete
HTConstants 2 25.0 293.15 0.8 0.8;

HTPattern AmbientBC 1 {
	HeatFluxBC -HTEntity 3 -face 2 -type ConvecAndRad -HTConstants 1;  #for line
	#HeatFluxBC -HTEntity 2 -face 4 -type ConvecAndRad -HTConstants 1; 
	#HeatFluxBC -HTEntity 1 -face 1 -type -ConvecAndRad -HTConstants 2; #secondary concrete beam
	# 16 is the top surface of concrete slab
}

#puts "first entity done"

FireModel standard 1; #-duration 3600; 
FireModel UserDefined 3 -file fire.dat; #-duration 3600; 
FireModel Localised 2 -origin 0.0 -3.0 0.0 -firePars 1.0 2E6 3.0 2 ;
#FireModel hydroCarbon 1; 

#puts "fire1" 

# fire for secondary concrete beam and slab
HTPattern fire 2 model 1 {
	#HeatFluxBC -HTEntity 1 -face 1 -type -Prescribed -HTConstants 2; #secondary concrete beam
	#HeatFluxBC -HTEntity 1 -face 1 -type -ConvecAndRad -HTConstants 2;
	HeatFluxBC -HTEntity 3 -face 1 -type -ConvecAndRad -HTConstants 2;
}


#puts "fire3"

#Record concrete beam temperatures
#HTRecorder -file temp2.out -xloc 0.0 -yloc [expr $cd/18] [expr 3*$cd/18] [expr 5*$cd/18] [expr 7*$cd/18] [expr 9*$cd/18] [expr 11*$cd/18] [expr 13*$cd/18] [expr 15*$cd/18] [expr 17*$cd/18];


# End of building the heat transfer models;
#HTNodeSet 1 -Entity 1 -Face 1  ;
HTRecorder -file temp.out Node 1 51 101 ;
#HTNodeSet 1 -Entity 1 -Face 1 -Locx 0.0 -Locy -0.05 -Locz 0.0;

	
HTAnalysis HeatTransfer
#HTAnalyze $numSteps $timeStep;
HTAnalyze 360 30;	

wipe;

	
	