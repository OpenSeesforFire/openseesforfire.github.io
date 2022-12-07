#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees

#Dimension
set slabt 0.001;
set eley [expr $slabt/10]


wipe;

HeatTransfer 3D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

#HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1.
#HTMaterial ConcreteEC2 2 0.0; 
HTMaterial GenericMaterial 1 10 10 1000;
#HTEntity Isection3d $tag $centreX $centreY $centreZ $Bf    $Hb   $Tw  $Tf $Len;
set NumSec 28;
set IterNum [expr $NumSec+1];
for {set i 0} {$i<$IterNum} {incr i} {
set Tag [expr $i+1];
HTEntity  Line $Tag 0.0 $slabt;

}

#HTEntity Block $tag $centreX $centreY $slabB $slabH;


#HTMesh $MeshTag $EntityTag  $MaterialTag -SecondMat 2
for {set i 0} {$i<$IterNum} {incr i} {
set Tag [expr $i+1];
set Loc [expr 2.333+$i*(18.666-2.333)/$NumSec]
HTMesh $Tag $Tag 1 -phaseChange 0 -SectionLoc 0.0 $Loc -MeshCtrls $eley; 
}


 
#HTMesh 2 2 1 -phaseChange 1 -MeshCtrls 0.02 0.02 

HTMeshAll;
puts "mesh done";

#RefineMesh 1 -seed 1 ;
SetInitialT 293.15;


HTConstants 1 0.1 293.15 0.85 0.85 ;
HTConstants 2 35.0 293.15 0.85 0.85;

#HTPattern AmbientBC 1 {
#	HeatFluxBC -HTEntity 1 -face  1 4 5 6 7 8 9  -type -ConvecAndRad -HTConstants 2;
#}

#FireModel standard 1; 
#FireModel Localised 2 -origin 0.0 -2.0 $Centrez -firePars 1.0 1127.0E3 2.4 2 ;
#FireModel LocalisedSFPE 3 -origin 0.0 -2.0 $Centrez -firePars 1.0 1127.0E3 2.40 2.0 2 ;
#FireModel LocalisedSFPE 4 -origin 0.0 -1.0 $Centrez -firePars 0.5 200E3 1.15 1.00 2 ;
FireModel NaturalFire 5 -firePars 1.0 600E3 2.42 2 122.9;
#SetFirePars firemodel 5 0.0 -2.7 4.3 600E3 1.2 158.0 5E3;

set q 2351E3;
set d [expr sqrt($q/400.0e3/3.1415926)*2]
puts $d;
SetFirePars firemodel 5 0.0 -2.42 10.5 $q $d 424.0 25.29E3;

#FireModel idealised 5 -origin 0.0 -0.6 $Centrez -Q 5000 -Linear 2.0 0.5;
puts "Fire done";
HTPattern fire 2 model 5 {

	for {set i 0} {$i<$IterNum} {incr i} {
		set Tag [expr $i+1];
 
		HeatFluxBC -HTEntity $Tag -face 1  -type -Prescribed -HTConstants 2;
		HeatFluxBC -HTEntity $Tag -face 1  -type -ConvecAndRad -HTConstants 2;
	}

}

#Recorder

HTNodeSet 10 -Locx 0.0 -Locy 0.0;


HTRecorder -file tempProbe.out -NodeSet 10;    #Nodes of beam on the beam-slab interface 

#HTRecorder -file temp1.out -NodeSet 2;    #Nodes of slab on the beam-slab interface 
#HTRecorder -file temp2.out -NodeSet 5;    #All the nodes along x=0.0

HTAnalysis HeatTransfer
HTAnalyze 30 10;	
#HTAnalyze 80 15;
wipeHT;


