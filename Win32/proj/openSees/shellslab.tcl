#This is written to test the new development for thermo-mechanical shell analysis;
#The development contains the following parts--
#nDMaterial: ElasticIsotropicThermal, J2PlasticityThermal, DruckerPragerThermal, CDPPlanestressThermal
#ndMaterialWrapper: PlateFiberThermal, Platefrom
# units: m,N
wipe;	
				
				
set ANALYSIS "HasLoad";
set TANALYSIS "HasThermo"; 

model BasicBuilder -ndm 3 -ndf 6;

source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl
set nx 6;
set ny 6;
#2.0m*1.0m*0.1m slab;
set UDL 1E4;
set slabT 0.05;

#Thermo-mechanical NDMaterial 1:ElasticIsotropicThermal(for concrete or steel, with options for stiffness softening)
nDMaterial ElasticIsotropicThermal 1 2e10 0.2 0 1.2e-5 -cSoft;

#Thermo-mechanical NDMaterial 2:DruckerPragerThermal(for concrete)
#---bulk modulus
set k 6233.333333e6;
#---shear modulus
set G 4675e6;
#---yield stress
set sigY 5.170506912e6;
#---failure surface and associativity
set rho 0.590737158;set rhoBar 0;
#---isotropic hardening
set Kinf 0;set K0 0;set delta1 0;
#---kinematic hardening
set H 6.23e9;set theta 1;
#---tension softening
set delta2 1;
set mDen 2400;
set sigT 6e6;
nDMaterial DruckerPragerThermal 2 $k $G $sigY $rho $rhoBar $Kinf $K0 $delta1 $delta2 $H $theta $mDen;

#Thermo-mechanical NDMaterial 3:J2PlasticityThermal(for steel)
nDMaterial J2PlasticityThermal 3 $k $G $sigY $sigY 0 0;

nDMaterial CDPPlaneStressThermal 4 3e10 0.2 3e6 30e6;

#Thermo-mechanical NDMaterial 4:PlateFiberThermal(plate fiber from 3d material )
nDMaterial PlateFiberThermal 10 1;
nDMaterial PlateFromPlaneStressThermal 11 4 1e10;

#Thermo-mechanical Platefiber section(membrane plate section ) # $secTag $matTag $thickness
section PlateFiberThermal 2 10 $slabT;

uniaxialMaterial SteelECThermal 21 EC2NH 2.8e8 2e11;
nDMaterial PlateRebarThermal 23 21 0;
nDMaterial PlateRebarThermal 24 21 90;

#Thermo-mechanical layered shell section( steel mesh to represent reinforcement)
section LayeredShellThermal  3  5  11 0.01 11 0.01 11 0.01 11 0.01 11 0.01 ;
puts "here0";
# ShellMITC4Thermal:Geo-nonlinear, ShellNLDKGQThermal: Geo and mat-nonlinear
block2D $nx $ny 1 1 ShellMITC4Thermal 3 { 
    1   0. 0. 0.
    2   2.0 0. 0.
    3  2.0 1.0 0.
    4   0. 1.0 0.
}

#fully simply supported
fixX 0  1 1 1 0 0 0 ;
fixX 2  0 0 1 0 0 0 ;



#fixY 0.2   0 1 0 1 0 1 ;

#fix 6   1 1 1 1 1 1 ;
#fix 11   1 0 1 1 1 1 ;
#fixX 3.0 0 0 1 0 0 0 
#fix 1 1 1 1 0 0 0;
#fix [expr 1+($nx+1)*$ny]  0 0 1 0 0 0 ;
#fix [expr ($nx+1)] 0 0 1 0 0 0 ;
#fix [expr ($nx+1)*($ny+1)] 0 0 1 0 0 0;
#fixY 0. 0 1 0 0 0 0;
#fixY 1 0 1 0 0 0 0;

recorder Node -file N_EndNode_DOF.out -time -node 24 -dof 1 2 3 disp;

puts "here";

#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 5;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 5;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

set CorEle 1;
set MidEle [expr 1+($nx)*$ny/2+$nx/2];
set SideEle [expr 1+($nx)*$ny/2];

puts "here1";
#print domain.out


if {$ANALYSIS == "HasLoad"} {
set UDLP [expr $UDL*2.0*1.0/$nx/$ny];

pattern Plain 1 Linear {
   set NumNodes [expr ($nx+1)*($ny+1)]
	for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		load $nodeID 0 0 $UDLP 0 0 0 ;
		}
   # set Load 100;	
	#for {set ID 0} {$ID<=$ny} {incr ID} {
	#	set nodeID [expr ($nx+1)*($ID+1)];
	#	load $nodeID  0 0 [expr $Load/($ny+1)] 0 0 0 ;
	#}
		
}
puts "Point";



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-4  100 1;
algorithm Newton;
integrator LoadControl 0.01;	
analysis Static;			
analyze 100;
loadConst -time 0.0
}



if {$TANALYSIS == "HasThermo"} {

puts "Thermal action to slab"

set NumEles [expr $nx*$ny];

#set StartNodeTag [expr ($nx+1)*($ny/2)+$nx/2+1]
#set MidNodeTag [expr ($nx+1)*($ny*13/20)+1+$nx*13/20]
#set SecMidNodeTag [expr ($nx+1)*($ny*17/20)+1+$nx*17/20]
#set EndNodeTag [expr ($nx+1)*($ny+1)]

set StartNodeTag 1
set MidNodeTag [expr ($nx+1)*($ny*7/20)+1+$nx*7/20]
set SecMidNodeTag [expr ($nx+1)*($ny*14/20)+1+$nx*14/20]
set EndNodeTag [expr ($nx+1)*($ny+1)]

#set StartNodeTag 1
#set MidNodeTag [expr ($nx+1)*($ny/2)+1+$nx/2]
#set EndNodeTag [expr ($nx+1)*($ny+1)]

set minusHalfD [expr -$slabT/2];
set HalfD [expr $slabT/2];



pattern Plain 3 Linear {
#set shellEndID [expr 100+$nx*$ny]
#for {set m 1} {$m<=$ny} {incr m} {
 # for {set n 1} {$n<=$nx} {incr n} {
   # set shellLoadID [expr 10000+($n-1)*100+$m-1];
	#set shellID [expr 1+($m-1)*$nx+$n-1];
    #set filepath "SlabTA/HTNodeRecorder";
    #set fileType ".dat";
    #append filepath $shellLoadID $fileType;
    #eleLoad -ele $shellID -type -shellThermal -source $filepath -50 50;
	#eleLoad -range 1 $NumEles -type -shellThermal 800 [expr -$slabT/2] 400 [expr $slabT/2];
	
#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc $StartNodeTag 0 $MidNodeTag 0.3 $MidNodeTag 0.7 $EndNodeTag 1; 
	eleLoad -range 1 $NumEles -type -shellThermal 400 [expr -$slabT/2] 100 [expr $slabT/2];
	
#}
#}
#eleLoad -ele 1 -type -beamThermal 600 -$HalfD 600 $HalfD;
}



#wipe;
constraints Plain;
numberer Plain;
system BandGeneral;
#test NormUnbalance 1.0e-4 10 1;
test NormDispIncr 1e-4  100 1;
algorithm Newton;
integrator LoadControl 0.05;	
analysis Static;			
analyze 20;

}
print domain.out;

