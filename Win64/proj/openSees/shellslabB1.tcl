
# units: m,N
wipe;	
				
				
set ANALYSIS "HasPoint";
set TANALYSIS "Has0Thermo"; 

model BasicBuilder -ndm 3 -ndf 6;

file mkdir ShellDataB1;

source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl

#8.72
#Concrete model
#these should both be even, number of elements per edge
set nx 20;
set ny 14;
set slabT 0.0682;
set slabB 1.829;
set slabL 2.745;
set UDL 1.0E3;


puts "here0";
#nDMaterial PlateFiberThermal 4 7;
set gt [expr 1.87e6/1.12e10*1.87e6*2];
set gc [expr 18.7e6/1.12e10*18.7e6*6];
nDMaterial  CDPPlaneStressThermal 7  1.12e10 0.2  1.87e6  18.7e6 $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    4   7    1e9


uniaxialMaterial SteelECThermal 1 EC2NH 4.50e8 2e11;
#uniaxialMaterial ElasticThermal 1 2e11 1.2e-5;
nDMaterial PlateRebarThermal 3 1 0;
nDMaterial PlateRebarThermal 5 1 90;

#nDMaterial  J2PlasticityThermal 10 22 2.0e11 0.3 4.5e8 5.45e8 0 0
#nDMaterial   PlateFiberThermal    44   10;

nDMaterial  J2PlaneStressThermal 10 21 2.06e11 0.3 4.5e8 5.45e8 0.1 0;
nDMaterial   PlateFromPlaneStressThermal    44   10   20e10;

# $secTag $matTag $thickness
#section PlateFiberThermal 2 4 $slabT;
#FOR NONLINEAR
#section LayeredShellThermal  2  14  4  0.01  4 0.009497345 3 0.000502655 5 0.000502655 4 0.009497345 4  0.01  4  0.01 4  0.01  4  0.01 4 0.009497345 5 0.000502655 3 0.000502655 4 0.009497345 4  0.01 ;
#section LayeredShellThermal  2  10  4  0.01  4 0.01 4 0.01 4  0.01  4  0.01 4  0.01  4  0.01 4  0.01  4 0.01 4 0.01 ;
#section LayeredShellThermal  2  11  4 0.008225  4 0.008095 24 0.00026 4 0.00685 4 0.00698  4 0.00698 4 0.00698 4 0.00685 24 0.00026 4 0.008095 4 0.008225 ;
#due to the section defition in ShellNLDKGQ, the layer order is reversed
section LayeredShellThermal  2  16  4 0.005457 4 0.005457 4 0.005197 3 0.00026 5 0.00026 4 0.00519 4 0.00545  4 0.00545 4 0.00545 4 0.00545 4 0.00519 5 0.00026 3 0.00026 4 0.006117 4 0.006377 4 0.006377;
section LayeredShellThermal  3 12  4 0.008185  4 0.008055 44 0.00026  4 0.00532 4 0.00545  4 0.00545 4 0.00545 4 0.00545 4 0.00532 44 0.00026 4 0.009435 4 0.009565 ;
puts "here1";

#section LayeredShellThermal  2  10  4  0.01 4 0.01  4 0.01  4  0.01  4 0.01 4  0.01  4  0.01 4  0.01 4  0.01 4  0.01 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
block2D $nx $ny 1 1 ShellNLDKGQThermal  3 { 
    1   0. 0. 0.
    2   2.745 0. 0.
    3  2.745 1.829 0.
    4   0. 1.829 0.
}

#fully simply supported
fixX 0  0 0 1 0 0 0 ;
fixX 2.745  0 0 1 0 0 0 ;
fixY 0 0 0 1 0 0 0 ;
fixY 1.829  0 0 1 0 0 0 ;
fix 1  1 1 1 0 0 1 ;



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

puts "here";

set midSlabsb [expr $nx/2+($nx+1)*$ny/2];
#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

set CorEle 1;
set MidEle [expr 1+($nx)*$ny/2+$nx/2];
set SideEle [expr 1+($nx)*$ny/2];


recorder Node -file ShellDataB1/DFreeCentrezT.out -time -node [expr $nx/2+($nx+1)*$ny/2] -dof 3 disp;

#------------------------------------

recorder Element -file ShellDataB1/EleForceSec1sigma.out -time -ele $MidEle  material 1 fiber 1 stress;	
recorder Element -file ShellDataB1/EleForceSec1Eps.out -time -ele $MidEle  material 1 fiber 1 strain;
recorder Element -file ShellDataB1/EleForceSec1Temp.out -time -ele $MidEle   material 1 fiber 1 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec2sigma.out -time -ele $MidEle  material 1 fiber 2 stress;	
recorder Element -file ShellDataB1/EleForceSec2Eps.out -time -ele $MidEle  material 1 fiber 2 strain;
recorder Element -file ShellDataB1/EleForceSec2Temp.out -time -ele $MidEle   material 1 fiber 2 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec3sigma.out -time -ele $MidEle  material 1 fiber 3 stress;	
recorder Element -file ShellDataB1/EleForceSec3Eps.out -time -ele $MidEle  material 1 fiber 3 strain;
recorder Element -file ShellDataB1/EleForceSec3Temp.out -time -ele $MidEle   material 1 fiber 3 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec4sigma.out -time -ele $MidEle  material 1 fiber 4 stress;	
recorder Element -file ShellDataB1/EleForceSec4Eps.out -time -ele $MidEle  material 1 fiber 4 strain;
recorder Element -file ShellDataB1/EleForceSec4Temp.out -time -ele $MidEle   material 1 fiber 4 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec5sigma.out -time -ele $MidEle  material 1 fiber 5 stress;	
recorder Element -file ShellDataB1/EleForceSec5Eps.out -time -ele $MidEle  material 1 fiber 5 strain;
recorder Element -file ShellDataB1/EleForceSec5Temp.out -time -ele $MidEle   material 1 fiber 5 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec6sigma.out -time -ele $MidEle  material 1 fiber 6 stress;	
recorder Element -file ShellDataB1/EleForceSec6Eps.out -time -ele $MidEle  material 1 fiber 6 strain;
recorder Element -file ShellDataB1/EleForceSec6Temp.out -time -ele $MidEle   material 1 fiber 6 TempAndElong;


recorder Element -file ShellDataB1/EleForceSec7sigma.out -time -ele $MidEle  material 1 fiber 7 stress;	
recorder Element -file ShellDataB1/EleForceSec7Eps.out -time -ele $MidEle  material 1 fiber 7 strain;
recorder Element -file ShellDataB1/EleForceSec7Temp.out -time -ele $MidEle   material 1 fiber 7 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec8sigma.out -time -ele $MidEle  material 1 fiber 8 stress;	
recorder Element -file ShellDataB1/EleForceSec8Eps.out -time -ele $MidEle  material 1 fiber 8 strain;
recorder Element -file ShellDataB1/EleForceSec8Temp.out -time -ele $MidEle   material 1 fiber 8 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec9sigma.out -time -ele $MidEle  material 1 fiber 9 stress;	
recorder Element -file ShellDataB1/EleForceSec9Eps.out -time -ele $MidEle  material 1 fiber 9 strain;
recorder Element -file ShellDataB1/EleForceSec9Temp.out -time -ele $MidEle   material 1 fiber 9 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec10sigma.out -time -ele $MidEle  material 1 fiber 10 stress;	
recorder Element -file ShellDataB1/EleForceSec10Eps.out -time -ele $MidEle  material 1 fiber 10 strain;
recorder Element -file ShellDataB1/EleForceSec10Temp.out -time -ele $MidEle   material 1 fiber 10 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec11sigma.out -time -ele $MidEle  material 1 fiber 11 stress;	
recorder Element -file ShellDataB1/EleForceSec11Eps.out -time -ele $MidEle  material 1 fiber 11 strain;
recorder Element -file ShellDataB1/EleForceSec11Temp.out -time -ele $MidEle   material 1 fiber 11 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec12sigma.out -time -ele $MidEle  material 1 fiber 12 stress;	
recorder Element -file ShellDataB1/EleForceSec12Eps.out -time -ele $MidEle  material 1 fiber 12 strain;
recorder Element -file ShellDataB1/EleForceSec12Temp.out -time -ele $MidEle   material 1 fiber 12 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec13sigma.out -time -ele $MidEle  material 1 fiber 13 stress;	
recorder Element -file ShellDataB1/EleForceSec13Eps.out -time -ele $MidEle  material 1 fiber 13 strain;
recorder Element -file ShellDataB1/EleForceSec13Temp.out -time -ele $MidEle   material 1 fiber 13 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec14sigma.out -time -ele $MidEle  material 1 fiber 14 stress;	
recorder Element -file ShellDataB1/EleForceSec14Eps.out -time -ele $MidEle  material 1 fiber 14 strain;
recorder Element -file ShellDataB1/EleForceSec14Temp.out -time -ele $MidEle   material 1 fiber 14 TempAndElong;

recorder Element -file ShellDataB1/EleForceSec15sigma.out -time -ele $MidEle  material 1 fiber 15 stress;	
recorder Element -file ShellDataB1/EleForceSec15Eps.out -time -ele $MidEle  material 1 fiber 15 strain;
recorder Element -file ShellDataB1/EleForceSec15Temp.out -time -ele $MidEle   material 1 fiber 15 TempAndElong;
puts "here1";
#print domain.out


if {$ANALYSIS == "HasPoint"} {
set UDLP [expr -$UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
    # for {set numx 1} {$numx <= $nx} {incr numx} {
	#for {set numy 1} {$numy <= $ny} {incr numy} {  
	#	set nodeID [expr $numx+ ($nx+1)*($numy-1)];
	#	load $nodeID 0 0 $UDLP 0 0 0 ;
	#	}
	#	}
 set NumNodes [expr ($nx+1)*($ny+1)]
	for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		load $nodeID 0 0 $UDLP 0 0 0 ;
		}
    #set Load 100;	
	#for {set ID 0} {$ID<=$ny} {incr ID} {
		#set nodeID [expr ($nx+1)*($ID+1)];
		#load $nodeID  0 0 [expr $Load/($ny+1)] 0 0 0 ;
	#}
		
}
puts "Point";



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-3  300 1;
algorithm Newton;
integrator LoadControl 0.5;	
analysis Static;			
analyze 200;
loadConst -time 0.0
}



if {$TANALYSIS == "HasThermo"} {

puts "Thermal action to slab"

set NumEles [expr $nx*$ny];

set StartNodeTag [expr ($nx+1)*($ny/2)+$nx/2+1]
set MidNodeTag [expr ($nx+1)*($ny*3/4)+1+$nx*3/4]
set EndNodeTag [expr ($nx+1)*($ny+1)]

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
	load $StartNodeTag -nodalThermal 800 $minusHalfD 200 $HalfD;
	load $MidNodeTag -nodalThermal 200 $minusHalfD 50 $HalfD
	load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD
	
	eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc $StartNodeTag 0 $MidNodeTag 0.5 $EndNodeTag 1; 
	#eleLoad -range 1 $NumEles -type -shellThermal 400 [expr -$slabT/2] 100 [expr $slabT/2];
	
#}
#}
#eleLoad -ele 1 -type -beamThermal 600 -$HalfD 600 $HalfD;
}

#wipe;
constraints Plain;
numberer Plain;
system BandGeneral;
#test NormUnbalance 1.0e-4 10 1;
test NormDispIncr 1e-4  300 1;
algorithm Newton;
integrator LoadControl 0.005;	
analysis Static;			
analyze 200;

}

wipe;
