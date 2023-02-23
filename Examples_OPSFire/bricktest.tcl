
# units: m,N
wipe;	
				
				
set ANALYSIS "HasPoint";
set TANALYSIS "HasThermo"; 

model BasicBuilder -ndm 3 -ndf 6;

file mkdir ShellData;

source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl

#8.72
#Concrete model
#these should both be even, number of elements per edge
set nx 10;
set ny 10
set nz 5;
set slabT 0.1;
set slabB 1.0;
set slabL 1.0;
set UDL 0.5E3;

#loaded nodes
set mid [expr (  ($nx+1)*($ny+1)+1 ) / 2 ];

#nDMaterial ElasticIsotropic3DThermal 2 1.92e10 0.2 0 1.4e-5;
#nDMaterial DruckerPragerThermal 2 $k $G $sigY $rho $rhoBar $Kinf $K0 $delta1 $delta2 $H $theta $mDen;
#nDMaterial PlateFiberThermal 4 2;

set gt [expr 3.0e6/1.79e10*3.0e6*2];
set gc [expr 30e6/1.79e10*30e6*6];
nDMaterial  CDPPlaneStressThermal 100  1.79e10 0.2  3.0e6  30e6 $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    4   100  10e9;


uniaxialMaterial SteelECThermal 1 EC2NH 3.45e8 2e11;
#uniaxialMaterial ElasticThermal 1 2e11 1.2e-5;
nDMaterial PlateRebarThermal 3 1 0;
nDMaterial PlateRebarThermal 5 1 90;
# $secTag $matTag $thickness
#section PlateFiberThermal 2 4 $slabT;
#FOR NONLINEAR
#section LayeredShellThermal  2  14  4  0.01  4 0.009497345 3 0.000502655 5 0.000502655 4 0.009497345 4  0.01  4  0.01 4  0.01  4  0.01 4 0.009497345 5 0.000502655 3 0.000502655 4 0.009497345 4  0.01 ;
#section LayeredShellThermal  2  10  4  0.01  4 0.01 4 0.01 4  0.01  4  0.01 4  0.01  4  0.01 4  0.01  4 0.01 4 0.01 ;

nDMaterial ElasticIsotropic3DThermal  12 2e10 0.0 0 1.0e-5;
#nDMaterial ElasticIsotropic3DThermal  12 1.79e10 0.2 0 1.4e-5 -CSoft;
nDMaterial PlateFiberThermal 14 12;
#nDMaterial ElasticIsotropic3DThermal  22 1.79e5 0.2 0 2.4e-5 -CSoft;
#nDMaterial PlateFiberThermal 24 22;

set Es 2e11; set v 0.3;
set K [expr $Es/2.0/(1-$v)];
set G [expr $Es/2.0/(1+$v)]; 
#nDMaterial  J2PlasticityThermal 23 22 $K $G 3.45e8 4.00e8 0.1 0;
#nDMaterial   PlateFiberThermal    24   23;

nDMaterial  J2PlaneStressThermal 23 21 2e11 0.3 3.45e8 4.00e8 0 0;
nDMaterial   PlateFromPlaneStressThermal    24   23   20e10;


#section LayeredShellThermal  2  14  4  0.01  4 0.009607301 3 0.000392699 5 0.000392699 4 0.009607301 4  0.01  4  0.01 4  0.01  4  0.01 4 0.009607301 3 0.000392699 5 0.000392699 4 0.009607301 4  0.01 ;
#section LayeredShellThermal  3  12  4  0.01  4 0.01  24 0.000392699 4 0.00960731 4  0.01  4  0.01 4  0.01  4  0.01  4 0.00960731  24 0.000392699 4 0.01 4 0.01;

section LayeredShellThermal  6  10  14  0.01  14 0.01 14 0.01 14  0.01  14  0.01 14  0.01  14  0.01 14  0.01  14 0.01 14 0.01 ;
section LayeredShellThermal  4  5  14  0.02 14 0.02  14 0.02  14 0.02 14 0.02;

section PlateFiberThermal 5 14 $slabT;
set offset 0.3;
section LayeredShellThermal  7  -offset $offset 10  14  0.01  14 0.01 14 0.01 14  0.01  14  0.01 14  0.01  14  0.01 14  0.01  14 0.01 14 0.01 ;


puts "here0";
#section LayeredShellThermal  2  10  4  0.005 4 0.005  4 0.005  4  0.005  4 0.005 4  0.005  4  0.005 4  0.005  4  0.005 4  0.005 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
brick3D $nx $ny $nz 1 1  ShellNLDKGQThermal 7 { 
    1   0. 0. 0.
    2   3 0. 0.
    3  3 0.4 0.
    4   0. 0.4 0.
}

#fully simply supported
fixX 0  0 0 1 0 0 0 ;
#fixX 3  0 0 1 1 0 1 ;
#fixX 6.0   0 1 0 1 0 1 ;
fixX 3.0   0 0 1 0 0 0 ;
#fixY 0  0 1 0 1 0 1 ;
#fixY 0.2   0 1 0 1 0 1 ;

#fix 6   1 1 1 1 1 1 ;
#fix 11   1 0 1 1 1 1 ;
#fixX 3.0 0 0 1 0 0 0 
fix 1 1 1 0 0 0 1;
#fix [expr 1+($nx+1)*$ny/2]  0 1 0 0 0 0 ;
#fix [expr ($nx+1)] 0 0 1 0 0 0 ;
#fix [expr ($nx+1)*($ny+1)] 0 0 1 0 0 0;
#fixY 0. 0 1 0 0 0 0;
#fixY 1 0 1 0 0 0 0;

set midSlabsb [expr $nx/2+($nx+1)*$ny/2];
#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.01;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

set MidEle [expr 1+($nx)*$ny/2+$nx/2];
set MidNode [expr ($nx+1)*$ny/2+$nx/2+1]
set EndEle [expr 1+$nx];
recorder Node -file ShellData/DFreeSlabxT.out -time -nodeRange [expr 1+($nx+1)*$ny/2] [expr ($ny/2+1)*($nx+1)] -dof 1  disp;
recorder Node -file ShellData/DFreeSlabyT.out -time -nodeRange [expr 1+($nx+1)*$ny/2] [expr ($ny/2+1)*($nx+1)] -dof 2 disp;
recorder Node -file ShellData/DFreeSlabzT.out -time -node $MidNode -dof 1 2 3 disp;

recorder Node -file ShellData/DFreeSlabzTUP.out -time -nodeRange 1 [expr ($nx+1)] -dof 3 disp;
recorder Node -file ShellData/DFreeSlabzTDo.out -time -nodeRange [expr 1+($nx+1)*$ny] [expr ($ny+1)*($nx+1)] -dof 3 disp;

recorder Element -file ShellData/EleForce.out -time -ele $MidEle force;	
recorder Element -file ShellData/EleMatForce.out -time -ele $MidEle material 2 force;	
recorder Element -file ShellData/EleMatDeform.out -time -ele $MidEle material 2 deformation;

recorder Element -file ShellData/EleForceSec1sigma.out -time -ele $MidEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/EleForceSec1Eps.out -time -ele $MidEle  material 1 fiber 1 strain;
recorder Element -file ShellData/EleForceSec1Temp.out -time -ele $MidEle   material 1 fiber 1 TempAndElong;

recorder Element -file ShellData/EleForceSec2sigma.out -time -ele $MidEle  material 1 fiber 2 stress;	
recorder Element -file ShellData/EleForceSec2Eps.out -time -ele $MidEle  material 1 fiber 2 strain;
recorder Element -file ShellData/EleForceSec2Temp.out -time -ele $MidEle   material 1 fiber 2 TempAndElong;

recorder Element -file ShellData/EleForceSec3sigma.out -time -ele $MidEle  material 1 fiber 3 stress;	
recorder Element -file ShellData/EleForceSec3Eps.out -time -ele $MidEle  material 1 fiber 3 strain;
recorder Element -file ShellData/EleForceSec3Temp.out -time -ele $MidEle   material 1 fiber 3 TempAndElong;

recorder Element -file ShellData/EleForceSec4sigma.out -time -ele $MidEle  material 1 fiber 4 stress;	
recorder Element -file ShellData/EleForceSec4Eps.out -time -ele $MidEle  material 1 fiber 4 strain;
recorder Element -file ShellData/EleForceSec4Temp.out -time -ele $MidEle   material 1 fiber 4 TempAndElong;

recorder Element -file ShellData/EleForceSec5sigma.out -time -ele $MidEle  material 1 fiber 5 stress;	
recorder Element -file ShellData/EleForceSec5Eps.out -time -ele $MidEle  material 1 fiber 5 strain;
recorder Element -file ShellData/EleForceSec5Temp.out -time -ele $MidEle   material 1 fiber 5 TempAndElong;

recorder Element -file ShellData/EleForceSec6sigma.out -time -ele $MidEle  material 1 fiber 6 stress;	
recorder Element -file ShellData/EleForceSec6Eps.out -time -ele $MidEle  material 1 fiber 6 strain;
recorder Element -file ShellData/EleForceSec6Temp.out -time -ele $MidEle   material 1 fiber 6 TempAndElong;


recorder Element -file ShellData/EleForceSec7sigma.out -time -ele $MidEle  material 1 fiber 7 stress;	
recorder Element -file ShellData/EleForceSec7Eps.out -time -ele $MidEle  material 1 fiber 7 strain;
recorder Element -file ShellData/EleForceSec7Temp.out -time -ele $MidEle   material 1 fiber 7 TempAndElong;

recorder Element -file ShellData/EleForceSec8sigma.out -time -ele $MidEle  material 1 fiber 8 stress;	
recorder Element -file ShellData/EleForceSec8Eps.out -time -ele $MidEle  material 1 fiber 8 strain;
recorder Element -file ShellData/EleForceSec8Temp.out -time -ele $MidEle   material 1 fiber 8 TempAndElong;

recorder Element -file ShellData/EleForceSec9sigma.out -time -ele $MidEle  material 1 fiber 9 stress;	
recorder Element -file ShellData/EleForceSec9Eps.out -time -ele $MidEle  material 1 fiber 9 strain;
recorder Element -file ShellData/EleForceSec9Temp.out -time -ele $MidEle   material 1 fiber 9 TempAndElong;

recorder Element -file ShellData/EleForceSec10sigma.out -time -ele $MidEle  material 1 fiber 10 stress;	
recorder Element -file ShellData/EleForceSec10Eps.out -time -ele $MidEle  material 1 fiber 10 strain;
recorder Element -file ShellData/EleForceSec10Temp.out -time -ele $MidEle   material 1 fiber 10 TempAndElong;

recorder Element -file ShellData/EleForceSec11sigma.out -time -ele $MidEle  material 1 fiber 11 stress;	
recorder Element -file ShellData/EleForceSec11Eps.out -time -ele $MidEle  material 1 fiber 11 strain;
recorder Element -file ShellData/EleForceSec11Temp.out -time -ele $MidEle   material 1 fiber 11 TempAndElong;

recorder Element -file ShellData/EleForceSec12sigma.out -time -ele $MidEle  material 1 fiber 12 stress;	
recorder Element -file ShellData/EleForceSec12Eps.out -time -ele $MidEle  material 1 fiber 12 strain;
recorder Element -file ShellData/EleForceSec12Temp.out -time -ele $MidEle   material 1 fiber 12 TempAndElong;

recorder Element -file ShellData/EleForceSec13sigma.out -time -ele $MidEle  material 1 fiber 13 stress;	
recorder Element -file ShellData/EleForceSec13Eps.out -time -ele $MidEle  material 1 fiber 13 strain;
recorder Element -file ShellData/EleForceSec13Temp.out -time -ele $MidEle   material 1 fiber 13 TempAndElong;

recorder Element -file ShellData/EleForceSec14sigma.out -time -ele $MidEle  material 1 fiber 14 stress;	
recorder Element -file ShellData/EleForceSec14Eps.out -time -ele $MidEle  material 1 fiber 14 strain;
recorder Element -file ShellData/EleForceSec14Temp.out -time -ele $MidEle   material 1 fiber 14 TempAndElong;

recorder Element -file ShellData/EleForceSec1Lefsigma.out -time -ele $EndEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/EleForceSec1LeftEps.out -time -ele $EndEle  material 1 fiber 1 strain;
recorder Element -file ShellData/EleForceSec1LeftTemp.out -time -ele $EndEle   material 1 fiber 1 TempAndElong;
#print domain.out


if {$ANALYSIS == "HasPoint"} {
#set UDLP [expr -$UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
   # set NumNodes [expr ($nx+1)*($ny+1)]
	#for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		#load $nodeID 0 0 $UDLP 0 0 0 ;
		#}
    set Load -2000;	
	for {set ID 0} {$ID<=$ny} {incr ID} {
		set nodeID [expr ($nx+1)*($ID+1)-$nx/2];
		load $nodeID  0 0 [expr $Load/($ny+1)] 0 0 0 ;
	}
		
}
puts "Point";



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-3 600 2;
algorithm Newton;
integrator LoadControl 0.2;	
analysis Static;			
analyze 5;
loadConst -time 0.0
}


if {$TANALYSIS == "HasThermo"} {

puts "Thermal action to slab"

set NumEles [expr $nx*$ny];
set MidNodeTag [expr ($nx/2+1)*($ny/2+1)]
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
	load 1 -nodalThermal 800 $minusHalfD 400 $HalfD;
	load $MidNodeTag -nodalThermal 200 $minusHalfD 100 $HalfD
	load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD
	
	#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $MidNodeTag 0.5 $EndNodeTag 1; 
	eleLoad -range 1 $NumEles -type -shellThermal 500 [expr -$HalfD-$offset] 100 [expr $HalfD-$offset]
	
#}
#}
#eleLoad -ele 1 -type -beamThermal 600 -$HalfD 100 $HalfD;
}

#wipe;
constraints Plain;
numberer Plain;
system BandGeneral;
test NormUnbalance 0.1 1000 1;
#test NormDispIncr 1e-3  600 1;
algorithm Newton;
integrator LoadControl 0.01;	
analysis Static;			
analyze 100;

}

wipe;
