
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
set ny 10;
set slabT 0.1;
set slabB 3.150;
set slabL 4.150;
set UDL 5.4E3;

nDMaterial ElasticIsotropic3DThermal 12 2.2e10 0.2 0 1.4e-5;
#nDMaterial DruckerPragerThermal 2 $k $G $sigY $rho $rhoBar $Kinf $K0 $delta1 $delta2 $H $theta $mDen;
#nDMaterial PlateFiberThermal 4 12;
# $secTag $matTag $thickness
#section PlateFiberThermal 2 4 0.05;

#nDMaterial Damage2p $matTag $fcc 
#nDMaterial Damage2p 2 30e6 3e6; 
#nDMaterial ElasticIsotropic3DThermal 2 2e11 0.3 0 1.2e-5 ;
#nDMaterial Damage2p 2 30e6 3e6; 
#nDMaterial ElasticIsotropic3DThermal 2 2e11 0.3 0 1.2e-5 ;

#set k 2.13e10;set G 1.6e10;set sigY 6.05e6;set rho 0.174;set rhoBar 0.174;set Kinf 0;set K0 0;set delta1 1;set H 1.75e10;set theta 1;set delta2 0;set mDen 2400;
#nDMaterial DruckerPragerThermal 2 $k $G $sigY $rho $rhoBar $Kinf $K0 $delta1 $delta2 $H $theta $mDen $SigT;
#          CapPlasticity $tag $ndm  $rho     $G    $K  $X  $D  $W   $R  $lambda  $theta  $beta   $alpha    $T       $tol
#nDMaterial CapPlasticityThermal 2    3    2400  $G  $k  0.5032e8 4.6412e-10 0.42 4.43 7.9979e6 0.11 6.3816e-8 2.6614e7 -2.0684e6 1.0e-8
#nDMaterial J2Plasticity 2 $k $G $sigY $sigY 0 0
#nDMaterial Damage2p 2 30;
puts "3dhere0";
#nDMaterial PlateFiberThermal 4 12;
set gt [expr 4e6/2.2e10*4e6*2];
set gc [expr 37.0e6/2.2e10*37.0e6*6];
nDMaterial  CDPPlaneStressThermal 100 2.2e10 0.2 4e6 37e6 $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    4   100   10e9;

set gc1 [expr 37.0e6/2.2e10*37.0e6*100];
nDMaterial  CDPPlaneStressThermal 101 2.2e10 0.2 4e6 80e6 $gt $gc1;
nDMaterial   PlateFromPlaneStressThermal    34   101   10e9;


#nDMaterial DruckerPragerSteelThermal 23 $Bar_k $Bar_G $Bar_sigY $Bar_rho $Bar_rhoBar $Bar_Kinf $Bar_K0 $Bar_delta1 $Bar_delta2 $Bar_H $Bar_theta $Bar_mDen;
nDMaterial ElasticIsotropic3DThermal 23 2.2e10 0.2 0 1.2e-5 -cSoft ;
nDMaterial PlateFiberThermal 24 23;

uniaxialMaterial SteelECThermal 1 EC2NH 4.68e8 2e11;
#uniaxialMaterial ElasticThermal 1 20e11 1.2e-5;
nDMaterial PlateRebarThermal 3 1 0;
nDMaterial PlateRebarThermal 5 1 90;

nDMaterial  J2PlaneStressThermal 10 2e11 0.3 4.68e8 5.45e8 0 2e11;
nDMaterial   PlateFromPlaneStressThermal    44   10   20e10;

# $secTag $matTag $thickness
#section PlateFiberThermal 2 4 $slabT;

#D147Mesh
section LayeredShellThermal   1  13 4  0.01 4  0.01  4 0.005 3 0.000198 5 0.000198 4 0.004604  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
section LayeredShellThermal   2  12 4  0.01 4  0.01  4 0.0047175 44 0.000565 4 0.0047175  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
#661Mesh
section LayeredShellThermal   3  13 4  0.01 4  0.01  4 0.005 3 0.000198 5 0.000198 4 0.004604  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
section LayeredShellThermal   4  12 4  0.01 4  0.01  4 0.0047175 44 0.000565 4 0.0047175  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
#HD12 bars
section LayeredShellThermal   5  13 4  0.01 4  0.01  4 0.005 3 0.000565 5 0.000565 4 0.00387  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
section LayeredShellThermal   6  12 4  0.01 4  0.01  4 0.0047175 44 0.000565 4 0.0047175  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;

puts "here0";

#section LayeredShellThermal  2  10  4  0.01 4 0.01  4 0.01  4  0.01  4 0.01 4  0.01  4  0.01 4  0.01 4  0.01 4  0.01 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
block2D $nx $ny 1 1 ShellNLDKGQThermal  6 { 
    1   0. 0. 0.
    2   4.15 0. 0.
    3  4.15 3.15 0.
    4   0. 3.15 0.
}

#fully simply supported
fixX 0  0 0 1 0 0 0 ;
fixX 4.15  0 0 1 0 0 0 ;
fixY 0  0 0 1 0 0 0 ;
fixY 3.15  0 0 1 0 0 0 ;
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
recorder Node -file ShellData/DFreeSlabxT.out -time -nodeRange [expr 1+($nx+1)*$ny/2] [expr ($ny/2+1)*($nx+1)] -dof 1  disp;
recorder Node -file ShellData/DFreeSlabyT.out -time -nodeRange [expr 1+($nx+1)*$ny/2] [expr ($ny/2+1)*($nx+1)] -dof 2 disp;
recorder Node -file ShellData/DFreeSlabzT.out -time -nodeRange [expr 1+($nx+1)*$ny/2] [expr ($ny/2+1)*($nx+1)] -dof 3 disp;
recorder Node -file ShellData/DFreeSlabzTCent.out -time -node [expr ($ny/2)*($nx+1)+$nx/2] -dof 3 disp;
#------------------------------------
recorder Element -file ShellData/dtdcElerange110_12.out -time -eleRange 1 10  material 1 fiber 12 TempAndElong;
recorder Element -file ShellData/dtdcElerange110_9.out -time -eleRange 1 10  material 1 fiber 9 TempAndElong;

recorder Element -file ShellData/dtdcElerange2.out -time -eleRange 51 60  material 1 fiber 12 TempAndElong;

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

recorder Element -file ShellData/EleForceSec15sigma.out -time -ele $MidEle  material 1 fiber 15 stress;	
recorder Element -file ShellData/EleForceSec15Eps.out -time -ele $MidEle  material 1 fiber 15 strain;
recorder Element -file ShellData/EleForceSec15Temp.out -time -ele $MidEle   material 1 fiber 15 TempAndElong;

#-----------------------------------------------------------

recorder Element -file ShellData/CornerEleForceSec1sigma.out -time -ele $CorEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/CornerEleForceSec1Eps.out -time -ele $CorEle  material 1 fiber 1 strain;
recorder Element -file ShellData/CornerEleForceSec1Temp.out -time -ele $CorEle   material 1 fiber 1 TempAndElong;

recorder Element -file ShellData/CornerEleForceSec12sigma.out -time -ele $CorEle  material 1 fiber 12 stress;	
recorder Element -file ShellData/CornerEleForceSec12Eps.out -time -ele $CorEle  material 1 fiber 12 strain;
recorder Element -file ShellData/CornerEleForceSec12Temp.out -time -ele $CorEle   material 1 fiber 12 TempAndElong;

recorder Element -file ShellData/CornerEleForceSec11sigma.out -time -ele $CorEle  material 1 fiber 11 stress;	
recorder Element -file ShellData/CornerEleForceSec11Eps.out -time -ele $CorEle  material 1 fiber 11 strain;
recorder Element -file ShellData/CornerEleForceSec11Temp.out -time -ele $CorEle   material 1 fiber 11 TempAndElong;

recorder Element -file ShellData/CornerEleForceSec4sigma.out -time -ele $CorEle  material 1 fiber 4 stress;	
recorder Element -file ShellData/CornerEleForceSec4Eps.out -time -ele $CorEle  material 1 fiber 4 strain;
recorder Element -file ShellData/CornerEleForceSec4Temp.out -time -ele $CorEle   material 1 fiber 4 TempAndElong;

#------------------------------------

recorder Element -file ShellData/SideEleForceSec1sigma.out -time -ele $SideEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/SideEleForceSec1Eps.out -time -ele $SideEle  material 1 fiber 1 strain;
recorder Element -file ShellData/SideEleForceSec1Temp.out -time -ele $SideEle   material 1 fiber 1 TempAndElong;

recorder Element -file ShellData/SideEleForceSec12sigma.out -time -ele $SideEle  material 1 fiber 12 stress;	
recorder Element -file ShellData/SideEleForceSec12Eps.out -time -ele $SideEle  material 1 fiber 12 strain;
recorder Element -file ShellData/SideEleForceSec12Temp.out -time -ele $SideEle   material 1 fiber 12 TempAndElong;

recorder Element -file ShellData/SideEleForceSec11sigma.out -time -ele $SideEle  material 1 fiber 11 stress;	
recorder Element -file ShellData/SideEleForceSec1Eps.out -time -ele $SideEle  material 1 fiber 11 strain;
recorder Element -file ShellData/SideEleForceSec11Temp.out -time -ele $SideEle   material 1 fiber 11 TempAndElong;

recorder Element -file ShellData/SideEleForceSec4sigma.out -time -ele $SideEle  material 1 fiber 4 stress;	
recorder Element -file ShellData/SideEleForceSec4Eps.out -time -ele $SideEle  material 1 fiber 4 strain;
recorder Element -file ShellData/SideEleForceSec4Temp.out -time -ele $SideEle   material 1 fiber 4 TempAndElong;

#------------------------------------
puts "here1";


if {$ANALYSIS == "HasPoint"} {
set UDLP [expr $UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
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
integrator LoadControl 0.1;	
analysis Static;			
analyze 10;
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
	
	#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc $StartNodeTag 0 $MidNodeTag 0.5 $EndNodeTag 1; 
	eleLoad -range 1 $NumEles -type -shellThermal -source "slab.dat" [expr -$slabT/2] [expr $slabT/2];
	#eleLoad -range 1 $NumEles -type -shellThermal 878.293 -0.05 620.924 -0.0375 449.884 -0.025 331.591 -0.0125 246.801 0 184.89 0.0125 141.852 0.025 116.437 0.0375 104.777 0.5;
	#eleLoad -range 1 $NumEles -type -shellThermal 1000 [expr -$slabT/2] 100 [expr $slabT/2];
	
#}
#}
#eleLoad -ele 1 -type -beamThermal 600 -$HalfD 600 $HalfD;
}

#wipe;
constraints Plain;
numberer Plain;
system BandGeneral;
#test NormUnbalance 1.0e-4 10 1;
test NormDispIncr 1e-3  500 1;
algorithm Newton;
integrator LoadControl 30;	
analysis Static;			
analyze 120;

}

wipe;
