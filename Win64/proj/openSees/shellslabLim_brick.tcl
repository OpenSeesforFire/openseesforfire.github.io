# Three flat slab in Lim Test can be analysed. (D147,661 and HD 12)
# units: m,N

wipe;					
				
set ANALYSIS "HasPoint"; # ambient loads
set TANALYSIS "StandardThermo"; 

model BasicBuilder -ndm 3 -ndf 3;

file mkdir BrickData;

# source DisplayPlane.tcl
# source DisplayModel2D.tcl
# source DisplayModel3D.tcl


set nx 20;
set ny 20;
set nz 10 
set slabT 0.1;
set slabL 4.0;
set slabB 3.0;
set UDL 5.4E3;
set fyD147 5.65e8;
set typeD147 EC2NC;

set RebarMesh 1
set SectionType $RebarMesh

set Specimen D147

#nDMaterial PlateFiberThermal 4 12;
set gt [expr 3.7e6/2.2e10*3.7e6*2];
set gc [expr 37.0e6/2.2e10*37.0e6*6];

# concrete
# nDMaterial  CDPPlaneStressThermal 100 2.2e10 0.2 3.7e6 37e6 $gt $gc;
# nDMaterial   PlateFromPlaneStressThermal    4   100   1e9;

nDMaterial ElasticIsotropic3DThermal  12 2e10 0.2 0 1.4e-5 -CSoft;
nDMaterial PlateFiberThermal 4 12;
# if {$Specimen == "HD12"} {
	# uniaxialMaterial SteelECThermal 1 $typeHD12 $fyHD12  2e11;
# } elseif {$Specimen == "D147"} {
	# uniaxialMaterial SteelECThermal 1 $typeD147 $fyD147 2e11;
# } elseif {$Specimen == "661"} {
	# uniaxialMaterial SteelECThermal 1 $type661 $fy661 2e11;
# } 
  
#nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)  
 # nDMaterial PlateRebarThermal  3        1                0;
 # nDMaterial PlateRebarThermal  5        1                90;

#nDMaterial  J2PlasticityThermal 23 $K $G 3.45e8 4.00e8 0.1 0;
#nDMaterial   PlateFiberThermal    44   23;

#D147 reinfircement J2
# nDMaterial  J2PlaneStressThermal 11 22 2e11 0.3 5.65e8 5.85e8 0.1 2e8;
# nDMaterial   PlateFromPlaneStressThermal    44   11   20e10;

#661 reinfircement J2
# nDMaterial  J2PlaneStressThermal 12 22 2e11 0.3 5.68e8 6.45e8 0.1 2e8;
# nDMaterial   PlateFromPlaneStressThermal    45   12   20e10;

#HD12 reinfircement J2
# nDMaterial  J2PlaneStressThermal 13 21 2e11 0.3 4.68e8 6.00e8 0.1 2e8;
# nDMaterial   PlateFromPlaneStressThermal    46   13   20e10;
# if {$Specimen == "D147"} {
# #D147Mesh
 # #Rebar mesh
	# #section LayeredShellThermal   $RebarMesh  13 4  0.005 4  0.01  4 0.009802 3 0.000198 5 0.000198 4 0.004802  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
    # section LayeredShellThermal   $RebarMesh   4  0.01  4  0.01  4  0.01  4  0.01  4  0.01  4  0.01  4  0.01  4  0.01  4  0.01  4  0.01;
 # #J2 steel layer
	# section LayeredShellThermal   $SteelLayer  12 4  0.005 4  0.01  4 0.01 44 0.000198 4 0.004802  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
# } elseif {$Specimen == "661"} {
	# #661Mesh
	# #Rebar mesh
	# section LayeredShellThermal   $RebarMesh  13 4  0.01 4  0.01  4 0.005 3 0.000295 5 0.000295 4 0.00441  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
	# #J2 steel layer
	# section LayeredShellThermal   $SteelLayer  12 4  0.01 4  0.01  4 0.005 45 0.000295 4 0.004705  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
# } elseif {$Specimen == "HD12"} {
	# #HD12 bars
	# #Rebar mesh
	# section LayeredShellThermal   $RebarMesh  13 4  0.01 4  0.01  4 0.005 3 0.000565 5 0.000565 4 0.00387  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
	# #J2 steel layer
	# section LayeredShellThermal   $SteelLayer  12 4  0.01 4  0.01  4 0.005 46 0.000565 4 0.004435 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
# }

#section LayeredShellThermal  2  10  4  0.01 4 0.01  4 0.01  4  0.01  4 0.01 4  0.01  4  0.01 4  0.01 4  0.01 4  0.01 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
# block2D $nx $ny 1 1 ShellNLDKGQThermal   $SectionType { 
	# 1   0. 0. 0.
    # 2   3.15 0. 0.
    # 3   3.15 3.15 0.
    # 4   0. 3.15 0.
# }
block3D $nx $ny $nz 1 1 BrickThermal  12 { 
	1      0.    0.   0.
    2      4.0   0.   0.
    3      4.0  3.0   0.
    4      0.   3.0   0.
	5      0.    0.  0.1
    6      4.0   0.  0.1
    7      4.0  3.0  0.1
    8      0.   3.0  0.1
}
	

#fully simply supported
# fixX 0  0 0 1 0 0 0 ;
# fixX $slabL 0 0 1 0 0 0 ;
# fixY 0  0 0 1 0 0 0 ;
# fixY $slabB  0 0 1 0 0 0 ;
# fix 1  1 1 0 0 0 1 ;

fixX 0  0 0 1;
fixX $slabL 0 0 1;
fixY 0  0 0 1 ;
fixY $slabB  0 0 1;
fix 1  1 1 0;

set MidNode [expr ($nx+1)*($ny/2)+$nx/2+1];
#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 1;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

set CorEle 1;
set MidEle [expr ($nx)*($ny/2)+$nx/2+1];
set SideEle [expr ($nx)*($ny/4)+$nx/2+1];

#recorder Node -file BrickData/Disp3Center.out -time -node $MidNode -dof 3 disp;

recorder Node -file BrickData/MIDU3.out -time -nodeRange [expr ($nx+1)*($ny/2)+1] [expr ($nx+1)*($ny/2)+$nx+1] -dof 3 disp;

# recorder Element -file BrickData/brickstrain1.out -time -ele $SideEle material 2 strain;
# recorder Element -file BrickData/brickstrain2.out -time -ele [expr $SideEle + 1*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain3.out -time -ele [expr $SideEle + 2*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain4.out -time -ele [expr $SideEle + 3*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain5.out -time -ele [expr $SideEle + 4*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain6.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain7.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain8.out -time -ele [expr $SideEle + 7*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain9.out -time -ele [expr $SideEle + 8*($nx*$ny)] material 2 strain;
# recorder Element -file BrickData/brickstrain10.out -time -ele [expr $SideEle + 9*($nx*$ny)] material 2 strain;

# recorder Element -file BrickData/brickstress1.out -time -ele $SideEle material 2 stress;
# recorder Element -file BrickData/brickstress2.out -time -ele [expr $SideEle + 1*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress3.out -time -ele [expr $SideEle + 2*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress4.out -time -ele [expr $SideEle + 3*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress5.out -time -ele [expr $SideEle + 4*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress6.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress7.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress8.out -time -ele [expr $SideEle + 7*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress9.out -time -ele [expr $SideEle + 8*($nx*$ny)] material 2 stress;
# recorder Element -file BrickData/brickstress10.out -time -ele [expr $SideEle + 9*($nx*$ny)] material 2 stress;
	
# recorder Element -file BrickData/brickTempElong1.out -time -ele $SideEle  material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong2.out -time -ele [expr $SideEle + 1*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong3.out -time -ele [expr $SideEle + 2*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong4.out -time -ele [expr $SideEle + 3*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong5.out -time -ele [expr $SideEle + 4*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong6.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong7.out -time -ele [expr $SideEle + 5*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong8.out -time -ele [expr $SideEle + 7*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong9.out -time -ele [expr $SideEle + 8*($nx*$ny)] material 2 TempAndElong;
# recorder Element -file BrickData/brickTempElong10.out -time -ele [expr $SideEle + 9*($nx*$ny)] material 2 TempAndElong;
						
puts "here1";

#######################################################################################################

if {$ANALYSIS == "HasPoint"} {
set UDLP [expr -$UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
   set NumNodes [expr ($nx+1)*($ny+1)]
	for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		load $nodeID 0 0 $UDLP;
		}	
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


if {$TANALYSIS == "StandardThermo"} {

set minusHalfD [expr -$slabT/2];
set HalfD [expr $slabT/2];
set NumEle [expr $nx*$ny*$nz]

pattern Plain 3 Linear {
	   #eleLoad -range firstNode lastNode 		-type -shellThermal -source "fileName" lowestLayer  $highestLayer
		eleLoad -range 1 	 $NumEle -type -brickThermal -source "slab.dat" 0  $slabT
}

set ok 0
set stepSize 30
set currentTime 0.0
set tFinal 7200

constraints Plain    		
numberer Plain			
system BandGeneral		
test NormUnbalance 100 2000 0;
#test NormDispIncr 1e-2 1000 1	
#algorithm KrylovNewton	
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static
analyze 240;

# while {$currentTime < $tFinal && $ok == 0} {
	# puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	# set ok [analyze 1];
	# if {$ok == 0} {
		# set currentTime [getTime]
		# puts "Analysis step successful. Current time is $currentTime"
	# } else {
		# puts "Analysis failed at time = $currentTime"
		
		# return
	# }
# }
# }
