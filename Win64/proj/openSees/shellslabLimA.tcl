# Three flat slab in Lim Test can be analysed. (D147,661 and HD 12)
# units: m,N

wipe;					
				
set ANALYSIS "Has0Point"; # ambient loads
#if TANALYSIS = StandardThermo, standard thermal load will be applied; if = LocalisedThermo, localised thermal load will be applied  

# set TANALYSIS "LocalisedThermo"; 
set TANALYSIS "StandardThermo"; 

model BasicBuilder -ndm 3 -ndf 6;

file mkdir ShellDataLim147;

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
set fyD147 5.65e8;
set fy661  5.68e8;
set fyHD12 4.68e8;
set typeD147 EC2NC;
set type661 EC2NC;
set typeHD12 EC2NH;

set RebarMesh 1
set SteelLayer 2
set SectionType $RebarMesh
# set SectionType $SteelLayer

set Specimen D147
# set Specimen HD12
# set Specimen "661"


#nDMaterial PlateFiberThermal 4 12;
set gt [expr 3.7e6/2.2e10*3.7e6*2];
set gc [expr 37.0e6/2.2e10*37.0e6*6];

# concrete
nDMaterial  CDPPlaneStressThermal 100 2.2e10 0.2 3.7e6 37e6 $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    4   100   1e9;




# Reinforcement Rebar mesh
# uniaxialMaterial SteelECThermal 1  EC2NC  $fy  $Es

if {$Specimen == "HD12"} {
	uniaxialMaterial SteelECThermal 1 $typeHD12 $fyHD12  2e11;
} elseif {$Specimen == "D147"} {
	uniaxialMaterial SteelECThermal 1 $typeD147 $fyD147 2e11;
} elseif {$Specimen == "661"} {
	uniaxialMaterial SteelECThermal 1 $type661 $fy661 2e11;
} 
  
#nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)  
 nDMaterial PlateRebarThermal  3        1                0;
 nDMaterial PlateRebarThermal  5        1                90;

#nDMaterial  J2PlasticityThermal 23 $K $G 3.45e8 4.00e8 0.1 0;
#nDMaterial   PlateFiberThermal    44   23;

#D147 reinfircement J2
puts "first difference here."
nDMaterial  J2PlaneStressThermal 11 22 2e11 0.3 5.65e8 5.85e8 0.1 2e8;
nDMaterial   PlateFromPlaneStressThermal    44   11   20e10;

#661 reinfircement J2
nDMaterial  J2PlaneStressThermal 12 22 2e11 0.3 5.68e8 6.45e8 0.1 2e8;
nDMaterial   PlateFromPlaneStressThermal    45   12   20e10;

#HD12 reinfircement J2
nDMaterial  J2PlaneStressThermal 13 21 2e11 0.3 4.68e8 6.00e8 0.1 2e8;
nDMaterial   PlateFromPlaneStressThermal    46   13   20e10;
if {$Specimen == "D147"} {
#D147Mesh
 #Rebar mesh
 puts "second difference here."
	section LayeredShellThermal   $RebarMesh  13 4  0.008 4  0.008  4 [expr 0.009 - 0.000198] 3 0.000198 5 0.000198 4 [expr 0.005 - 0.000198]  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01;
 #J2 steel layer
	section LayeredShellThermal   $SteelLayer  12 4  0.008 4  0.008  4 [expr 0.009 - 0.000198*0.5] 44 0.000198 4 [expr 0.005 - 0.000198*0.5]  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
	
} elseif {$Specimen == "HD12"} {
	#HD12 bars
	#Rebar mesh
	section LayeredShellThermal   $RebarMesh  13 4  0.008 4  0.008  4 [expr 0.009 - 0.000565]  3 0.000565 5 0.000565 4 [expr 0.005 - 0.000565] 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01;
	# J2 steel layer
	section LayeredShellThermal   $SteelLayer  12 4  0.008 4  0.008  4 [expr 0.009 - 0.000565*0.5] 46 0.000565 4 [expr 0.005 - 0.000565*0.5] 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
	
} elseif {$Specimen == "661"} {
	#661Mesh
	#Rebar mesh
	section LayeredShellThermal   $RebarMesh  13 4  0.008 4  0.008  4 [expr 0.009 - 0.000295]  3 0.000295 5 0.000295 4 [expr 0.005 - 0.000295] 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01;
	#J2 steel layer
	section LayeredShellThermal   $SteelLayer  12 4  0.008 4  0.008  4 [expr 0.009 - 0.000295*0.5] 45 0.000295 4 [expr 0.005 - 0.000295*0.5]  4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 4  0.01 ;
}

#section LayeredShellThermal  2  10  4  0.01 4 0.01  4 0.01  4  0.01  4 0.01 4  0.01  4  0.01 4  0.01 4  0.01 4  0.01 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
block2D $nx $ny 1 1 ShellNLDKGQThermal  $SectionType { 
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

recorder Node -file ShellDataLim147/Disp3Center.out -time -node [expr ($ny/2)*($nx+1)+$nx/2+1] -dof 3 disp;
recorder Node -file ShellDataLim147/U3.dat -time -nodeRange 1 [expr ($nx+ 1)*($ny + 1)] -dof 3 disp;
#------------------------------------
recorder	Element	-file	ShellDataLim147/BottomTempAndDt.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	1	TempAndElong;


if {$SectionType == $RebarMesh} {
	recorder	Element	-file	ShellDataLim147/LongSpanSteelStrain.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	4	strain;
	recorder	Element	-file	ShellDataLim147/ShortSpanSteelStrain.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	5	strain;
	recorder	Element	-file	ShellDataLim147/LongSpanSteelStress.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	4	stress;
	recorder	Element	-file	ShellDataLim147/ShortSpanSteelStress.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	5	stress;
	recorder	Element	-file	ShellDataLim147/TopTempAndDt.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	13	TempAndElong;
} elseif {$SectionType == $SteelLayer} {
	recorder	Element	-file	ShellDataLim147/SteelLayerStrain.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	4	strain;
	recorder	Element	-file	ShellDataLim147/SteelLayerStress.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	4	stress;
	recorder	Element	-file	ShellDataLim147/TopTempAndDt.out	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	12	TempAndElong;
}

recorder	Element	-file	ShellDataLim147/Pt1StressFiber1.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	1	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber2.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	2	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber3.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	3	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber4.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	4	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber5.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	5	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber6.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	6	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber7.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	7	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber8.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	8	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber9.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	9	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber10.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	10	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber11.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	11	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber12.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	12	stress;
recorder	Element	-file	ShellDataLim147/Pt1StressFiber13.dat	-time	-eleRange	1	[expr $nx*$ny]	material	1	fiber	13	stress;
												
recorder	Element	-file	ShellDataLim147/Pt2StressFiber1.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	1	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber2.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	2	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber3.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	3	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber4.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	4	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber5.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	5	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber6.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	6	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber7.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	7	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber8.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	8	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber9.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	9	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber10.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	10	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber11.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	11	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber12.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	12	stress;
recorder	Element	-file	ShellDataLim147/Pt2StressFiber13.dat	-time	-eleRange	1	[expr $nx*$ny]	material	2	fiber	13	stress;
												
recorder	Element	-file	ShellDataLim147/Pt3StressFiber1.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	1	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber2.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	2	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber3.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	3	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber4.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	4	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber5.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	5	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber6.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	6	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber7.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	7	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber8.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	8	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber9.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	9	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber10.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	10	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber11.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	11	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber12.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	12	stress;
recorder	Element	-file	ShellDataLim147/Pt3StressFiber13.dat	-time	-eleRange	1	[expr $nx*$ny]	material	3	fiber	13	stress;
												
recorder	Element	-file	ShellDataLim147/Pt4StressFiber1.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	1	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber2.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	2	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber3.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	3	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber4.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	4	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber5.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	5	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber6.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	6	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber7.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	7	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber8.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	8	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber9.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	9	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber10.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	10	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber11.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	11	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber12.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	12	stress;
recorder	Element	-file	ShellDataLim147/Pt4StressFiber13.dat	-time	-eleRange	1	[expr $nx*$ny]	material	4	fiber	13	stress;



#------------------------------------
puts "here1";

#######################################################################################################

if {$ANALYSIS == "HasPoint"} {
set UDLP [expr -$UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
   set NumNodes [expr ($nx+1)*($ny+1)]
	for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		load $nodeID 0 0 $UDLP 0 0 0 ;
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


pattern Plain 3 Linear {
	   #eleLoad -range firstNode lastNode 		-type -shellThermal -source "fileName" lowestLayer  $highestLayer
		eleLoad -range 1 		 [expr $nx*$ny] -type -shellThermal -source "slab.dat" $minusHalfD 		$HalfD
}

set ok 0
set stepSize 30
set currentTime 0.0
set tFinal 3600

constraints Plain    		
numberer Plain			
system BandGeneral		
test NormDispIncr 1e-3 1000 1	
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static
analyze 120
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
}


if {$TANALYSIS == "LocalisedThermo"} {

puts "Thermal action to slab"

set NumEles [expr $nx*$ny];

set StartNodeTag [expr ($nx+1)*($ny/2)+$nx/2+1]
set MidNodeTag [expr ($nx+1)*($ny*3/4)+1+$nx*3/4]
set EndNodeTag [expr ($nx+1)*($ny+1)]
puts $StartNodeTag","$MidNodeTag","$EndNodeTag;
set minusHalfD [expr -$slabT/2];
set HalfD [expr $slabT/2];

pattern Plain 3 Linear {
	# For localised distributed shell thermal action using appointed temperature
	load $StartNodeTag -nodalThermal -source "slab01.dat"  $minusHalfD $HalfD;
	load $MidNodeTag -nodalThermal -source "slab02.dat"  $minusHalfD $HalfD;
	load $EndNodeTag -nodalThermal -source "slab03.dat"  $minusHalfD $HalfD;
	eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc $StartNodeTag 0 $MidNodeTag 0.5 $EndNodeTag 1; 
	#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc $StartNodeTag 0 $EndNodeTag 1;
	

}
set ok 0
set stepSize 30
set currentTime 0.0
set tFinal 3600

constraints Plain    		
numberer Plain			
system BandGeneral		
test NormDispIncr 1e-3 1000 1
# test RelativeNormUnbalance 0.01 500	1
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static
analyze 100
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

}
# wipe;
