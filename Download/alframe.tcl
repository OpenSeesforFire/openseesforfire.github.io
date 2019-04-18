wipe;
#unit: N , mm, s
#					y                  |--------------
#					|                  |--------------
#					|	               |-\-/-\-/--\-/-
#					|_________	x      |-\-/-\-/--\-/-
#				   /				   |-\-/-\-/--\-/-
#				  /                    |--------------
#				 /z                    |--------------

#Load parameters
set Load -200;  			#Point Load on column top: 100N
set UDL -0.005;  			#Uniform Distributed Load on beams: 0.005N/mm
set ThermalLoad 650; 		#Heated part: °C
set HeatTransfer 100; 		#Unheated part: °C

#Model parameters
set nx 10;         			#number of elements along one beam
set ny 6;          			#number of elements for one strorey 
set modelb 250.0;    		#breadth of the frame model
set modelsh 50.0;    		#storey height
set numStorey 10;  			#number of storeys
file mkdir demoData;    	#create a data directory

#Model Domain
model BasicBuilder -ndm 2 -ndf 3;

#Material definition-------------------------------------------//
uniaxialMaterial ElasticThermal 1 68000 2.36e-5;

#Section definition--------------------------------------------//
#Beam
set NumFibersy 10;
set NumFibersz 2;
set b 25.0;
set d 2.0;
set fiberArea [expr ($d/$NumFibersy)*($b/$NumFibersz)];
puts "section 1 for beam";
section fiberSecThermal 1 {
	for {set fiberIDy 1} {$fiberIDy<=[expr $NumFibersy]} {incr fiberIDy} {
		set fiberLocy [expr ($d/$NumFibersy)*($fiberIDy-0.5)-$d/2];
		set fiberLocy1 [expr ($d/$NumFibersy)*($fiberIDy-0.5)];
		for {set fiberIDz 1} {$fiberIDz<=[expr $NumFibersz]} {incr fiberIDz} {
			set fiberLocz [expr ($b/$NumFibersz)*($fiberIDz-0.5)-$b/2];
			#Bottom-Solid Material
			fiber $fiberLocy $fiberLocz $fiberArea 1;
			#fiber $Locy $Locz $fiberArea $MatTag;
			puts "$fiberLocy, $fiberLocz, $fiberArea, 1";
		}
	}
}

#column
set NumFibersy 10;
set NumFibersz 2;
set bc 25.0; 				#local z-direction
set dc 1.0; 				#local y-direction
set fiberArea [expr ($dc/$NumFibersy)*($bc/$NumFibersz)];
puts "section 2 for column";
section fiberSecThermal 2 {
	for {set fiberIDy 1} {$fiberIDy<=$NumFibersy} {incr fiberIDy} {
		set fiberLocy [expr ($dc/$NumFibersy)*($fiberIDy-0.5)-$dc/2];		
		for {set fiberIDz 1} {$fiberIDz<=$NumFibersz} {incr fiberIDz} {
			set fiberLocz [expr ($bc/$NumFibersz)*($fiberIDz-0.5)-$bc/2];		
			fiber $fiberLocy $fiberLocz $fiberArea 1;
			#fiber $Locy $Locz $fiberArea $MatTag;
			puts "$fiberLocy, $fiberLocz, $fiberArea, 1";
		}
	}
}

#Geometric Transformation--------------------------------------------//
#define coordinate transforamtion: geomTransf $type $TransfTag;
#three transformation types can be chosen: Linear, PDelta, Corotational)
#Note: in 2D anlayis local-z direciton is aligned to global Z-driection
geomTransf Corotational 1 ; 

#Generate beam nodes---------------------------------------------//
set EleLenx [expr $modelb/$nx]
for {set storeyID 1} {$storeyID <= $numStorey} {incr storeyID} {
	for {set nxID 0} {$nxID <= $nx} {incr nxID} {
		set locX [expr $nxID*$EleLenx];
		set locY [expr $storeyID*$modelsh];
		set NodeTag [expr $storeyID*100+$nxID+1];
		node $NodeTag $locX $locY;
	}
}
#Beam node tag: 100*storeyID+ 1~(nx+1)
#Generate column nodes---------------------------------------------//
set EleLeny [expr $modelsh/$ny];
set numEley [expr $ny*$numStorey];
for {set nodeID 0} {$nodeID <= $numEley} {incr nodeID} {
	set locX 0;
	set locY [expr $nodeID*$EleLeny];
	set NodeTag [expr 5000+$nodeID+1];
	node $NodeTag $locX $locY;	
}
#Column node tag: 5000+ 1~(ny*No. of storey+1)
puts "nodes generated";

#Generate beam elements---------------------------------------------//
for {set storeyID 1} {$storeyID <= $numStorey} {incr storeyID} {
	for {set nxID 1} {$nxID <=$nx} {incr nxID} {
		set eleID [expr $storeyID*100+$nxID];
		set NodeTag0 [expr $storeyID*100+$nxID];
		set NodeTag1 [expr $storeyID*100+$nxID+1];
		element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 5 1 1;
		#element dispBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag 
	}
}

#Generate column elements---------------------------------------------//
#Column elements are tagged as 5000+1~ny*num_of_storey
for {set nyID 1} {$nyID <=$numEley } {incr nyID} {
		set eleID [expr 5000+$nyID];
		set NodeTag0 [expr 5000+$nyID];
		set NodeTag1 [expr 5000+$nyID+1];
		element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 5 2 1;
		#element dispBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag 
}
puts "elements generated";

#Generate boundary conditions---------------------------------------------//
for {set storeyID 1} {$storeyID <= $numStorey} {incr storeyID} {
	set rightNodeID [expr $storeyID*100+$nx+1];
    fix $rightNodeID 1 1 0;
	
	set leftBeamNodeID [expr $storeyID*100+1];
	set leftColumnNodeID [expr 5000+$storeyID*$ny+1];
	equalDOF $leftBeamNodeID $leftColumnNodeID 1 2;
}
fix 5001 1 1 0;  #Column bottom end


#display 2D deformation shape
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
set  xPixels 450;	# height of graphical window in pixels
set  yPixels 600;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.5;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel2D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

#Data recorder
recorder Node -file demoData/DFreeJntFloor4.out -time -node [expr 5000+$ny*4+1] -dof 1 2 3  disp;
recorder Node -file demoData/DFreeJntFloor5.out -time -node [expr 5000+$ny*5+1] -dof 1 2 3  disp;
recorder Node -file demoData/DFreeJntFloor6.out -time -node [expr 5000+$ny*6+1] -dof 1 2 3  disp;
recorder Node -file demoData/DFreeMidFloor4.out -time -node [expr 400+$nx/2+1] -dof  2   disp;
recorder Node -file demoData/DFreeMidFloor5.out -time -node [expr 500+$ny/2+1] -dof  2  disp;
recorder Node -file demoData/DFreeMidFloor6.out -time -node [expr 600+$ny/2+1] -dof  2  disp;
recorder Element -file demoData/EleSecColumnBot.out -time -ele 5001 section 1  fiber -1 10 stressStrainTangent;
recorder Element -file demoData/EleSecBeamFireExp.out -time -ele [expr 500+$nx/2] section 1  fiber -1 2 stressStrainTangent;
recorder Element -file demoData/EleSecBeamFireUnExp.out -time -ele [expr 700+$nx/2] section 1  fiber -1 2 stressStrainTangent;

#Elements and Nodes
print -file demoData/MeshDomain.out

#Log File
logFile "demoData/Analysis.log"

#Mechanical load
puts "Now applying uniform load"
pattern Plain 1 Linear {
	for {set storeyID 1} {$storeyID < $numStorey} {incr storeyID} {
		for {set nxID 1} {$nxID<= $nx} {incr nxID} {
			set eleID [expr $storeyID*100+$nxID];;
			eleLoad -ele $eleID -type -beamUniform $UDL 0 0;
		}
	}
    set columnTopNode [expr 5000+$numEley+1];
   load	$columnTopNode 0  $Load 0;
};
constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-4 10 2 ;
algorithm Newton;
integrator LoadControl 0.1;
analysis Static;
analyze 10;
loadConst -time 0.0

#Thermal Load Locations
set QuatD [expr $d/4];
set HalfD [expr $d/2];
puts "$QuatD, $HalfD";


puts "Now applying Temperature loading"
# Define Thermal Load
pattern Plain 2 Linear {
eleLoad -range [expr 4*100+1] [expr 4*100+$nx]  -type -beamThermal $ThermalLoad -$HalfD $ThermalLoad -$QuatD $ThermalLoad -0.00001 $HeatTransfer 0.00001 $HeatTransfer $HalfD;
eleLoad -range [expr 5*100+1] [expr 5*100+$nx]  -type -beamThermal $ThermalLoad -$HalfD $ThermalLoad -$QuatD $ThermalLoad -0.00001 $HeatTransfer 0.00001 $HeatTransfer $HalfD;
eleLoad -range [expr 6*100+1] [expr 6*100+$nx]  -type -beamThermal $ThermalLoad -$HalfD $ThermalLoad -$QuatD $ThermalLoad -0.00001 $HeatTransfer 0.00001 $HeatTransfer $HalfD;
}
constraints Plain;
numberer Plain;
system BandGeneral;
test NormUnbalance 1.0e-4 100 1;
algorithm Newton;
integrator LoadControl 0.005;
analysis Static;
analyze 200;
loadConst -time 0.0