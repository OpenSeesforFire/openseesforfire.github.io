wipe;
#units are mm, N, MPa, ton, and C

# Column coordinates
set E1x -9000.
set E1y 0.
set E2x -9000.
set E2y 6000.
set F2x 0.
set F2y 6000.
set F1x 0.
set F1y 0.
set F12x 0.
set F12y 3000.
set E12x -9000.
set E12y 3000.
set EF12x -3900.
set EF12y 7500.
set xlim -10000.
set ylim 7500.
set xlim2 -14400.

set h_rib 70.
set h_slab 70.
set h_UB356 355.
set h_UB305 303.4
set h_UB610 602.6 

#model builder
model BasicBuilder -ndm 3 -ndf 6;
#nodes####################################################################################################################

#horizontal elements
set rib_space 300
set horiz_elem_length $rib_space
set vert_elem_length 500
set num_elems_horiz [expr int(abs($F1x - $xlim2)/$horiz_elem_length)]
set num_elems_vert [expr int(abs($F1y - $EF12y)/$vert_elem_length) ]
set num_elems_ribs [expr int($num_elems_vert*$num_elems_horiz)] 

#FE1, FE12, FE2
for {set i 1} {$i <= [expr $num_elems_horiz+1]} {incr i 1} {
	set x [expr -($i - 1.0)*($horiz_elem_length)];
	set y $F12y
	set y2 $F2y
	node $i $x 0 [expr -0.5*$h_UB356 - $h_rib -0.5*$h_slab];
	node [expr 100 + $i] $x $y [expr -0.5*$h_UB305 - $h_rib -0.5*$h_slab];
	node [expr 200 + $i] $x $y2 [expr -0.5*$h_UB610 - $h_rib -0.5*$h_slab];
};

#F12, E12x
for {set i 1} {$i <= [expr $num_elems_vert+1]} {incr i 1} {
	set y [expr ($i - 1.0)*($vert_elem_length )];
	set x $F1x
	set x1 $E1x
	node [expr 300 + $i] 0 $y [expr -0.5*$h_UB356 - $h_rib -0.5*$h_slab];
	node [expr 400 + $i] $x1 $y [expr -0.5*$h_UB356 - $h_rib -0.5*$h_slab];
};

#Ribs NEED to make sure that E12 doesn't get a rib through it
for {set i 1} {$i < [expr $num_elems_horiz]} {incr i 1} {
	
	if {abs([expr $i*$horiz_elem_length]) >= abs($E1x)} {
	set x [expr -($i+1)*($horiz_elem_length)];
	} else {
	set x [expr -($i)*($horiz_elem_length)];
	};
	for {set j 1} {$j <= [expr $num_elems_vert+1]} {incr j 1} {
		set y [expr ($j - 1)*($vert_elem_length )];
		node [expr $i*1000 + $j] $x $y [expr - 0.5*$h_rib -0.5*$h_slab];
	};
};

#rigid links
# horizontal beams and ribs
set type bar
for {set i 2} {$i < $num_elems_horiz} {incr i 1} {
	# EF1
	set masterNode1 [expr int($i)];
	set slaveNode1 [expr int(($i-1)*1000 + $F1y/$vert_elem_length + 1)]
	rigidLink $type $masterNode1 $slaveNode1;
	
	#EF12
	set masterNode2 [expr int($i+100)]
	set slaveNode2 [expr int(($i-1)*1000 + $F12y/$vert_elem_length + 1)]
	rigidLink $type $masterNode2 $slaveNode2;
	
	#EF2
	set masterNode3 [expr int($i+200)];
	set slaveNode3 [expr int(($i-1)*1000 + $F2y/$vert_elem_length + 1)]
	rigidLink $type $masterNode3 $slaveNode3;
}

# slab and ribs




#Material props
set fpc -35
set epsc0 -0.003
set fpcu [expr $fpc*0.05];
set epsU -0.02
set lambda 0.1
set ft 3
set epstu 0.002
set Ets [expr $ft/$epstu];
set Ec [expr abs(1.5*$fpc/$epsc0)]

uniaxialMaterial Steel01Thermal 1 400 2e5 0.01;  								# S355 - range of 371 to 413 MPa
uniaxialMaterial Steel01Thermal 2 300 2e5 0.01;									# S275 - range of 291 - 318 MPa
uniaxialMaterial Concrete02Thermal 3 $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets;	# Grade 35 LW Concrete
uniaxialMaterial Steel01Thermal 4 350 2e5 0.01;									# Decking steel - minimum of 350 MPa		 						

set gt [expr abs($ft/$Ec*$ft*6)];
set gc [expr abs($fpc/$Ec*$fpc*6)];

nDMaterial   CDPPlaneStressThermal 5 $Ec 0.2 $ft [expr abs($fpc)] $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    6   5   10e3;

#nDMaterial J2Plasticity 		 $matTag $K  	$G  		$sig0   $sigInf $delta   $H
nDMaterial  J2PlaneStressThermal 7 	     22  	2e5  0.3 	485 	550 	10  	 2e4; # Must check with liming how this works
nDMaterial  PlateFromPlaneStressThermal    8   7   2.0e5; 							      # Steel mesh - unknown strength for the A142 anti-crack mesh						       
nDMaterial PlateRebarThermal 9 4 0; # Decking steel - minimum of 350 MPa

	
#slab section
#section LayeredShell         $sectionTag $nLayers $matTag1 $thickness1...$matTagn $thicknessn
section LayeredShellThermal   5  		  9	   9 0.5   6 10	  6 10   8 0.141   6 10   6 10    6 10   6 10   6 10; # MUST double check slab thickness!!

set e1 100000; # first element
set n1 100000; # first node
block2D $num_elems_horiz $num_elems_vert $e1 $n1 ShellNLDKGQThermal  5 { 
    1   0. 			0.	 	 0.
    2   0. 			7500.    0.
    3   -14400.  	7500. 	 0.
    4   -14400. 	0.		 0.
};

#Sections 
#UB356x171x51
set MatTag 1 
set SecTag 1 
section FiberThermal $SecTag { 
#fiber       $yLoc   $zLoc          $A  		  $matTag
 fiber		 0.		 175.583300     657.41670000  $MatTag
 fiber 		 0. 	 171.75000000   657.41670000  $MatTag
 fiber 		 0. 	 167.91670000   657.41670000  $MatTag
 fiber 		 0. 	 132.80000000   491.36000000  $MatTag
 fiber 		 0. 	 66.40000000    491.36000000  $MatTag
 fiber 		 0. 	 0.00000000     491.36000000  $MatTag
 fiber 		 0. 	 -66.40000000   491.36000000  $MatTag
 fiber 		 0. 	 -132.80000000  491.36000000  $MatTag
 fiber 		 0. 	 -167.91670000  657.41670000  $MatTag
 fiber 		 0. 	 -171.75000000  657.41670000  $MatTag
 fiber 		 0. 	 -175.58330000  657.41670000  $MatTag
}

#UB305x165x40
set MatTag 2
set SecTag 2
section FiberThermal $SecTag { 
#fiber       $yLoc   $zLoc          $A  			   $matTag
 fiber 		 0.  	150.00000000    561.00000000       $MatTag
 fiber		 0.  	146.60000000    561.00000000       $MatTag
 fiber 		 0.  	143.20000000    561.00000000       $MatTag
 fiber 		 0. 	113.20000000    339.60000000       $MatTag
 fiber 		 0.  	0.00000000      339.60000000       $MatTag
 fiber 		 0. 	-56.60000000    339.60000000       $MatTag
 fiber 		 0. 	-113.20000000   339.60000000       $MatTag
 fiber 		 0. 	-143.20000000   561.00000000       $MatTag
 fiber 		 0. 	-146.60000000   561.00000000       $MatTag
 fiber 		 0. 	-150.00000000   561.00000000       $MatTag
}

#UB610x229x101
set MatTag 2
set SecTag 3
section FiberThermal $SecTag { 
#fiber       $yLoc   $zLoc          $A  		   $matTag
 fiber 		 0.  	 298.83330000 	1122.82670000  $MatTag
 fiber 		 0.  	 293.90000000 	1122.82670000  $MatTag
 fiber 		 0.  	 288.96670000 	1122.82670000  $MatTag
 fiber 		 0.  	 229.20000000 	1203.30000000  $MatTag
 fiber 		 0.  	 114.60000000 	1203.30000000  $MatTag
 fiber 		 0.  	 0.00000000   	1203.30000000  $MatTag
 fiber 		 0. 	 -114.60000000  1203.30000000  $MatTag
 fiber 		 0. 	 -229.20000000  1203.30000000  $MatTag
 fiber 		 0. 	 -288.96670000  1122.82670000  $MatTag
 fiber 		 0. 	 -293.90000000  1122.82670000  $MatTag
 fiber 		 0. 	 -298.83330000  1122.82670000  $MatTag
}

#ribs 300 mm wide, 70 mm deep
set MatTag 3
set SecTag 4
set yI -150
set yJ 150
set zI -35
set zJ 35
set t 0.9

section FiberThermal $SecTag {
   #patch rect $matTag $numSubdivY $numSubdivZ  $yI $zI $yJ $zJ
	patch rect $MatTag 1 		   5			$yI  $zI $yJ $zJ
	#fiber       $yLoc   $zLoc          $A  		   					$matTag # steel deck
	fiber 		 0 		 $zI 			[expr $t*(abs($yJ)+abs($yI))]   4
    fiber 		 $yI 	 0 			    [expr $t*(abs($zJ)+abs($zI))]   4
	fiber 		 $yJ 	 0 			    [expr $t*(abs($zJ)+abs($zI))]   4	
};

#Transformation
geomTransf Corotational 1 1 0 0;
geomTransf Corotational 2 0 1 0;
set VertTransform 1
set HoriTransform 2

#ELES FOR BEAMS

for {set i 1} {$i <= $num_elems_horiz} {incr i 1} {

	# Beams FE1, ED1
	set secTag 1
	set node1 $i
	set node2 [expr $i+1]
	set eleID $i	
	element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $HoriTransform;
	
	#Beams FE12, ED12
	set secTag 2
	set node1 [expr 100+$i]
	set node2 [expr 100+$i+1]
	set eleID [expr 100+$i+1]	
	element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $HoriTransform;
	
	#Beams FE2, ED2
	set secTag 3
	set node1 [expr 200+$i]
	set node2 [expr 200+$i+1]
	set eleID [expr 200+$i+1]	
	element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $HoriTransform;
}

for {set i 1} {$i <= $num_elems_vert} {incr i 1} {
	
	#Beams F 12 and F 23
	set secTag 1
	set node1 [expr 300+$i]
	set node2 [expr 300+$i+1]
	set eleID [expr 300+$i+1]	
	element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $VertTransform;
	
	#Beams E 12 and E 23
	set secTag 1
	set node1 [expr 400+$i]
	set node2 [expr 400+$i+1]
	set eleID [expr 400+$i+1]	
	element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $VertTransform;
}

#Ribs
set secTag 4
for {set i 1} {$i < [expr $num_elems_horiz]} {incr i 1} {
		for {set j 1} {$j <= [expr $num_elems_vert]} {incr j 1} {			
			set node1 [expr $i*1000 + $j]
			set node2 [expr $i*1000 + $j + 1]
			set eleID [expr $i*1000 + $j]	
			element dispBeamColumnThermal $eleID $node1 $node2 5 $secTag $VertTransform
		};
};

source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl

set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

print nodes.txt -node
print elements.txt -ele

if {0} {



#reiforcement is 2 no.10 rebars located 20 mm from top and bottom of the beam

set cc 0.02;
set rebar_diameter 0.01; 

#reinforcement layer properties
set rebars_area [expr $n*$rebar_diameter*$rebar_diameter*3.14/4.0];
set rebar_layer_t [expr $rebars_area/$b]; 

#cover layer properties
set cc_area [expr $cc*$b];

#mesh refinement
set no_of_elements 30; 
set no_of_fibers 3;
set num_int_pts 5;

#load
set ThermalAnalysis 0;
set P_load 2.5e3;
set T_top 0.0;  
set T_bot 1000;

#material properties
#reinforcing steel
set fy 345e6;
set E0 2.06e11;

#concrete
set fpc -30.0e6; #concrete compressive strength at 28 days (compression is negative)
set epsc0 -0.003; #concrete strain at maximum strength
set fpcu [expr $fpc*0.05]; #concrete crushing strength
set epsU -0.02; #concrete strain at crushing strength*
set lambda 0.1; #ratio between unloading slope at $epscu and initial slope
set ft 3.0e6; #tensile strength
set Ets [expr -$ft/0.005]; #tension softening stiffness (absolute value) (slope of the linear tension softening branch)






#materials################################################################################################################
#reinforcing steel
#uniaxial rebar material
uniaxialMaterial SteelECThermal 1 EC2NH $fy $E0; #eurocode 2 (vs. eurocode 3) with no hardening
#creating a plate fiber section from the uniaxial steel rebar material
#nDmaterial PlateRebar 		  $newmatTag $matTag $sita http://www.luxinzheng.net/download/OpenSEES/En_THUShell_OpenSEES.htm
nDMaterial 	PlateRebarThermal 90 		 1 		 90; #longitudinal steel; this can now be used in the multilayered shell section

#concrete
#                                                               fc       ft      fcu     epsc0     epscu  epstu  stc
#nDMaterial PlaneStressUserMaterial   2       40         7    $fpc   	 $ft     $fpcu   $epsc0    $epsU  0.001  0.08
set gt [expr 3.0e6/1.79e10*3.0e6*2];
set gc [expr 30e6/1.79e10*30e6*6];
nDMaterial  CDPPlaneStressThermal 100  1.79e10 0.2  3.0e6  30e6 $gt $gc;
nDMaterial   PlateFromPlaneStressThermal    2   100  10e9;

#section definition########################################################################################################
set section_conc_depth [expr ($t)-2*$cc-2*$rebar_diameter];
set layer_depth [expr $section_conc_depth/$no_of_fibers];
set fiber_area [expr $layer_depth*$b]
#section LayeredShell $sectionTag $nLayers $matTag1 $thickness1...$matTagn $thicknessn
section LayeredShellThermal 1 [expr 4+$no_of_fibers] 2 $cc 90 $rebar_layer_t 2 $layer_depth 2 $layer_depth 2 $layer_depth 90 $rebar_layer_t 2 $cc

#elements##################################################################################################################
for {set i 1} {$i <= $no_of_elements} {incr i} {
   #element ShellMITC4Thermal $eleTag $iNode $jNode 		$kNode 			$lNode 			$secTag
	element ShellNLDKGQThermal $i  	  $i 	 [expr $i + 1] 	[expr $i + 101]	[expr $i + 100]	1;		
}

#recorder##################################################################################################################
recorder Node -file EndDisp.out -time -node [expr $no_of_elements+1] -dof 2 disp;

#nodal_load################################################################################################################
pattern Plain 1 Linear {
load [expr $no_of_elements+1] 0. [expr $P_load/2.0] 0. 0. 0. 0.
load [expr $no_of_elements+101] 0. [expr $P_load/2.0] 0. 0. 0. 0.
}

if {$ThermalAnalysis} { 
pattern Plain 2 Linear {	
for {set i 1} {$i <=$no_of_elements} {incr i} { 
	eleLoad -ele $i -type -beamThermal $T_bot [expr -$t/2] $T_top [expr $t/2];
	}
}
}

#analysis_object###########################################################################################################
#recorder plot EndDisp.out End_Displacement 10 10 400 400 -columns 1 2
record

constraints Transformation;     		
numberer Plain;			
system BandGeneral;		
test NormDispIncr 1e-2 1000 1;	
algorithm Newton;					
integrator LoadControl 0.01;				
analysis Static;
analyze 100;					
loadConst
wipe}
