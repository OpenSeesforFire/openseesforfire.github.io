# simply supported composite slab - validation exercise for section offset
# Simply supported composite slab subjected to two point loads
# Units are mm, s, N, MPa, ton, N.mm, and C
# Gravity in negative z 
# Prepared by Mhd Anwar Orabi

# Problem description:

# A ribbed slab on 0.8 mm steel deck with dimensions of 3 m x 0.65 m x 0.138 m. Pinned on one end and rollered on the other.
# loaded until failure at L/4 from both supports. 
# 50 mm flat concrete layer on top, and no reinforcement other than the decking.  
# Equivalent rectangular rib is 100 mm wide and 88 mm deep. Rib located in center of the slab and 200 mm from center of other ribs.
# This means five central elements are 100 mm each, and two side (flat) elements are 72.5 mm wide.   


# define model parameters

# general
set l 3000
set w 645
set t 138
set dt 0.8
set flat 50
set offset [expr -(0.5*$flat - 0.5*$t)]

# mesh 
set nx 20
set ny 7
set elemx [expr $l/$nx]
set elemy 100
set elemyi 72.5

# steel properties
set fy 375
set fu 400
set Es 210.0e3
set vs 0.2

set Ks [expr $Es/(3*(1 - 2*$vs))]
set Gs [expr $Es/(2*(1 + $vs))]

# concrete properties
set fc 55
set ft [expr 0.1*$fc]
set v 0.2
set epscu 0.0025
set Ec [expr 1.5*$fc/$epscu]
set gt [expr $ft/$Ec*$ft*2]
set gc [expr $fc/$Ec*$fc*6]

####################################################################################

# Start up OpenSees model builder and create output folder

wipe
model BasicBuilder -ndm 3 -ndf 6
file mkdir output
file mkdir output/stress
file mkdir output/strain

# Load display engine
source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl

####################################################################################

# Create nodes
for {set i 0} {$i <= $nx} {incr i 1} {
	set x [expr $i*$elemx]
	for {set j 0} {$j <= $ny} {incr j 1} {
		if {$j <= 1} {
			set y [expr $j*$elemyi]
		} elseif {$j > 1 && $j < $ny} {
			set y [expr $elemyi + ($j - 1)*$elemy]
		} elseif {$j == $ny} {
			set y $w
		}
		set nodeID [expr int($i*100 + $j + 1)]
	   #node nodeTag x  y  z
		node $nodeID $x $y 0
	}
}

# set any node tags you want to use later here:
set loadpt1 [expr int(($l/4)/$elemx)*100]
set loadpt2 [expr int((3*$l/4)/$elemx)*100]
set midpt1 [expr int(($l/2)/$elemx)*100 + 4]
set midpt2 [expr int(($l/2)/$elemx)*100 + 5]
####################################################################################

# Create constraints

#fixX x u1  u2  u3 u11 u22 u33
fixX 0  1 	1 	1	1 	0 	1; # pin
fixX $l 0 	0 	1	1 	0 	1; # roller


 
####################################################################################

# Create geometric transformation
# no geometric transformation requried for shell elements

####################################################################################

# Create material models

# Reinforcement
#uniaxialMaterial SteelECThermal matTag EC2NC yieldStrength YoungsModulus
   # uniaxialMaterial Steel01Thermal 1  $fy $Es 0.01;
   # nDMaterial PlateRebarThermal 3 1 0;
 nDMaterial  J2PlaneStressThermal 2 22 $Es 0.3 $fy $fu 0.01 50;
 nDMaterial   PlateFromPlaneStressThermal    3   2   20e10;
 
#nDMaterial ElasticIsotropic3DThermal $matTag $E0 $Poisson_ratio $Density $Thermal_expansion_ratio <-cSoft/-sSoft>;
 nDMaterial ElasticIsotropic3DThermal 22 	  $Es 0.3 			 0 		  0
 # nDMaterial J2Plasticity 2 $Ks $Gs $fy $fu 0.1 0.01
 
  nDMaterial PlateFiberThermal 23 22
 
# Concrete
#nDMaterial  CDPPlaneStressThermal matTag Ec  v  ft  fc  gt  gc
 nDMaterial  CDPPlaneStressThermal 4 	  $Ec $v $ft $fc $gt $gc
#nDMaterial   PlateFromPlaneStressThermal    matTag   	nDMatTag    forStability
 nDMaterial   PlateFromPlaneStressThermal    5   		4   		1e9
 
  # nDMaterial ElasticIsotropic3DThermal 4 	  $Ec 0.2 			 0 		  0
  # nDMaterial PlateFiberThermal 5 4
 

#nDMaterial DruckerPragerThermal 2 $k $G $sigY $rho $rhoBar $Kinf $K0 $delta1 $delta2 $H $theta $mDen;

####################################################################################

# Create sections
# flat part (50 mm)
set flatSec 1
#section LayeredShellThermal   secTag  		   numoflayers mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
section  LayeredShellThermal   $flatSec  	   10 		   3  		[expr 0.35*$dt]	  5	 [expr 10 - 0.35*$dt]	5  5.0 5 5.0  5 5.0  5 5.0	5  5.0	5 5.0 5  5.0	5 5.0 

# rib (138)
set ribSec 2
#section LayeredShellThermal   secTag  		   numoflayers -offset $offset	mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
section  LayeredShellThermal   $ribSec  	   19 		   -offset $offset 	3  		 $dt 		   5		[expr 10 - $dt]	5 10	5 10	5 10	5 10	5 10	5 10	5 10	5 10	5 10	5  5.0 5 5.0  5 5.0  5 5.0	5  5.0	5 5.0 5  4.0	5 4.0;
####################################################################################

# Create elements
for {set i 0} {$i < $nx} {incr i 1} {
	for {set j 0} {$j < $ny} {incr j 1} {
	
		set elemID [expr int($i*100 + $j + 1)]
		
		set node1 [expr int($i*100 + $j + 1)]
		set node2 [expr int(($i + 1)*100 + $j + 1)]
		set node3 [expr $node2 + 1]
		set node4 [expr $node1 + 1]
		
		if {[expr !fmod($j + 1, 2)]} {
			# element ShellNLDKGQThermal	 $eleTag $iNode $jNode $kNode $lNode $secTag
			element ShellNLDKGQThermal	 $elemID $node1 $node2 $node3 $node4 $ribSec 
		 } else {
		 	# element ShellNLDKGQThermal	 $eleTag $iNode $jNode $kNode $lNode $secTag   $ribSec   $flatSec
			element ShellNLDKGQThermal	 $elemID $node1 $node2 $node3 $node4 $flatSec
		}
	}
}

# set any element tags you want to use later here:
set midrib [expr int(($l/2)/$elemx)*100 + 4]
set midflat [expr int(($l/2)/$elemx)*100 + 3]

####################################################################################

# Call graphical engine

set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 1.0e-5;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels


####################################################################################

# Create recorders
recorder Node -file output/MidU3.out -time -node $midpt1 -dof 3 disp;
recorder Element -file output/deform_rib.out -time -ele $midrib  strains;
recorder Element -file output/deform_flat.out -time -ele $midflat  strains;

recorder Element -file output/stress/rib_STEEL_stress.out -time -ele $midrib  material 1 fiber 1 stress;	
recorder Element -file output/stress/rib_BOT_stress.out -time -ele $midrib  material 1 fiber 2 stress;
recorder Element -file output/stress/rib_03_stress.out -time -ele $midrib  material 1 fiber 3 stress;
recorder Element -file output/stress/rib_04_stress.out -time -ele $midrib  material 1 fiber 4 stress;
recorder Element -file output/stress/rib_05_stress.out -time -ele $midrib  material 1 fiber 5 stress;
recorder Element -file output/stress/rib_06_stress.out -time -ele $midrib  material 1 fiber 6 stress;
recorder Element -file output/stress/rib_07_stress.out -time -ele $midrib  material 1 fiber 7 stress;
recorder Element -file output/stress/rib_08_stress.out -time -ele $midrib  material 1 fiber 8 stress;
recorder Element -file output/stress/rib_09_stress.out -time -ele $midrib  material 1 fiber 9 stress;
recorder Element -file output/stress/rib_10_stress.out -time -ele $midrib  material 1 fiber 10 stress;
recorder Element -file output/stress/rib_11_stress.out -time -ele $midrib  material 1 fiber 11 stress;
recorder Element -file output/stress/rib_12_stress.out -time -ele $midrib  material 1 fiber 12 stress;
recorder Element -file output/stress/rib_13_stress.out -time -ele $midrib  material 1 fiber 13 stress;
recorder Element -file output/stress/rib_14_stress.out -time -ele $midrib  material 1 fiber 14 stress;
recorder Element -file output/stress/rib_TOP_stress.out -time -ele $midrib  material 1 fiber 15 stress;	

recorder Element -file output/stress/flat_STEEL_stress.out -time -ele $midflat  material 1 fiber 1 stress;	
recorder Element -file output/stress/flat_BOT_stress.out -time -ele $midflat  material 1 fiber 2 stress;
recorder Element -file output/stress/flat_03_stress.out -time -ele $midflat  material 1 fiber 3 stress;
recorder Element -file output/stress/flat_04_stress.out -time -ele $midflat  material 1 fiber 4 stress;
recorder Element -file output/stress/flat_05_stress.out -time -ele $midflat  material 1 fiber 5 stress;
recorder Element -file output/stress/flat_06_stress.out -time -ele $midflat  material 1 fiber 6 stress;
recorder Element -file output/stress/flat_07_stress.out -time -ele $midflat  material 1 fiber 7 stress;
recorder Element -file output/stress/flat_08_stress.out -time -ele $midflat  material 1 fiber 8 stress;
recorder Element -file output/stress/flat_09_stress.out -time -ele $midflat  material 1 fiber 9 stress;
recorder Element -file output/stress/flat_TOP_stress.out -time -ele $midflat  material 1 fiber 10 stress;


recorder Element -file output/strain/rib_STEEL_strain.out -time -ele $midrib  material 1 fiber 1 strain;	
recorder Element -file output/strain/rib_BOT_strain.out -time -ele $midrib  material 1 fiber 2 strain;
recorder Element -file output/strain/rib_03_strain.out -time -ele $midrib  material 1 fiber 3 strain;
recorder Element -file output/strain/rib_04_strain.out -time -ele $midrib  material 1 fiber 4 strain;
recorder Element -file output/strain/rib_05_strain.out -time -ele $midrib  material 1 fiber 5 strain;
recorder Element -file output/strain/rib_06_strain.out -time -ele $midrib  material 1 fiber 6 strain;
recorder Element -file output/strain/rib_07_strain.out -time -ele $midrib  material 1 fiber 7 strain;
recorder Element -file output/strain/rib_08_strain.out -time -ele $midrib  material 1 fiber 8 strain;
recorder Element -file output/strain/rib_09_strain.out -time -ele $midrib  material 1 fiber 9 strain;
recorder Element -file output/strain/rib_10_strain.out -time -ele $midrib  material 1 fiber 10 strain;
recorder Element -file output/strain/rib_11_strain.out -time -ele $midrib  material 1 fiber 11 strain;
recorder Element -file output/strain/rib_12_strain.out -time -ele $midrib  material 1 fiber 12 strain;
recorder Element -file output/strain/rib_13_strain.out -time -ele $midrib  material 1 fiber 13 strain;
recorder Element -file output/strain/rib_14_strain.out -time -ele $midrib  material 1 fiber 14 strain;
recorder Element -file output/strain/rib_TOP_strain.out -time -ele $midrib  material 1 fiber 15 strain;	

recorder Element -file output/strain/flat_STEEL_strain.out -time -ele $midflat  material 1 fiber 1 strain;	
recorder Element -file output/strain/flat_BOT_strain.out -time -ele $midflat  material 1 fiber 2 strain;
recorder Element -file output/strain/flat_07_strain.out -time -ele $midflat  material 1 fiber 7 strain;
recorder Element -file output/strain/flat_08_strain.out -time -ele $midflat  material 1 fiber 8 strain;
recorder Element -file output/strain/flat_09_strain.out -time -ele $midflat  material 1 fiber 9 strain;
recorder Element -file output/strain/flat_TOP_strain.out -time -ele $midflat  material 1 fiber 10 strain;


####################################################################################

# Create ambient loads
pattern Plain 1 Linear {
	   #load nodeTag	f1 	f2 	f3 	 f11  f22  f33
	   	load  [expr $loadpt1 + 3] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt1 + 4] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt1 + 5] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt1 + 6] 		0   0   -250.0  0 	  0    0
		
		load  [expr $loadpt2 + 3] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt2 + 4] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt2 + 5] 		0   0   -250.0  0 	  0    0
		load  [expr $loadpt2 + 6] 		0   0   -250.0  0 	  0    0
}

####################################################################################

# Analyse ambient loads

set ok 0
set stepSize 0.1
set currentTime 0.0
set tFinal 40.0

constraints Plain    		
numberer Plain			
system BandGeneral		
#test RelativeNormUnbalance 0.05 1000 2
test NormDispIncr 0.01 1000 2	
algorithm ModifiedNewton	
#integrator ArcLength 1 0.01			
integrator LoadControl $stepSize 1000 0.001 0.2				
analysis Static
analyze 400

# integrator LoadControl 0.001			
# analysis Static
# analyze 1000

# integrator DisplacementControl $midpt1 3 -0.2 	
# analysis Static
# analyze 500


# while {$currentTime < $tFinal && $ok == 0} {
	# puts ""; puts ""; 
	# puts "-------------------------------------------------------"
	# puts "-------------------------------------------------------"
	# puts "";
	# puts "";
	# puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	# puts ""; puts "";
	# set ok [analyze 1];
	# if {$ok == 0} {
		# set currentTime [getTime]
		# puts "Analysis step successful. Current load is $currentTime kN."
		# set centralDisp [nodeDisp $midpt1 3]
		# puts "Central displacement is $centralDisp mm"
		# puts "-------------------------------------------------------"
		# set ribFib15 [eleResponse $midrib material 1 fiber 15 stress]
		# set ribFib02 [eleResponse $midrib material 1 fiber 2 stress]
		# set ribSteels [eleResponse $midrib material 1 fiber 1 stress]
		
		# set ribFib15strain [eleResponse $midrib material 1 fiber 15 strain]
		# set ribFib02strain [eleResponse $midrib material 1 fiber 2 strain]
		# set ribSteelstrain [eleResponse $midrib material 1 fiber 1 strain]
		# puts "_____________________________RIB_____________________________"
		# puts "Stress in the rib's top conc. : ($ribFib15) MPa"
		# puts "Strain in the rib's top conc. : ($ribFib15strain)"
		# puts "";
		# puts "Stress in the rib's bot. conc : ($ribFib02) MPa"
		# puts "Strain in the rib's bot. conc : ($ribFib02strain)"
		# puts "";
		# puts "Stress in the rib's decking is: ($ribSteels) MPa"
		# puts "Strain in the rib's decking is: ($ribSteelstrain)"
		# puts "";
		
		# set flatFib06 [eleResponse $midflat material 1 fiber 6 stress]
		# set flatFib02 [eleResponse $midflat material 1 fiber 2 stress]
		# set flatSteels [eleResponse $midflat material 1 fiber 1 stress]
		
		# set flatFib06strain [eleResponse $midflat material 1 fiber 6 strain]
		# set flatFib02strain [eleResponse $midflat material 1 fiber 2 strain]
		# set flatSteelstrain [eleResponse $midflat material 1 fiber 1 strain]
		# puts "_____________________________FLAT_____________________________"
		# puts "Stress in the flat's top conc. : ($flatFib06) MPa"
		# puts "Strain in the flat's top conc. : ($flatFib06strain)"
		# puts ""
		# puts "Stress in the flat's bot. conc : ($flatFib02) MPa"
		# puts "Strain in the flat's bot. conc : ($flatFib02strain)"
		# puts ""
		# puts "Stress in the flat's decking is: ($flatSteels) MPa"
		# puts "Strain in the flat's decking is: ($flatSteelstrain)"
	# } else {
		# puts "Analysis failed at load = $currentTime kN"
		# break
	# }
# }
# if {$currentTime == $tFinal} {
	# puts "Analysis completed."
# }

loadConst -time 0.0

#file copy ?-force? ?- -? source target
####################################################################################



# constraints Plain    		
# numberer Plain			
# system BandGeneral		
# algorithm ModifiedNewton							
# analysis Static


# set ok 0.0; 			# indicating state of the analysis
# set currentTime 0.0; 	# time within the fire
# set deltaT 0.05;			# time step for the fire
# set deltaTt $deltaT;	# time step for the fire

# integrator LoadControl $deltaT;
# test RelativeNormUnbalance 0.05 1000 1

# while {$ok == 0 && $currentTime < $tFinal && $deltaT > 0.0001} {
	# puts "Attempting analysis for time [expr $currentTime + $deltaT], with time step $deltaT."
	# set ok [analyze 1]
	# if {$ok != 0} {
		# while {$ok != 0 && $deltaTt >= $deltaT*pow(0.5,5)} { 
			# puts "Analysis did not converge with deltaTt = $deltaTt. Halving timestep to [expr $deltaTt/2.0]."
			# set deltaTt [expr $deltaTt/2.0]
			# integrator LoadControl $deltaTt;
			# puts "Attempting analysis for time [expr $currentTime + $deltaTt]"
			# set ok [analyze 1]			
		# }
		# set deltaT $deltaTt
	# }
	# if {$ok != 0 || $deltaT <= 0.0001} {
		# if {$deltaT <= 0.0001} {
			# puts "Analysis failed at time [getTime] because increment is too small at $deltaT. Last successful step was $currentTime"
		# } else {
			# puts "Analysis failed at time [getTime]. Last successful step was $currentTime"
		# }
	# }	elseif {$ok == 0} {
			# set currentTime [getTime]
			# puts "Analysis for this step succeeded. Current time is $currentTime"
	# };
		
# };
# if {$currentTime == $tFinal} {
# puts "Analysis succeeded - current time is $currentTime"
# };
wipe;