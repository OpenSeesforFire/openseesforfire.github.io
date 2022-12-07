# Lim Test 4(Hi-bond composite slabs)
# Simply supported composited stell slab subjected to experimental temperatures
# Units are mm, N, Pa, ton, and C
# Gravity in negative z 
# Script by Qiu Jin

# Problem description
# this model is built by using Method3;
# that is :shell elemetns represent both flat and rib parts. 0.075 m flat shells are centered at z = 0, while 0.13 m Rib shells offset by-0.0275 m.
# A rectangular 4.15 m x 3.15 m comosite slab with 0.13 m total thickness
# The continuous concrete cover has a thickness of 0.075 m , and 198 mm2/m bottom reinforcement is located at 0.02 m above the ribs
# The trapezoidal steel deck has a thickness of 0.00075 m; the width of the upper and lower flange of the rib are 0.132 m and 0.12 m respectively and has a height of 0.0275 m 
# Equivalent rectangular rib is 0.144 m wide and 0.055 m deep. Rib is 0.3 m from center of other ribs.
# There are 11 ribs in total distributing from the center of the rib in the Y direction so the rib span along X direction；
# Supported in z direction on all edges, and x and y translations and z rotation in its four corners
# Load is a uniform load of 5.52 kN/m2 followed by thermal gradient over 3 hour as found in file temp.dat

# define model parameters

# general
set l 4.15  
set w 3.15
set t 0.13
set ct 0.075
set dt 0.00075
set rt 0.055
set lw 0.13;#the lower width of the rib
set uw 0.182;#the upper width of the rib
#set rw [expr 0.5*($lw+$uw)];# the width of the idelised rectangular
#198 mm2/m translates to layer with thickness of 0.198e-3 m
set ro 0.000198
set UDL -5430
set HalfT [expr 0.5*$ct]
#set HalfR [expr 0.5*$rt]

# reinforcement properties
set rfy 565e6
set rEs 2.06e11

# steel decking properties
set dfy 550e6
set dfu 586e6
#set dfu 628e6
set dEs 2.06e11
set vs 0.2

# concrete properties
set fc 32.1e6
set fpc -$fc
set epsc0 -0.0025
#set fpcu [expr 0.05*$fpc] 
#set epscu [expr 10*$epsc0]
set ft [expr 0.1*$fc]
#set epstu 0.002
set Ec [expr 1.5*abs($fpc/$epsc0)]
#set Ets [expr $ft/$epstu]
set v 0.2
#set J [expr $rt*pow($rw,3)/3.0]; #torsional moment
#set G [expr $Ec/(2.0*(1+$v))]; #shear modulus 
#set JG [expr $J*$G] 
set gt [expr $ft/$Ec*$ft*2]
set gc [expr $fc/$Ec*$fc*6]
set lambda 0.15

#mesh
#flat part
set nx 50
# set endNode1 43; #has 43 nodes in the Y-direction
set ny 40
set seleml 0.308; #the length of one composite slab element  
set elemx  [expr $l/$nx]
set elemy1 [expr 0.5*$uw]
set elemy2 [expr 0.5*($seleml-$uw)]
set elemy3 0.098
#set elemy4 $rw
set flatH 0
#puts "$elemx + $elemy1 + $elemy2 + $elemy3"
#set ribH [expr -0.5*$rt- 0.5*$ct]

####################################################################################

# Start up OpenSees model builder and create output folder

wipe
model BasicBuilder -ndm 3 -ndf 6
file mkdir output2

# Load display engine
 source DisplayPlane.tcl
 source DisplayModel2D.tcl
 source DisplayModel3D.tcl

####################################################################################

# Create nodes

# flat part nodes
for {set i 0} {$i <= $nx} {incr i 1} {
	set x [expr $i*$elemx]
	for {set j 0} {$j <= $ny} {incr j 1} {
	    if {$j<=1} {
		set y [expr $j*$elemy3]
	}   elseif {$j==$ny} {
	    set y $w
	}   elseif {$j%4 == 2} {
	    set y [expr $elemy3+(($j-2)/4)*(2*$elemy2+2*$elemy1)+$elemy1]
	}   elseif {$j%4 == 3} {
	    set y [expr $elemy3+(($j-2)/4)*(2*$elemy2+2*$elemy1)+2*$elemy1]
	}   elseif  {$j%4==0} {
	    set y [expr $elemy3+(($j-2)/4)*(2*$elemy2+2*$elemy1)+2*$elemy1+$elemy2]
	}   elseif {$j%4==1} { 
		set y [expr $elemy3+(($j-2)/4)*(2*$elemy2+2*$elemy1)+2*($elemy1+$elemy2)]
	}
		set nodeID [expr int(($j+1)*100+$i+1)]
	   #node nodeTag x  y  z
		node $nodeID $x $y $flatH
		#puts "$nodeID + $x + $y + $flatH"
    }
}

# set any node tags you want to use later here:
set flatmidNode [expr int(($ny/2+1)*100+$nx*0.5+1)];#the mid node in the thick part
#puts "flatmidNode = $flatmidNode"
#puts "flatmidNodeinthin = $flatmidNode1"
#return

####################################################################################

# Create constraints

#fixX x u1  u2  u3 u11 u22 u33
fixX 0  0 	0 	1	0 	0 	0
fixX $l 0 	0 	1	0 	0 	0
fixY 0  0 	0 	1	0 	0 	0
fixY $w 0 	0 	1	0 	0 	0

#fix nodeTag u1  u2  u3 u11 u22 u33
fix  2101    1 	 1 	 0 	 0 	 0 	 1
fix  2151    0 	 0 	 0 	 0 	 0 	 1
#fix  2101    0 	 0 	 0 	 0 	 0 	 1
#fix  2101    0 	 0 	 0 	 0 	 0 	 1
 
####################################################################################

# Create geometric transformation
#geomTransf Corotational $transfTag $vecxzX $vecxzY $vecxzZ 
#geomTransf Corotational    1         0       -1         0

####################################################################################

# Create material models

# Reinforcement
 #uniaxialMaterial SteelECThermal 1  EC2NC  $fy  $Es
  uniaxialMaterial SteelECThermal 1 EC2NC $rfy $rEs; 
  #uniaxialMaterial Steel01Thermal 1      $rfy $rEs 0.1;
 #nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
  nDMaterial PlateRebarThermal 2 	 1 				0
  nDMaterial PlateRebarThermal 3 	 1 				90
 
# Flat part Concrete 
 #nDMaterial  CDPPlaneStressThermal matTag Ec  v  ft  fc  gt  gc
  nDMaterial  CDPPlaneStressThermal 4 	  $Ec $v $ft $fc $gt $gc
 #nDMaterial   PlateFromPlaneStressThermal    matTag   	nDMatTag    forStability
  nDMaterial   PlateFromPlaneStressThermal    5   		4   		1e9
  
# Steel deck 
 #uniaxialMaterial Steel01Thermal $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
  uniaxialMaterial Steel01Thermal 6      $dfy $dEs 0.1;   # Steel deck in the ribs
 #nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
  nDMaterial PlateRebarThermal 7 	  6 				0
  nDMaterial PlateRebarThermal 8 	  6 				90
  
 nDMaterial  J2PlaneStressThermal 9 21 $dEs $vs $dfy $dfu 0.1
 nDMaterial   PlateFromPlaneStressThermal    20   9   2e10;
 
 
####################################################################################

# Create sections

# Flat part section 
set flatsec 1
#section LayeredShellThermal   secTag  numoflayers mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
section  LayeredShellThermal   $flatsec  	   11 		20  $dt	  5 0.01    5 0.01   3 $ro  2 $ro	 5 [expr 0.01-$dt-2*$ro]   5 0.01   5 0.01   5 0.01	  5 0.01	5  0.005; #original&all bottom
#section  LayeredShellThermal   $flatsec  	   11 		20  $dt	 5 10  	5 10   3 $ro  2 [expr 2.4*$ro]	 5 10   5 10   5 10   5 10	  5 10	 5 [expr 5-$dt-3.4*$ro];#only one layer rebar
#section  LayeredShellThermal   $flatsec  	   11 		20  [expr 0.3+$dt]	 5 10  	5 10   3 $ro  2 [expr 2.4*$ro]	 5 [expr 10-($dt+0.3)-3.4*$ro]   5 10   5 10   5 10	  5 10	5 5;#half of web was given to top
# Rib section
set ribsec 2
#section LayeredShellThermal   secTag    -offset           y             numoflayers  mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
#section  LayeredShellThermal   $ribsec   -offset  [expr ($t*0.5-$ct*0.5)]   17 	     7   	 $dt	5  	   0.01   5   0.01   5  0.01   5    0.01   5   0.01   5  0.005  5  0.01  5  0.01   3  $ro 	2 $ro	 5 [expr 0.01-$dt-2*$ro]  5  0.01  5  0.01  5  0.01	 5  0.01	5  0.005;
#section  LayeredShellThermal   $ribsec    21   20 0.54 -92.125   5 7.44 -86.75  8 0.124 -81.75   5 7.9 -76.75  7 0.124 -71.75   5 8.48 -66.75   7 0.124 -61.75    5 9 -56.75     7 0.124 -51.75  5 9.52 -46.75   5 4.2 -39.625   5 10.75 -32.125   5 10 -21.75   3 $ro -16.651   2 $ro -16.453   5 10 -11.354  5 10 -1.354  5 10 8.646   5 10 18.646   5 10 28.646   5 3.854 35.573; #no translation
#section  LayeredShellThermal   $ribsec    -Rib  0.075  0.055  0.0   21   20 0.00054 -0.092125   5 0.00744 -0.08675  8 0.000124 -0.08175   5 0.0079 -0.07675  7 0.000124 -0.07175   5 0.00848 -0.06675   7 0.000124 -0.06175    5 0.009 -0.05675     7 0.000124 -0.05175  5 0.00952 -0.04675   5 0.005 -0.03925   5 0.01 -0.03175   5 0.01 -0.02175   3 $ro -0.016651   2 $ro -0.016453   5 0.01 -0.011354  5 0.01 -0.001354  5 0.01 0.008646   5 0.01 0.018646   5 0.01 0.028646   5 0.003854 0.035573; #original
section  LayeredShellThermal   $ribsec    -Rib  0.075  0.055  0.0   21   20 0.00054 -0.092125   5 0.00744 -0.08675  8 0.000124 -0.08175   5 0.00796 -0.07675  7 0.000124 -0.07175   5 0.00848 -0.06675   7 0.000124 -0.06175    5 0.009 -0.05675     7 0.000124 -0.05175  5 0.00952 -0.04675   5 0.0042 -0.039625   5 0.01075 -0.032125   5 0.01 -0.02175   3 $ro -0.016651   2 $ro -0.016453   5 0.009604 -0.011552  5 0.01 -0.00175  5 0.01 0.00825   5 0.01 0.01825   5 0.01 0.02825   5 0.00425 0.035375; #original
#section  LayeredShellThermal   $ribsec    -Rib  0.075  0.020  0.0   13   20 0.0006 -0.0572   5 0.0094 -0.0522   5 0.01 -0.0425   5 0.01075 -0.032125   5 0.01 -0.02175   3 $ro -0.016651   2 $ro -0.016453   5 0.009604 -0.011552  5 0.01 -0.00175  5 0.01 0.00825   5 0.01 0.01825   5 0.01 0.02825   5 0.00425 0.035375; #original
#section  LayeredShellThermal   $ribsec    -Rib  0.075  0.055  0.0   19   20 0.00054 -0.092125  8 0.000248 -0.091505  7 0.000248 -0.090815  5 0.00744 -0.08675     5 0.00796 -0.07675     5 0.00848 -0.06675     5 0.009 -0.05675     5 0.00952 -0.04675   5 0.0042 -0.039625   5 0.01075 -0.03125   5 0.01 -0.02175   3 $ro -0.016651   2 $ro -0.016453   5 0.009604 -0.011552  5 0.01 -0.00175  5 0.01 0.00825   5 0.01 0.01825   5 0.01 0.02825   5 0.00425 0.035375; #original
#section  LayeredShellThermal   $ribsec    -Rib  75.0  55.0  0.0   20   20 0.54 -92.125   5 7.44 -86.75  8 0.124 -81.75   5 7.9 -76.75  7 0.124 -71.75   5 8.48 -66.75   7 0.124 -61.75    5 9 -56.75     7 0.124 -51.75  5 9.52 -46.75   5 4.2 -39.625   5 10.75 -32.125   5 10 -21.75   3 $ro -16.651     5 10.198 -11.453  5 10 -1.354  5 10 8.646   5 10 18.646   5 10 28.646   5 3.854 35.573; #only one layer
#section  LayeredShellThermal   $ribsec     20 	     20  0.54 -92.12   8 0.124 -91.6775   7 0.124 -91.3325    7 0.124 -90.9875   7 0.124 -90.6425    5 7.44 -86.75   5 7.9 -76.75    5 8.48 -66.75   5 9 -56.75    5 9.52 -46.75   5 4.2 -39.625   5 10.75 -32.125   5 10 -22.75   3 $ro -16.651    5 10.198 -11.453  5  10 -1.354   5 10 8.646   5 10 18.646   5 10 28.646   5 3.854 35.573; #all bottom

####################################################################################

# Create elements

# We can achieve this via the following nested loop
for {set i 1} {$i <= $nx} {incr i 1} {
	for {set j 1} {$j <= $ny} {incr j 1} {
	   if {$j%4==2 || $j%4==3} {
	    set sec $ribsec
      }  else {
        set sec $flatsec
		
      }		
		set node1 [expr int($j*100+$i) ]
		set node2 [expr $node1 +1]
		set node3 [expr int(($j+1)*100+$i+1)]
		set node4 [expr $node3 - 1]	
		
	    set elemID $node1
		#set elemID [expr ]
		
		#element ShellNLDKGQThermal	 $eleTag $iNode $jNode $kNode $lNode $secTag
		 element ShellNLDKGQThermal	 $elemID $node1 $node2 $node3 $node4 $sec
		#puts "$elemID + $node1 + $node2 + $node3 + $node4 + $sec"
	}
}
#return
####################################################################################

# Call graphical engine

set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
 DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

####################################################################################

# Create recorders

recorder Node -file output2/M3midNodeU3.out -time -node $flatmidNode -dof 3 disp

#recorder Element	-file	output2/TempAndE.out	-time	-eleRange	101	4050  material	1	fiber	11	TempAndElong; #flat slab
recorder Element	-file	output2/TempAndE1.out	-time	-eleRange	101	150	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE2.out	-time	-eleRange	201	350	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE3.out	-time	-eleRange	401	550	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE4.out	-time	-eleRange	601	750	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE5.out	-time	-eleRange	801	950	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE6.out	-time	-eleRange	1001	1150	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE7.out	-time	-eleRange	1201	1350	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE8.out	-time	-eleRange	1401	1550	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE9.out	-time	-eleRange	1601	1750	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE10.out	-time	-eleRange	1801	1950	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE11.out	-time	-eleRange	2001	2150	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE12.out	-time	-eleRange	2201	2350	material	1	fiber	21  TempAndElong;
recorder Element	-file	output2/TempAndE13.out	-time	-eleRange	2401	2550	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE14.out	-time	-eleRange	2601	2750	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE15.out	-time	-eleRange	2801	2950	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE16.out	-time	-eleRange	3001	3150	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE17.out	-time	-eleRange	3201	3350	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE18.out	-time	-eleRange	3401	3550	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE19.out	-time	-eleRange	3601	3750	material	1	fiber	11	TempAndElong;
recorder Element	-file	output2/TempAndE20.out	-time	-eleRange	3801	3950	material	1	fiber	21	TempAndElong;
recorder Element	-file	output2/TempAndE21.out	-time	-eleRange	4001	4050	material	1	fiber	11	TempAndElong;

# recorder Element	-file	output2/BottomS1.out	-time	-eleRange	101	4050  material	1	fiber	1	stress;
# recorder Element	-file	output2/BottomS2.out	-time	-eleRange	101	4050  material	1	fiber	2	stress;
# recorder Element	-file	output2/BottomTE.out	-time	-eleRange	101	4050  material	1	fiber	2	TempAndElong;

recorder Element -file output2/RlayerS1.out -time -ele 1926   material 1 fiber 1 TempAndElong;
recorder Element -file output2/RStressS1.out -time	-ele 1926	material  1	fiber 1	stress;
recorder Element -file output2/RStrainS1.out -time	-ele 1926	material  1	fiber 1	strain;
recorder Element -file output2/RlayerC1.out -time -ele 1926   material 1 fiber 2 TempAndElong;
recorder Element -file output2/RStressC1.out -time	-ele 1926	material  1	fiber 2	stress;
recorder Element -file output2/RstrainC1.out -time	-ele 1926	material  1	fiber 2	strain;
recorder Element -file output2/RlayerS2.out -time -ele 1926   material 1 fiber 3 TempAndElong;
recorder Element -file output2/RStressS2.out -time -ele 1926   material 1 fiber 3 stress;
recorder Element -file output2/RstrainS2.out -time -ele 1926   material 1 fiber 3 strain;
recorder Element -file output2/RlayerC2.out -time -ele 1926   material 1 fiber 4 TempAndElong;
recorder Element -file output2/RStressC2.out -time -ele 1926   material 1 fiber 4 stress;
recorder Element -file output2/RstrainC2.out -time -ele 1926   material 1 fiber 4 strain;
recorder Element -file output2/RlayerS3.out -time -ele 1926   material 1 fiber 5 TempAndElong;
recorder Element -file output2/RStressS3.out -time -ele 1926   material 1 fiber 5 stress;
recorder Element -file output2/RStrainS3.out -time -ele 1926   material 1 fiber 5 strain;
recorder Element -file output2/RlayerC3.out -time -ele 1926   material 1 fiber 6 TempAndElong;
recorder Element -file output2/RStressC3.out -time -ele 1926   material 1 fiber 6 stress;
recorder Element -file output2/RstrainC3.out -time -ele 1926   material 1 fiber 6 strain;
recorder Element -file output2/RlayerS4.out -time -ele 1926   material 1 fiber 7 TempAndElong;
recorder Element -file output2/RStressS4.out -time -ele 1926   material 1 fiber 7 stress;
recorder Element -file output2/RStrainS4.out -time -ele 1926   material 1 fiber 7 strain;
recorder Element -file output2/RlayerC4.out -time -ele 1926   material 1 fiber 8 TempAndElong;
recorder Element -file output2/RStressC4.out -time -ele 1926   material 1 fiber 8 stress;
recorder Element -file output2/RstrainC4.out -time -ele 1926   material 1 fiber 8 strain;
recorder Element -file output2/RlayerS5.out -time -ele 1926   material 1 fiber 9 TempAndElong;
recorder Element -file output2/RStressS5.out -time -ele 1926   material 1 fiber 9 stress;
recorder Element -file output2/RStrainS5.out -time -ele 1926   material 1 fiber 9 strain;
recorder Element -file output2/RlayerC5.out -time -ele 1926   material 1 fiber 10 TempAndElong;
recorder Element -file output2/RStressC5.out -time -ele 1926   material 1 fiber 10 stress;
recorder Element -file output2/RstrainC5.out -time -ele 1926   material 1 fiber 10 strain;
recorder Element -file output2/RlayerC6.out -time -ele 1926   material 1 fiber 11 TempAndElong;
recorder Element -file output2/RStressC6.out -time -ele 1926   material 1 fiber 11 stress;
recorder Element -file output2/RstrainC6.out -time -ele 1926   material 1 fiber 11 strain;
recorder Element -file output2/RlayerC7.out -time -ele 1926   material 1 fiber 12 TempAndElong;
recorder Element -file output2/RStressC7.out -time	-ele 1926	material  1	fiber 12	stress;
recorder Element -file output2/RstrainC7.out -time	-ele 1926	material  1	fiber 12	strain;
recorder Element -file output2/RlayerC8.out -time -ele 1926   material 1 fiber 13 TempAndElong;
recorder Element -file output2/RStressC8.out -time	-ele 1926	material  1	fiber 13	stress;
recorder Element -file output2/RstrainC8.out -time	-ele 1926	material  1	fiber 13	strain;
recorder Element -file output2/RlayerRV.out -time -ele 1926   material 1 fiber 14 TempAndElong;
recorder Element -file output2/RStressRV.out -time -ele 1926   material 1 fiber 14 stress;
recorder Element -file output2/RStrainRV.out -time -ele 1926   material 1 fiber 14 strain;
recorder Element -file output2/RlayerRL.out -time -ele 1926   material 1 fiber 15 TempAndElong;
recorder Element -file output2/RStressRL.out -time -ele 1926   material 1 fiber 15 stress;
recorder Element -file output2/RStrainRL.out -time -ele 1926   material 1 fiber 15 strain;
recorder Element -file output2/RlayerC9.out -time -ele 1926   material 1 fiber 16 TempAndElong;
recorder Element -file output2/RStressC9.out -time	-ele 1926	material  1	fiber 16	stress;
recorder Element -file output2/RstrainC9.out -time	-ele 1926	material  1	fiber 16	strain;
recorder Element -file output2/RlayerC10.out -time -ele 1926   material 1 fiber 17 TempAndElong;
recorder Element -file output2/RStressC10.out -time	-ele 1926	material  1	fiber 17	stress;
recorder Element -file output2/RstrainC10.out -time	-ele 1926	material  1	fiber 17	strain;
recorder Element -file output2/RlayerC11.out -time -ele 1926   material 1 fiber 18 TempAndElong;
recorder Element -file output2/RStressC11.out -time	-ele 1926	material  1	fiber 18	stress;
recorder Element -file output2/RstrainC11.out -time	-ele 1926	material  1	fiber 18	strain;
recorder Element -file output2/RlayerC12.out -time -ele 1926   material 1 fiber 19 TempAndElong;
recorder Element -file output2/RStressC12.out -time	-ele 1926	material  1	fiber 19	stress;
recorder Element -file output2/RstrainC12.out -time	-ele 1926	material  1	fiber 19	strain;
recorder Element -file output2/RlayerC13.out -time -ele 1926   material 1 fiber 20 TempAndElong;
recorder Element -file output2/RStressC13.out -time	-ele 1926	material  1	fiber 20	stress;
recorder Element -file output2/RstrainC13.out -time	-ele 1926	material  1	fiber 20	strain;
recorder Element -file output2/RlayerC14.out -time -ele 1926   material 1 fiber 21 TempAndElong;
recorder Element -file output2/RStressC14.out -time	-ele 1926	material  1	fiber 21	stress;
recorder Element -file output2/RstrainC14.out -time	-ele 1926	material  1	fiber 21	strain;



# recorder Element -file output2/ERlayerS1.out -time -ele 1013   material 1 fiber 1 TempAndElong;
# recorder Element -file output2/ERStressS1.out -time	-ele 1013	material  1	fiber 1	stress;
# # recorder Element -file output2/ERStrainS1.out -time	-ele 1013	material  1	fiber 1	strain;
# recorder Element -file output2/ERlayerC2.out -time -ele 1013   material 1 fiber 2 TempAndElong;
# recorder Element -file output2/ERStressC2.out -time	-ele 1013	material  1	fiber 2	stress;
# recorder Element -file output2/ERlayerS2.out -time -ele 1013   material 1 fiber 3 TempAndElong;
# recorder Element -file output2/ERStressS2.out -time -ele 1013   material 1 fiber 3 stress;
# # recorder Element -file output2/ERStrainS2.out -time -ele 1013   material 1 fiber 3 strain;
# recorder Element -file output2/ERlayerC4.out -time -ele 1013   material 1 fiber 4 TempAndElong;
# recorder Element -file output2/ERStressC4.out -time -ele 1013   material 1 fiber 4 stress;
# recorder Element -file output2/ERlayerS3.out -time -ele 1013   material 1 fiber 5 TempAndElong;
# recorder Element -file output2/ERStressS3.out -time -ele 1013   material 1 fiber 5 stress;
# # recorder Element -file output2/ERStrainS3.out -time -ele 1013   material 1 fiber 5 strain;
# recorder Element -file output2/ERlayerS4.out -time -ele 1013   material 1 fiber 7 TempAndElong;
# recorder Element -file output2/ERStressS4.out -time -ele 1013   material 1 fiber 7 stress;
# # recorder Element -file output2/ERStrainS4.out -time -ele 1013   material 1 fiber 7 strain;
# recorder Element -file output2/ERlayerS5.out -time -ele 1013   material 1 fiber 9 TempAndElong;
# recorder Element -file output2/ERStressS5.out -time -ele 1013   material 1 fiber 9 stress;
# # recorder Element -file output2/ERStrainS5.out -time -ele 1013   material 1 fiber 9 strain;
# recorder Element -file output2/ERlayerC12.out -time -ele 1013   material 1 fiber 12 TempAndElong;
# recorder Element -file output2/ERStressC12.out -time	-ele 1013	material  1	fiber 12	stress;
# recorder Element -file output2/ERlayerC13.out -time -ele 1013  material 1 fiber 13 TempAndElong;
# recorder Element -file output2/ERStressC13.out -time	-ele 1013	material  1	fiber 13	stress;
# recorder Element -file output2/ERlayerC14.out -time -ele 1013   material 1 fiber 14 TempAndElong;
# recorder Element -file output2/ERStressC14.out -time	-ele 1013	material  1	fiber 14	stress;
# recorder Element -file output2/ERlayerC15.out -time -ele 1013   material 1 fiber 15 TempAndElong;
# recorder Element -file output2/ERStressC15.out -time	-ele 1013	material  1	fiber 15	stress;
# recorder Element -file output2/ERlayerC16.out -time -ele 1013   material 1 fiber 16 TempAndElong;
# recorder Element -file output2/ERStressC16.out -time	-ele 1013	material  1	fiber 16	stress;
# recorder Element -file output2/ERlayerC17.out -time -ele 1013   material 1 fiber 17 TempAndElong;
# recorder Element -file output2/ERStressC17.out -time	-ele 1013	material  1	fiber 17	stress;
# recorder Element -file output2/ERlayerC18.out -time -ele 1013   material 1 fiber 18 TempAndElong;
# recorder Element -file output2/ERStressC18.out -time	-ele 1013	material  1	fiber 18	stress;
# recorder Element -file output2/ERlayerC19.out -time -ele 1013   material 1 fiber 19 TempAndElong;
# recorder Element -file output2/ERStressC19.out -time	-ele 1013	material  1	fiber 19	stress;
# recorder Element -file output2/ERlayerC20.out -time -ele 1013  material 1 fiber 20 TempAndElong;
# recorder Element -file output2/ERStressC20.out -time	-ele 1013	material  1	fiber 20	stress;
# recorder Element -file output2/ERlayerC21.out -time -ele 1013   material 1 fiber 21 TempAndElong;
# recorder Element -file output2/ERStressC21.out -time	-ele 1013	material  1	fiber 21	stress;
# recorder Element -file output2/ERlayerRV.out -time -ele 1013   material 1 fiber 14 TempAndElong;
# recorder Element -file output2/ERStressRV.out -time -ele 1013   material 1 fiber 14 stress;
# recorder Element -file output2/ERlayerRL.out -time -ele 1013   material 1 fiber 15 TempAndElong;
# recorder Element -file output2/ERStressRL.out -time -ele 1013   material 1 fiber 15 stress;

# recorder Element -file output2/ERlayerS1.out -time -ele 3026   material 1 fiber 1 TempAndElong;
# recorder Element -file output2/ERStressS1.out -time	-ele 3026	material  1	fiber 1	stress;
# recorder Element -file output2/ERStrainS1.out -time	-ele 3026	material  1	fiber 1	strain;
# recorder Element -file output2/ERlayerC2.out -time -ele 3026   material 1 fiber 2 TempAndElong;
# recorder Element -file output2/ERStressC2.out -time	-ele 3026	material  1	fiber 2	stress;
# recorder Element -file output2/ERstrainC2.out -time	-ele 3026	material  1	fiber 2	strain;
# recorder Element -file output2/ERlayerS2.out -time -ele 3026   material 1 fiber 3 TempAndElong;
# recorder Element -file output2/ERStressS2.out -time -ele 3026   material 1 fiber 3 stress;
# recorder Element -file output2/ERStrainS2.out -time -ele 3026   material 1 fiber 3 strain;
# recorder Element -file output2/ERlayerC4.out -time -ele 3026   material 1 fiber 4 TempAndElong;
# recorder Element -file output2/ERStressC4.out -time -ele 3026   material 1 fiber 4 stress;
# recorder Element -file output2/ERstrainC4.out -time -ele 3026   material 1 fiber 4 strain;
# recorder Element -file output2/ERlayerS3.out -time -ele 3026   material 1 fiber 5 TempAndElong;
# recorder Element -file output2/ERStressS3.out -time -ele 3026   material 1 fiber 5 stress;
# recorder Element -file output2/ERStrainS3.out -time -ele 3026   material 1 fiber 5 strain;
# recorder Element -file output2/ERlayerS4.out -time -ele 3026   material 1 fiber 7 TempAndElong;
# recorder Element -file output2/ERStressS4.out -time -ele 3026   material 1 fiber 7 stress;
# recorder Element -file output2/ERStrainS4.out -time -ele 3026   material 1 fiber 7 strain;
# recorder Element -file output2/ERlayerS5.out -time -ele 3026   material 1 fiber 9 TempAndElong;
# recorder Element -file output2/ERStressS5.out -time -ele 3026   material 1 fiber 9 stress;
# recorder Element -file output2/ERStrainS5.out -time -ele 3026   material 1 fiber 9 strain;
# recorder Element -file output2/ERlayerC12.out -time -ele 3026   material 1 fiber 12 TempAndElong;
# recorder Element -file output2/ERStressC12.out -time	-ele 3026	material  1	fiber 12	stress;
# recorder Element -file output2/ERstrainC12.out -time	-ele 3026	material  1	fiber 12	strain;
# recorder Element -file output2/ERlayerC13.out -time -ele 3026  material 1 fiber 13 TempAndElong;
# recorder Element -file output2/ERStressC13.out -time	-ele 3026	material  1	fiber 13	stress;
# recorder Element -file output2/ERstrainC13.out -time	-ele 3026	material  1	fiber 13	strain;
# recorder Element -file output2/ERlayerC14.out -time -ele 3026   material 1 fiber 14 TempAndElong;
# recorder Element -file output2/ERStressC14.out -time	-ele 3026	material  1	fiber 14	stress;
# recorder Element -file output2/ERstrainC14.out -time	-ele 3026	material  1	fiber 14	strain;
# recorder Element -file output2/ERlayerC15.out -time -ele 3026   material 1 fiber 15 TempAndElong;
# recorder Element -file output2/ERStressC15.out -time	-ele 3026	material  1	fiber 15	stress;
# recorder Element -file output2/ERstrainC15.out -time	-ele 3026	material  1	fiber 15	strain;
# recorder Element -file output2/ERlayerC16.out -time -ele 3026   material 1 fiber 16 TempAndElong;
# recorder Element -file output2/ERStressC16.out -time	-ele 3026	material  1	fiber 16	stress;
# recorder Element -file output2/ERstrainC16.out -time	-ele 3026	material  1	fiber 16	strain;
# recorder Element -file output2/ERlayerC17.out -time -ele 3026   material 1 fiber 17 TempAndElong;
# recorder Element -file output2/ERStressC17.out -time	-ele 3026	material  1	fiber 17	stress;
# recorder Element -file output2/ERstrainC17.out -time	-ele 3026	material  1	fiber 17	strain;
# recorder Element -file output2/ERlayerC18.out -time -ele 3026   material 1 fiber 18 TempAndElong;
# recorder Element -file output2/ERStressC18.out -time	-ele 3026	material  1	fiber 18	stress;
# recorder Element -file output2/ERstrainC18.out -time	-ele 3026	material  1	fiber 18	strain;
# recorder Element -file output2/ERlayerC19.out -time -ele 3026   material 1 fiber 19 TempAndElong;
# recorder Element -file output2/ERStressC19.out -time	-ele 3026	material  1	fiber 19	stress;
# recorder Element -file output2/ERstrainC19.out -time	-ele 3026	material  1	fiber 19	strain;
# recorder Element -file output2/ERlayerC20.out -time -ele 3026  material 1 fiber 20 TempAndElong;
# recorder Element -file output2/ERStressC20.out -time	-ele 3026	material  1	fiber 20	stress;
# recorder Element -file output2/ERstrainC20.out -time	-ele 3026	material  1	fiber 20	strain;
# recorder Element -file output2/ERlayerC21.out -time -ele 3026   material 1 fiber 21 TempAndElong;
# recorder Element -file output2/ERStressC21.out -time	-ele 3026	material  1	fiber 21	stress;
# recorder Element -file output2/ERSstrainC21.out -time	-ele 3026	material  1	fiber 21	strain;
# recorder Element -file output2/ERlayerRV.out -time -ele 3026   material 1 fiber 14 TempAndElong;
# recorder Element -file output2/ERStressRV.out -time -ele 3026   material 1 fiber 14 stress;
# recorder Element -file output2/ERstrainRV.out -time -ele 3026   material 1 fiber 14 strain;
# recorder Element -file output2/ERlayerRL.out -time -ele 3026   material 1 fiber 15 TempAndElong;
# recorder Element -file output2/ERStressRL.out -time -ele 3026   material 1 fiber 15 stress;
# recorder Element -file output2/ERstrainRL.out -time -ele 3026   material 1 fiber 15 strain;

# recorder Element -file output2/FlayerS1.out -time -ele 2026   material 1 fiber 1 TempAndElong;
# recorder Element -file output2/FStressS1.out -time	-ele 2026	material  1	fiber 1	stress;
# recorder Element -file output2/FstrainS1.out -time	-ele 2026	material  1	fiber 1	strain;
# recorder Element -file output2/FlayerC1.out -time -ele 2026   material 1 fiber 2 TempAndElong;
# recorder Element -file output2/FStressC1.out	-time	-ele 2026	material  1	fiber 2	stress;
# recorder Element -file output2/FstrainC1.out -time	-ele 2026	material  1	fiber 2	strain;
# recorder Element -file output2/FlayerC2.out -time -ele 2026   material 1 fiber 3 TempAndElong;
# recorder Element -file output2/FStressC2.out	-time	-ele 2026	material  1	fiber 3	stress;
# recorder Element -file output2/FstrainC2.out	-time	-ele 2026	material  1	fiber 3	strain;
# recorder Element -file output2/FlayerRV.out -time -ele 2026   material 1 fiber 4 TempAndElong;
# recorder Element -file output2/FStressRV.out	-time	-ele 2026	material  1	fiber 4	stress;
# recorder Element -file output2/FstrainRV.out	-time	-ele 2026	material  1	fiber 4	strain;
# recorder Element -file output2/FlayerRL.out -time -ele 2026   material 1 fiber 5 TempAndElong;
# recorder Element -file output2/FStressRL.out	-time	-ele 2026	material  1	fiber 5	stress;
# recorder Element -file output2/FstrainRL.out	-time	-ele 2026	material  1	fiber 5	strain;
# recorder Element -file output2/FlayerC3.out -time -ele 2026   material 1 fiber 6 TempAndElong;
# recorder Element -file output2/FStressC3.out	-time	-ele 2026	material  1	fiber 6	stress;
# recorder Element -file output2/FstrainC3.out	-time	-ele 2026	material  1	fiber 6	strain;
# recorder Element -file output2/FlayerC4.out -time -ele 2026   material 1 fiber 7 TempAndElong;
# recorder Element -file output2/FStressC4.out	-time	-ele 2026	material  1	fiber 7	stress;
# recorder Element -file output2/FstrainC4.out	-time	-ele 2026	material  1	fiber 7	strain;
# recorder Element -file output2/FlayerC5.out -time -ele 2026   material 1 fiber 8 TempAndElong;
# recorder Element -file output2/FStressC5.out	-time	-ele 2026	material  1	fiber 8	stress;
# recorder Element -file output2/FstrainC5.out	-time	-ele 2026	material  1	fiber 8	strain;
# recorder Element -file output2/FlayerC6.out -time -ele 2026   material 1 fiber 9 TempAndElong;
# recorder Element -file output2/FStressC6.out	-time	-ele 2026	material  1	fiber 9	stress;
# recorder Element -file output2/FstrainC6.out	-time	-ele 2026	material  1	fiber 9	strain;
# recorder Element -file output2/FlayerC7.out -time -ele 2026   material 1 fiber 10 TempAndElong;
# recorder Element -file output2/FStressC7.out	-time	-ele 2026	material  1	fiber 10 stress;
# recorder Element -file output2/FstrainC7.out	-time	-ele 2026	material  1	fiber 10 strain;
# recorder Element -file output2/FlayerC8.out -time -ele 2026   material 1 fiber 11 TempAndElong;
# recorder Element -file output2/FStressC8.out	-time	-ele 2026	material  1	fiber 11 stress;
# recorder Element -file output2/FstrainC8.out	-time	-ele 2026	material  1	fiber 11 strain;

# recorder Element -file output2/EFlayerS1.out -time -ele 2013   material 1 fiber 1 TempAndElong;
# recorder Element -file output2/EFStressS1.out -time	-ele 2013	material  1	fiber 1	stress;
# recorder Element -file output2/EFstrainS1.out -time	-ele 2013	material  1	fiber 1	strain;
# recorder Element -file output2/EFlayerC2.out -time -ele 2013   material 1 fiber 2 TempAndElong;
# recorder Element -file output2/EFStressC2.out	-time	-ele 2013	material  1	fiber 2	stress;
# recorder Element -file output2/EFstrainC2.out	-time	-ele 2013	material  1	fiber 2	strain;
# recorder Element -file output2/EFlayerC3.out -time -ele 2013   material 1 fiber 3 TempAndElong;
# recorder Element -file output2/EFStressC3.out	-time	-ele 2013	material  1	fiber 3	stress;
# recorder Element -file output2/EFstrainC3.out	-time	-ele 2013	material  1	fiber 3	strain;
# recorder Element -file output2/EFlayerRV.out -time -ele 2013   material 1 fiber 4 TempAndElong;
# recorder Element -file output2/EFStressRV.out	-time	-ele 2013	material  1	fiber 4	stress;
# recorder Element -file output2/EFstrainRV.out	-time	-ele 2013	material  1	fiber 4	strain;
# recorder Element -file output2/EFlayerRL.out -time -ele 2013   material 1 fiber 5 TempAndElong;
# recorder Element -file output2/EFStressRL.out	-time	-ele 2013	material  1	fiber 5	stress;
# recorder Element -file output2/EFstrainRL.out	-time	-ele 2013	material  1	fiber 5	strain;
# recorder Element -file output2/EFlayerC6.out -time -ele 2013   material 1 fiber 6 TempAndElong;
# recorder Element -file output2/EFStressC6.out	-time	-ele 2013	material  1	fiber 6	stress;
# recorder Element -file output2/EFstrainC6.out	-time	-ele 2013	material  1	fiber 6	strain;
# recorder Element -file output2/EFlayerC7.out -time -ele 2013   material 1 fiber 7 TempAndElong;
# recorder Element -file output2/EFStressC7.out	-time	-ele 2013	material  1	fiber 7	stress;
# recorder Element -file output2/EFstrainC7.out	-time	-ele 2013	material  1	fiber 7	strain;
# recorder Element -file output2/EFlayerC8.out -time -ele 2013   material 1 fiber 8 TempAndElong;
# recorder Element -file output2/EFStressC8.out	-time	-ele 2013	material  1	fiber 8	stress;
# recorder Element -file output2/EFstrainC8.out	-time	-ele 2013	material  1	fiber 8	strain;
# recorder Element -file output2/EFlayerC9.out -time -ele 2013   material 1 fiber 9 TempAndElong;
# recorder Element -file output2/EFStressC9.out	-time	-ele 2013	material  1	fiber 9	stress;
# recorder Element -file output2/EFstrainC9.out	-time	-ele 2013	material  1	fiber 9	strain;
# recorder Element -file output2/EFlayerC10.out -time -ele 2013   material 1 fiber 10 TempAndElong;
# recorder Element -file output2/EFStressC10.out	-time	-ele 2013	material  1	fiber 10 stress;
# recorder Element -file output2/EFstrainC10.out	-time	-ele 2013	material  1	fiber 10 strain;
# recorder Element -file output2/EFlayerC11.out -time -ele 2013   material 1 fiber 11 TempAndElong;
# recorder Element -file output2/EFStressC11.out	-time	-ele 2013	material  1	fiber 11 stress;
# recorder Element -file output2/EFstrainC11.out	-time	-ele 2013	material  1	fiber 11 strain;
####################################################################################

# Create ambient loads
set f3r [expr $UDL*$elemx*$elemy1] 
set f3f [expr $UDL*$elemx*$elemy2]
set f3f2 [expr $UDL*$elemx*$elemy3]
pattern Plain 1 Linear {
	for {set i 0} {$i <= $nx} {incr i 1} {
	  for {set j 0} {$j<=$ny} {incr j 1} {
	  if {$j==0 || $j==$ny} {
	  set f3 [expr $f3f2/2]
	 } elseif {$j==1 || $j==39} {
	  set f3 [expr 0.5*($f3r+$f3f2)]
	 } elseif {$j%4==2}  {
      set f3 $f3r
     } elseif {$j%4==1 || $j%4==3} {
	  set f3 [expr 0.5*($f3r+$f3f)]
	 } else {
	  set f3 $f3f
	 }
      set nodeID [expr ($j+1)*100+$i+1]	 
	   #load nodeTag	f1 	f2 	f3 	 f11  f22  f33
		load  $nodeID 	0   0   $f3  0 	  0    0
	}
}
}

####################################################################################

# Analyse ambient loads

set ok 0
set stepSize 0.1
set currentTime 0.0
set tFinal 1

constraints Plain;
numberer Plain;
system BandGeneral;		
test NormUnbalance 10 100 1	
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static
analyze 10;
loadConst -time 0.0

####################################################################################

# Create thermal loads


pattern Plain 2 Linear {

   for {set i 1} {$i <= $nx} {incr i 1} {
	  for {set j 1} {$j<=($ny)} {incr j 1} {
	  if {$j%4==2 || $j%4==3} {
	  set sec "rib"
	 } else {
      set sec "flat"
	 }
	 if {$sec=="rib"} { 	  
	#eleLoad   -ele       $eleTag	     -type -shellThermal -source    "fileName"           $Y1     $Y2   
	eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "M3ribtemp.dat"        -0.0925    0.0375
	#eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "M3ribtemp.dat"        -0.0575    0.0375
	} elseif {$sec=="flat"} { 
    #eleLoad -ele        $eleTag	     -type -shellThermal -source    "fileName"           $Y1     $Y2
	eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "M3flattemp.dat"      -$HalfT   $HalfT
	#eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "75.dat"      -$HalfT   $HalfT
    }
	 set elemd [expr int($j*100+$i)]
	}
}	
}

####################################################################################

# Analyse Thermal loads

set ok 0
set stepSize 10
set currentTime 0.0
set tFinal 10800
#set tFinal 360

constraints Plain;
numberer Plain;
system BandGeneral;			
#test NormUnbalance 8000 1000 1	
test NormDispIncr 1e-3 800 1	
#algorithm KrylovNewton	
algorithm Newton	
integrator LoadControl 15 720 0.001 30;		
#analysis VariableStatic
analysis Static 
#analyze 720 15 0.1 15 4
analyze 285 15;

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
wipe