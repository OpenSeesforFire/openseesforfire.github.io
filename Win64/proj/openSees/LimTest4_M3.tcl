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
set l 4150  
set w 3150
set t 130
set ct 75
set dt 0.75
set rt 55
set lw 130;#the lower width of the rib
set uw 182;#the upper width of the rib
#set rw [expr 0.5*($lw+$uw)];# the width of the idelised rectangular
#198 mm2/m translates to layer with thickness of 0.198e-3 m
set ro 0.198
set UDL -5.43e-3
set HalfT [expr 0.5*$ct]
#set HalfR [expr 0.5*$rt]

# reinforcement properties
set rfy 565
set rEs 2.1e5

# steel decking properties
set dfy 550
#set dfu 586
set dfu 628
set dEs 2.1e5
set vs 0.3


# concrete properties
set fc 32.1
set fpc -$fc
set epsc0 -0.0025
#set fpcu [expr 0.05*$fpc] 
#set epscu [expr 10*$epsc0]
set ft [expr 0.1*$fc]
#set ft 3.3
#set epstu 0.002
#set Ec [expr 1.5*abs($fpc/$epsc0)]
set Ec 2.2e4
#set Ets [expr $ft/$epstu]
set v 0.2
#set J [expr $rt*pow($rw,3)/3.0]; #torsional moment
#set G [expr $Ec/(2.0*(1+$v))]; #shear modulus 
#set JG [expr $J*$G] 
set gt [expr $ft/$Ec*$ft*4]
set gc [expr $fc/$Ec*$fc*6]
set lambda 0.15


#mesh
#flat part
set nx 50
# set endNode1 43; #has 43 nodes in the Y-direction
set ny 40
set seleml 308; #the length of one composite slab element  
set elemx  [expr $l/$nx]
set elemy1 [expr 0.5*$uw]
set elemy2 [expr 0.5*($seleml-$uw)]
set elemy3 98.0
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
# We want to number our nodes in this order:
# 4101 4102	 4103	4104	4105   4106	 4107   4108  4109  4110  4121
#  .	.	  .	      .	      .	    .	  .       .	    .     .	    .
#  .	.	  .	      .	      .	    .	  .       .	    .     .	    .
#  .	.	  .	      .	      .	    .	  .       .	    .     .	    .
#  .	.	  .	      .	      .	    .	  .       .	    .     .	    .
#  101	102	 103	104 	105	   106	 107    108	  109   110	  121

# We can achieve this via the following nested loop
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
set flatmidNode [expr int((($ny+1)/2+1)*100+($nx*0.5+1))];#the mid node in the thick part
# set flatmidNode1 [expr int(($endNode1/2+1)*100+($nx*0.5+1)-2)];#the mid node in the thin  part
puts "flatmidNode = $flatmidNode"
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
fix  101   	 1 	 1 	 0 	 0 	 0 	 1
# fix  4151    1 	 1 	 0 	 0 	 0 	 1
# fix  4101    1 	 1 	 0 	 0 	 0 	 1
# fix  151  	 1 	 1 	 0 	 0 	 0 	 1
 
####################################################################################

# Create geometric transformation
#geomTransf Corotational $transfTag $vecxzX $vecxzY $vecxzZ 
#geomTransf Corotational    1         0       -1         0

####################################################################################

# Create material models

# Reinforcement
 #uniaxialMaterial SteelECThermal 1  EC2NC  $fy  $Es
  uniaxialMaterial SteelECThermal 1 EC2NC $rfy $rEs
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
  uniaxialMaterial Steel01Thermal 6      $dfy $dEs 0;   # Steel deck in the ribs
 #nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
  nDMaterial PlateRebarThermal 7 	  6 				0
  nDMaterial PlateRebarThermal 8 	  6 				90
  
 nDMaterial  J2PlaneStressThermal 9 21 $dEs $vs $dfy $dfu 0.027 0;
 nDMaterial   PlateFromPlaneStressThermal    20   9   2e10;
 
 
####################################################################################

# Create sections

# Flat part section 
set flatsec 1
#section LayeredShellThermal   secTag  numoflayers mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
section  LayeredShellThermal   $flatsec  	   11 		20  $dt	  5 10    5 10   3 $ro  2 $ro	 5 10   5 10   5 10   5 10	  5 10	5  [expr 5-$dt-2*$ro]; #original&all bottom
#section  LayeredShellThermal   $flatsec  	   11 		20  $dt	 5 10  	5 10   3 $ro  2 [expr 2.4*$ro]	 5 10   5 10   5 10   5 10	  5 10	 5 [expr 5-$dt-3.4*$ro];#only one layer rebar
#section  LayeredShellThermal   $flatsec  	   11 		20  [expr 0.3+$dt]	 5 10  	5 10   3 $ro  2 [expr 2.4*$ro]	 5 [expr 10-($dt+0.3)-3.4*$ro]   5 10   5 10   5 10	  5 10	5 5;#half of web was given to top
# Rib section
set ribsec 2
#section LayeredShellThermal   secTag    -offset           y             numoflayers  mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
#section  LayeredShellThermal   $ribsec   -offset  [expr ($t*0.5-$ct*0.5)]   17 	     7   	 $dt	5  	   0.01   5   0.01   5  0.01   5    0.01   5   0.01   5  0.005  5  0.01  5  0.01   3  $ro 	2 $ro	 5 [expr 0.01-$dt-2*$ro]  5  0.01  5  0.01  5  0.01	 5  0.01	5  0.005;
#section  LayeredShellThermal   $ribsec    21    20 0.54 -92.125   5 7.44 -86.75  8 0.124 -81.75   5 7.9 -76.75  7 0.124 -71.75   5 8.48 -66.75   7 0.124 -61.75    5 9 -56.75     7 0.124 -51.75  5 9.52 -46.75   5 4.2 -39.625   5 10.75 -32.125   5 10 -21.75   3 $ro -16.651   2 $ro -16.453   5 10 -11.354  5 10 -1.354  5 10 8.646   5 10 18.646   5 10 28.646   5 3.854 35.573; #no translation
section  LayeredShellThermal   $ribsec    -Rib  75.0  55.0  0.0   21   20 0.54 -92.125   5 7.44 -86.75  8 0.124 -81.75   5 7.9 -76.75  7 0.124 -71.75   5 8.48 -66.75   7 0.124 -61.75    5 9 -56.75     7 0.124 -51.75  5 9.52 -46.75   5 4.2 -39.625   5 10.75 -32.125   5 10 -21.75   3 $ro -16.651   2 $ro -16.453   5 10 -11.354  5 10 -1.354  5 10 8.646   5 10 18.646   5 10 28.646   5 3.854 35.573; #original
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
				
		set node4 [expr int($j*100+$i) ]
		set node3 [expr $node4 + 100]
		set node2 [expr $node3 +1]
		set node1 [expr $node4 + 1]
	    set elemID $node4
		
		#element ShellNLDKGQThermal	 $eleTag $iNode $jNode $kNode $lNode $secTag
		 element ShellNLDKGQThermal	 $elemID $node1 $node2 $node3 $node4 $sec
		#puts "$elemID + $node1 + $node2 + $node3 + $node4 + $sec"
	}
}
# return
####################################################################################

# Call graphical engine

set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.00001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

####################################################################################

# Create recorders

recorder Node -file output2/M3midNodeU3.out -time -node $flatmidNode -dof 3 disp
recorder Element -file output2/eleRTlayer1.out -time -ele 2226   material 1 fiber 1 TempAndElong;
recorder Element -file output2/eleRTlayer2.out -time -ele 2226   material 1 fiber 2 TempAndElong;
recorder Element -file output2/eleRTlayerWebS.out -time -ele 2226   material 1 fiber 3 TempAndElong;
recorder Element -file output2/eleRTlayerC.out -time -ele 2226   material 1 fiber 4 TempAndElong;
recorder Element -file output2/eleRTlayerWenS.out -time -ele 2226   material 1 fiber 5 TempAndElong;
recorder Element -file output2/eleRTlayerSimi.out -time -ele 2226   material 1 fiber 6 TempAndElong;
recorder Element -file output2/eleRTlayerRV.out -time -ele 2226   material 1 fiber 14 TempAndElong;
recorder Element -file output2/eleRTlayerRL.out -time -ele 2226   material 1 fiber 15 TempAndElong;
recorder Element -file output2/eleRTlayer21.out -time -ele 2226   material 1 fiber 21 TempAndElong;

recorder Element -file output2/eleFTlayer1.out -time -ele 2526   material 1 fiber 1 TempAndElong;
recorder Element -file output2/eleFTlayer2.out -time -ele 2526   material 1 fiber 2 TempAndElong;
recorder Element -file output2/eleFTlayerRV.out -time -ele 2526   material 1 fiber 4 TempAndElong;
recorder Element -file output2/eleFTlayerRL.out -time -ele 2526   material 1 fiber 5 TempAndElong;

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
set stepSize 0.5
set currentTime 0.0
set tFinal 1

constraints Transformation    		
numberer RCM			
system SparseSYM		
test NormDispIncr 1e-2 100 1	
algorithm KrylovNewton				
integrator LoadControl $stepSize				
analysis Static

while {$currentTime < $tFinal && $ok == 0} {
	puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	set ok [analyze 1];
	if {$ok == 0} {
		set currentTime [getTime]
		puts "Analysis step successful. Current time is $currentTime"
	} else {
		puts "Analysis failed at time = $currentTime"
		return
	}
}
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
	eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "M3ribtemp.dat"        -92.501    37.501
	} elseif {$sec=="flat"} { 
    #eleLoad -ele        $eleTag	     -type -shellThermal -source    "fileName"           $Y1     $Y2
	eleLoad  -ele [expr int($j*100+$i)]  -type -shellThermal -source "M3flattemp.dat"     [expr -$HalfT-0.001]    [expr $HalfT+0.001]
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
#set tFinal 10800
set tFinal 7200

constraints Transformation     		
numberer RCM			
system SparseSYM			
test NormDispIncr 1e-1 2000 1
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static

while {$currentTime < $tFinal && $ok == 0} {
	puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	set ok [analyze 1];
	if {$ok == 0} {
		set currentTime [getTime]
		puts "Analysis step successful. Current time is $currentTime"
	} else {
		puts "Analysis failed at time = $currentTime"
		
		return
	}
}