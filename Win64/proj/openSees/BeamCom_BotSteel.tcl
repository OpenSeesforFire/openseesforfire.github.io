
# define model parameters

# general
set l 4.15  
set w 0.5
set t 0.13
set ct 0.075
set dt 0.00075
set rt 0.055
set lw 0.13;#the lower width of the rib
set uw 0.182;#the upper width of the rib
set fw [expr 0.5*($w-$uw)]
set ro 0.000198
set UDL -5430
set HalfT [expr 0.5*$ct]

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
set fpcu [expr 0.05*$fpc] 
set epsU [expr 10*$epsc0]
set ft [expr 0.1*$fc]
set epstu 0.002
set Ec [expr 1.5*abs($fpc/$epsc0)]
set Ets [expr $ft/$epstu]
set v 0.2
# set J [expr $rt*pow($uw,3)/3.0]; #torsional moment
# set G [expr $Ec/(2.0*(1+$v))]; #shear modulus 
# set JG [expr $J*$G] 
set gt [expr $ft/$Ec*$ft*4]
set gc [expr $fc/$Ec*$fc*6]
set lambda 0.15


####################################################################################

# Start up OpenSees model builder and create output folder

#wipe
model BasicBuilder -ndm 2 -ndf 3
file mkdir Tbeam

# Load display engine
# source DisplayPlane.tcl
# source DisplayModel2D.tcl
# source DisplayModel3D.tcl

####################################################################################

# Create nodes

set NumEles 50;
set BeamLen 4.15;

set EleLen [expr $BeamLen/$NumEles]

for {set NodeID 0} {$NodeID <= $NumEles} {incr NodeID} {
set locX [expr $NodeID*$EleLen];
set NodeTag [expr $NodeID+1];
      node $NodeTag $locX 0;
	  #puts "$NodeTag + $locX"
}
# set any node tags you want to use later here:
set flatmidNode [expr int($NumEles*0.5+1)];#the mid node in the thick part
#puts "flatmidNode = $flatmidNode"

#return

####################################################################################

# Create constraints
fix 1 1 1 0;
fix 51 0 1 0;
 
####################################################################################

# Create geometric transformation
geomTransf Corotational    1   

####################################################################################

# Create material models

# Reinforcement
 #uniaxialMaterial SteelECThermal 1  EC2NC  $fy  $Es
 uniaxialMaterial SteelECThermal 1 EC2NC $rfy $rEs;
  #uniaxialMaterial Steel01Thermal 1      $rfy $rEs 0.015;
 #nDsection PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
  # nDsection PlateRebarThermal 2 	 1 				0
  # nDsection PlateRebarThermal 3 	 1 				90
 
# Flat part Concrete 
 #nDsection  CDPPlaneStressThermal matTag Ec  v  ft  fc  gt  gc
  #nDsection  CDPPlaneStressThermal 4 	  $Ec $v $ft $fc $gt $gc
 #nDsection   PlateFromPlaneStressThermal    matTag   	nDMatTag    forStability
  #nDsection   PlateFromPlaneStressThermal    5   		4   		1e9
  uniaxialMaterial ConcreteECThermal 11 $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets

  
# Steel deck 
 #uniaxialMaterial Steel01Thermal $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
  uniaxialMaterial Steel01Thermal 6      $dfy $dEs 0.015;   # Steel deck in the ribs
 #nDsection PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
  #nDsection PlateRebarThermal 7 	  6 				0
  #nDsection PlateRebarThermal 8 	  6 				90
  
 #nDsection  J2PlaneStressThermal 9 21 $dEs $vs $dfy $dfu 0.015
 #nDsection   PlateFromPlaneStressThermal    20   9   2e10;
 
####################################################################################

# Create sections
set secTag 3
# section fiberSecThermal $secTag {
	# #fiber $locy $locx $area $matTag
	# fiber  -0.092125 0  0.0977625 6
	# fiber  -0.08675 0  1.3544  11
    # fiber  -0.07675 0  1.4489  11
	# fiber  -0.06675 0  1.5435  11
	# fiber  -0.05675 0  1.638  11
	# fiber  -0.04675 0  1.7325  11
	# fiber  -0.039625 0  0.7649575  11
	# fiber  -0.032125 0  978.25  11
	# fiber  -0.02175 0  910  11
	# fiber  -0.016651 0  18.018  1
	# fiber  -0.016453 0  18.018  1
	# fiber  -0.011552 0  873.964  11
    # fiber  -0.00175 0  910  11
	# fiber   0.00825 0  910  11
	# fiber   0.01825 0  910  11
	# fiber   0.02825 0  910  11
	# fiber   0.035375 0  386.75  11
	# }
	
section fiberSecThermal $secTag {
	#fiber $locy $locx $area $matTag
	fiber  -0.064625 0  0.0977625e-3 6
	fiber  -0.05925 0  1.3544e-3  11
    fiber  -0.04925 0  1.4489e-3  11
	fiber  -0.03925 0  1.5435e-3  11
	fiber  -0.02925 0  1.638e-3  11
	fiber  -0.01925 0  1.7325e-3  11
	fiber  -0.012125 0  0.7649575e-3  11
	fiber  -0.004625 0  5.375e-3  11
	fiber   0.00575 0  0.005  11
	fiber   0.010849 0  0.099e-3  1
	fiber   0.011047 0  0.099e-3  1
	fiber   0.015948 0  4.802e-3  11
    fiber   0.02575 0  0.005  11
	fiber   0.03575 0  0.005  11
	fiber   0.04575 0  0.005  11
	fiber   0.05575 0  0.005  11
	fiber   0.062875 0  2.125e-3  11
	}

####################################################################################

# Create elements

set NumInts 5;
for {set eleID 1} {$eleID <= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag 1;
	#element forceBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts $secTag 1;
	#puts "$eleID + $NodeTag0 + $NodeTag1 + $NumInts + $secTag"
}
#return

####################################################################################

# Create recorders

recorder Node -file Tbeam/M3midNodeU3.out -time -node $flatmidNode -dof 2 disp;

set Rmidtag 26
recorder Element -file Tbeam/RlayerS1.out -time -ele $Rmidtag   section 1 fiber 1 TempElong;
# recorder Element -file Tbeam/RStressS1.out -time	-ele $Rmidtag	section 1	fiber 1	stress;
# recorder Element -file Tbeam/RStrainS1.out -time	-ele $Rmidtag	section  1	fiber 1	strain;
# recorder Element -file Tbeam/RlayerC1.out -time -ele $Rmidtag   section 1 fiber 2 TempAndElong;
# recorder Element -file Tbeam/RStressC1.out -time	-ele $Rmidtag	section  1	fiber 2	stress;
# recorder Element -file Tbeam/RstrainC1.out -time	-ele $Rmidtag	section  1	fiber 2	strain;
# # recorder Element -file Tbeam/RlayerS2.out -time -ele $Rmidtag   section 1 fiber 3 TempAndElong;
# # recorder Element -file Tbeam/RStressS2.out -time -ele $Rmidtag   section 1 fiber 3 stress;
# # recorder Element -file Tbeam/RstrainS2.out -time -ele $Rmidtag   section 1 fiber 3 strain;
# recorder Element -file Tbeam/RlayerC2.out -time -ele $Rmidtag   section 1 fiber 3 TempAndElong;
# recorder Element -file Tbeam/RStressC2.out -time -ele $Rmidtag   section 1 fiber 3 stress;
# recorder Element -file Tbeam/RstrainC2.out -time -ele $Rmidtag   section 1 fiber 3 strain;
# # recorder Element -file Tbeam/RlayerS3.out -time -ele $Rmidtag  section 1 fiber 5 TempAndElong;
# # recorder Element -file Tbeam/RStressS3.out -time -ele $Rmidtag   section 1 fiber 5 stress;
# # recorder Element -file Tbeam/RStrainS3.out -time -ele $Rmidtag   section 1 fiber 5 strain;
# recorder Element -file Tbeam/RlayerC3.out -time -ele $Rmidtag   section 1 fiber 4 TempAndElong;
# recorder Element -file Tbeam/RStressC3.out -time -ele $Rmidtag   section 1 fiber 4 stress;
# recorder Element -file Tbeam/RstrainC3.out -time -ele $Rmidtag   section 1 fiber 4 strain;
# # recorder Element -file Tbeam/RlayerS4.out -time -ele $Rmidtag  section 1 fiber 7 TempAndElong;
# # recorder Element -file Tbeam/RStressS4.out -time -ele $Rmidtag   section 1 fiber 7 stress;
# # recorder Element -file Tbeam/RStrainS4.out -time -ele $Rmidtag   section 1 fiber 7 strain;
# recorder Element -file Tbeam/RlayerC4.out -time -ele $Rmidtag   section 1 fiber 5 TempAndElong;
# recorder Element -file Tbeam/RStressC4.out -time -ele $Rmidtag   section 1 fiber 5 stress;
# recorder Element -file Tbeam/RstrainC4.out -time -ele $Rmidtag   section 1 fiber 5 strain;
# # recorder Element -file Tbeam/RlayerS5.out -time -ele $Rmidtag   section 1 fiber 9 TempAndElong;
# # recorder Element -file Tbeam/RStressS5.out -time -ele $Rmidtag   section 1 fiber 9 stress;
# # recorder Element -file Tbeam/RStrainS5.out -time -ele $Rmidtag  section 1 fiber 9 strain;
# recorder Element -file Tbeam/RlayerC5.out -time -ele $Rmidtag   section 1 fiber 6 TempAndElong;
# recorder Element -file Tbeam/RStressC5.out -time -ele $Rmidtag   section 1 fiber 6 stress;
# recorder Element -file Tbeam/RstrainC5.out -time -ele $Rmidtag   section 1 fiber 6 strain;
# recorder Element -file Tbeam/RlayerC6.out -time -ele $Rmidtag   section 1 fiber 7 TempAndElong;
# recorder Element -file Tbeam/RStressC6.out -time -ele $Rmidtag   section 1 fiber 7 stress;
# recorder Element -file Tbeam/RstrainC6.out -time -ele $Rmidtag   section 1 fiber 7 strain;
# recorder Element -file Tbeam/RlayerC7.out -time -ele $Rmidtag  section 1 fiber 8 TempAndElong;
# recorder Element -file Tbeam/RStressC7.out -time	-ele $Rmidtag	section  1	fiber 8	stress;
# recorder Element -file Tbeam/RstrainC7.out -time	-ele $Rmidtag	section  1	fiber 8	strain;
# recorder Element -file Tbeam/RlayerC8.out -time -ele $Rmidtag   section 1 fiber 9 TempAndElong;
# recorder Element -file Tbeam/RStressC8.out -time	-ele $Rmidtag	section  1	fiber 9	stress;
# recorder Element -file Tbeam/RstrainC8.out -time	-ele $Rmidtag	section  1	fiber 9	strain;
# # recorder Element -file Tbeam/RlayerRV.out -time -ele $Rmidtag   section 1 fiber 14 TempAndElong;
# # recorder Element -file Tbeam/RStressRV.out -time -ele $Rmidtag   section 1 fiber 14 stress;
# # recorder Element -file Tbeam/RStrainRV.out -time -ele $Rmidtag   section 1 fiber 14 strain;
# # recorder Element -file Tbeam/RlayerRL.out -time -ele $Rmidtag   section 1 fiber 15 TempAndElong;
# # recorder Element -file Tbeam/RStressRL.out -time -ele $Rmidtag   section 1 fiber 15 stress;
# # recorder Element -file Tbeam/RStrainRL.out -time -ele $Rmidtag  section 1 fiber 15 strain;
# recorder Element -file Tbeam/RlayerC9.out -time -ele $Rmidtag   section 1 fiber 12 TempAndElong;
# recorder Element -file Tbeam/RStressC9.out -time	-ele $Rmidtag	section  1	fiber 12	stress;
# recorder Element -file Tbeam/RstrainC9.out -time	-ele $Rmidtag	section  1	fiber 12	strain;
# recorder Element -file Tbeam/RlayerC10.out -time -ele $Rmidtag   section 1 fiber 13 TempAndElong;
# recorder Element -file Tbeam/RStressC10.out -time	-ele $Rmidtag	section  1	fiber 13	stress;
# recorder Element -file Tbeam/RstrainC10.out -time	-ele $Rmidtag	section  1	fiber 13	strain;
# recorder Element -file Tbeam/RlayerC11.out -time -ele $Rmidtag   section 1 fiber 14 TempAndElong;
# recorder Element -file Tbeam/RStressC11.out -time	-ele $Rmidtag	section  1	fiber 14	stress;
# recorder Element -file Tbeam/RstrainC11.out -time	-ele $Rmidtag	section  1	fiber 14	strain;
# recorder Element -file Tbeam/RlayerC12.out -time -ele $Rmidtag   section 1 fiber 15 TempAndElong;
# recorder Element -file Tbeam/RStressC12.out -time	-ele $Rmidtag	section  1	fiber 15	stress;
# recorder Element -file Tbeam/RstrainC12.out -time	-ele $Rmidtag	section  1	fiber 15	strain;
# recorder Element -file Tbeam/RlayerC13.out -time -ele $Rmidtag   section 1 fiber 16 TempAndElong;
# recorder Element -file Tbeam/RStressC13.out -time	-ele $Rmidtag	section  1	fiber 16	stress;
# recorder Element -file Tbeam/RstrainC13.out -time	-ele $Rmidtag	section  1	fiber 16	strain;
# recorder Element -file Tbeam/RlayerC14.out -time -ele $Rmidtag  section 1 fiber 17 TempAndElong;
# recorder Element -file Tbeam/RStressC14.out -time	-ele $Rmidtag	section  1	fiber 17	stress;
# recorder Element -file Tbeam/RstrainC14.out -time	-ele $Rmidtag	section  1	fiber 17	strain;



####################################################################################

# Create ambient loads
set f3 [expr $UDL*$w*$EleLen]
pattern Plain 1 Linear {
	for {set i 0} {$i <= $NumEles} {incr i 1} {
	   set nodeID [expr $i+1] 
	   #load nodeTag	f1 	f2  f3 	
	    load  $nodeID 	0   $f3  0; 
	}
}


####################################################################################

# Analyse ambient loads

set ok 0
set stepSize 0.1
set currentTime 0.0
set tFinal 1

constraints Transformation    		
numberer Plain			
system SparseSYM		
test NormDispIncr 1e-4 100 1	
algorithm Newton				
integrator LoadControl $stepSize				
analysis Static
analyze 10;

loadConst -time 0.0

####################################################################################

# Create thermal loads


pattern Plain 2 Linear {

   for {set i 1} {$i <= $NumEles} {incr i 1} {
	#eleLoad   -ele   $eleTag	 -type -shellThermal -source    "fileName"           $Y1     $Y2   
	eleLoad    -ele    $i       -type -beamThermal  500   -0.065    100  0.065
    }	
}

####################################################################################

# Analyse Thermal loads

set ok 0
set stepSize 0.01
set currentTime 0.0
set tFinal 1
#set tFinal 360

constraints Transformation     		
numberer RCM			
system SparseSYM			
test NormDispIncr 1e-2 1000 1
#algorithm KrylovNewton		
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
wipe;

wipe;