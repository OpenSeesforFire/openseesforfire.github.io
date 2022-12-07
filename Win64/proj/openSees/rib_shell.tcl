# Units m, N/m
set MechANALYSIS "HasPoint";
set ANALYSIS "HasThermo"; 

wipe;
file mkdir rib_shell;
model BasicBuilder -ndm 3 -ndf 6;


#define geometry
set l 3.0;
set w 0.4;
set d 0.05;
set nx 10;
set ny 2;

#define material
set fc 30e6
set fpc -30e6
set epsc0 -0.0025
set fpcu [expr $fpc*0.2];
set epsU -0.02
set Ec [expr 2*abs($fpc/$epsc0)]
set lambda 0.15
set ft [expr 0.1*$fc]
set Ets [expr $ft/0.0010056];
set gt [expr $ft/$Ec*$ft*4]
set gc [expr $fc/$Ec*$fc*6]
set v 0.2

#uniaxialMaterial ConcreteECThermal 2 $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets

# elastic
nDMaterial ElasticIsotropic3DThermal 2 $Ec $v 0 1.0e-5;
nDMaterial PlateFiberThermal 4 2;

#plastic
set gt [expr 3.7e6/2.2e10*3.7e6*2];
set gc [expr 37.0e6/2.2e10*37.0e6*6];
#nDMaterial  CDPPlaneStressThermal 100 2.2e10 0.2 3.7e6 37e6 $gt $gc;
#nDMaterial   PlateFromPlaneStressThermal    4   100   1e9;

#nDMaterial  CDPPlaneStressThermal matTag Ec  v  ft  fc  gt  gc
#nDMaterial  CDPPlaneStressThermal 100  $Ec $v  $ft  $fc $gt $gc;
#nDMaterial   PlateFromPlaneStressThermal    4   100  10e9;

# define node
set elemx [expr $l/$nx]
set elemy [expr $w/$ny]
for {set i 0} {$i <= $nx} {incr i} {
    set x [expr $i*$elemx]
    for {set j 0} {$j <= $ny} {incr j} {
       set y [expr $j*$elemy]
       set NodeTag [expr ($j+1)*100+($i+1)];
       node $NodeTag $x $y 0;
	   #puts "$NodeTag + $x + $y "
    }
}	
#return

fixX 0 1 1 1 1 1 1;

set sec 4
section  LayeredShellThermal   $sec	   10  	4 0.005    4 0.005   4 0.005  4 0.005	 4 0.005   4 0.005   4 0.005   4 0.005	  4 0.005	4  0.005;

set offset 0.0
#section  LayeredShellThermal   $sec	 -offset $offset  10  	4 0.005    4 0.005   4 0.005  4 0.005	 4 0.005   4 0.005   4 0.005   4 0.005	  4 0.005	4  0.005;


for {set i 1} {$i <= $nx} {incr i} {
    for {set j 1} {$j <= $ny} {incr j} {
       set node1 [expr $j*100+$i]
	   set node2 [expr $node1+1]
	   set node3 [expr ($j+1)*100+$i+1]
	   set node4 [expr $node3-1]
       set eleID [expr $j*100+$i]
       element ShellNLDKGQThermal $eleID $node1 $node2 $node3 $node4 $sec
	   #puts "$eleID+ $node1 + $node2 + $node3 + $node4 + $sec "
    }
}
#return
#define output
set EndNodeTag2 [expr ($ny/2+1)*100+$nx+1]
set MidEle [expr ($ny/2+1)*100+($nx/2+1)] 

recorder Node -file rib_shell/EndnodeU3.out -time -node $EndNodeTag2 -dof 3 disp;
recorder Element -file rib_shell/REleForce.out -time -ele $MidEle force;

recorder Element -file rib_shell/Rlayer1.out -time -ele $MidEle   material 1 fiber 1 TempAndElong;
recorder Element -file rib_shell/Rstrain1.out -time -ele $MidEle  material 1 fiber 1 strain;
recorder Element -file rib_shell/Rstress.out -time -ele $MidEle  material 1 fiber 1 stress;

recorder Element -file rib_shell/Rlayer2.out -time -ele $MidEle   material 1 fiber 2 TempAndElong;
recorder Element -file rib_shell/Rstrain2.out -time -ele $MidEle  material 1 fiber 2 strain;
recorder Element -file rib_shell/Rlayer3.out -time -ele $MidEle   material 1 fiber 3 TempAndElong;
recorder Element -file rib_shell/Rstrain3.out -time -ele $MidEle  material 1 fiber 3 strain;
recorder Element -file rib_shell/Rlayer4.out -time -ele $MidEle   material 1 fiber 4 TempAndElong;
recorder Element -file rib_shell/Rstrain4.out -time -ele $MidEle  material 1 fiber 4 strain;
recorder Element -file rib_shell/Rlayer5.out -time -ele $MidEle   material 1 fiber 5 TempAndElong;
recorder Element -file rib_shell/Rstrain5.out -time -ele $MidEle  material 1 fiber 5 strain;
recorder Element -file rib_shell/Rlayer6.out -time -ele $MidEle   material 1 fiber 6 TempAndElong;
recorder Element -file rib_shell/Rstrain6.out -time -ele $MidEle  material 1 fiber 6 strain;
recorder Element -file rib_shell/Rlayer7.out -time -ele $MidEle   material 1 fiber 7 TempAndElong;
recorder Element -file rib_shell/Rstrain7.out -time -ele $MidEle  material 1 fiber 7 strain;
recorder Element -file rib_shell/Rlayer8.out -time -ele $MidEle   material 1 fiber 8 TempAndElong;
recorder Element -file rib_shell/Rstrain8.out -time -ele $MidEle  material 1 fiber 8 strain;
recorder Element -file rib_shell/Rlayer9.out -time -ele $MidEle   material 1 fiber 9 TempAndElong;
recorder Element -file rib_shell/Rstrain9.out -time -ele $MidEle  material 1 fiber 9 strain;
recorder Element -file rib_shell/Rlayer10.out -time -ele $MidEle   material 1 fiber 10 TempAndElong;
recorder Element -file rib_shell/Rstrain10.out -time -ele $MidEle  material 1 fiber 10 strain;
recorder Element -file rib_shell/Rstress10.out -time -ele $MidEle  material 1 fiber 10 stress;
#return

if {$MechANALYSIS == "HasPoint"} {

pattern Plain 1 Linear {
    set Load 100.0;

	# for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
		# eleLoad -ele $eleID -type -beamUniform $UDL 0 0;
	# }
   for {set ID 0} {$ID<=$ny} {incr ID} {
		set nodeID [expr ($ID+1)*100+($nx+1)];
		load $nodeID  0 0 [expr $Load/($ny+1)] 0 0 0 ;
    }
	

};

constraints Transformation;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-3 50 1 ;
algorithm KrylovNewton;
integrator LoadControl 0.1 
analysis Static
analyze 10
loadConst -time 0.0
}

if {$ANALYSIS == "HasThermo"} {

set HalfD [expr $d/2]

pattern Plain 2 Linear {

for {set i 1} {$i <= $nx} {incr i} {
    for {set j 1} {$j <= $ny} {incr j} {
	set eleID [expr $j*100+$i];
    #eleLoad -ele $eleID -type -shellThermal 500 -$HalfD  100 $HalfD
    eleLoad -ele $eleID -type -shellThermal 500 [expr -$HalfD-$offset] 100 [expr $HalfD-$offset]
    }
  }
}
set stepSize 0.01
set tFinal 1

constraints Transformation     		
numberer RCM			
system BandGeneral			
test NormDispIncr 1e-3 1000 1
algorithm KrylovNewton		
#algorithm Newton			
integrator LoadControl $stepSize				
analysis Static
analyze 100;

}
wipe;