#Korean Fire Tests on CHS
#Experimental Study on Limiting Temperatures of Circular Hollow Sections

# Units m, N/m
#		o-------------------------------------------------------o <---
#		1		2		3		4	5	6		7		8		9
#					y
#					|
#					|
#					|_________	x
#				   /
#				  /
#				 /z
set MechANALYSIS "HasPoint";
set ANALYSIS "HasThermo";

wipe;
file mkdir rib_beam;
model BasicBuilder -ndm 2 -ndf 3;


#define node
set NumEles 10;
set BeamLen 3.0;


set EleLen [expr $BeamLen/$NumEles]
#puts "$EleLen"
for {set NodeID 0} {$NodeID <= $NumEles} {incr NodeID} {
    set locX [expr $NodeID*$EleLen];
    set NodeTag [expr $NodeID+1];
    node $NodeTag $locX 0;
	#puts "$NodeTag + $locX "
}
#return
set EndNodeTag [expr $NumEles+1]
set MidNodeTag [expr $NumEles/2+1]

fix 1 1 1 1;
#fix 1 1 1 0;
#fix $EndNodeTag 0 1 0;


#define an elastic material with Tag=1 and E=2e11.
#uniaxialMaterial ElasticThermal 2 17900 1.4e-5;
#uniaxialMaterial ElasticThermal 3 206000 0.5e-5;
#uniaxialMaterial SteelECThermal 3 EC2NH 345 206000;
#uniaxialMaterial Steel02Thermal 3 280 200000 0.00001;
#define fibred section; Two fibres: fiber $yLoc $zLoc $A $matTag
set fpc -30e6
set epsc0 -0.0025
set fpcu [expr $fpc*0.2];
set epsU -0.02
set lambda 0.15
set ft 3.0e6
set Ets [expr $ft/0.0010056];

uniaxialMaterial ConcreteECThermal 2 $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets


set NumFibers 10;
set b 0.4;
set d 0.05;
set secTag 1;
set NumInts 5;
section fiberSecThermal $secTag {
	for {set fiberID 1} {$fiberID<=$NumFibers} {incr fiberID} {
		set fiberLocy [expr ($d/$NumFibers)*($fiberID-0.5)-$d/2];
		set fiberArea [expr ($d/$NumFibers)*$b*0.5];
		fiber $fiberLocy [expr -$b/4] $fiberArea 2;
		fiber $fiberLocy [expr $b/4] $fiberArea 2;
		puts "$fiberLocy, $fiberArea";
		}
	#fiber -25 -50 56.5 3;
	#fiber -25 50 56.5 3;
	#fiber 30 -50 78.54 3;
	#fiber 30 50 78.54 3;
	}


#define coordinate transforamtion: geomTransf $type $TransfTag;
#three transformation types can be chosen: Linear, PDelta, Corotational)
geomTransf Corotational 1 ;
#geomTransf PDelta 1 ;

for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	#element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
}

#define output
set MidEle [expr $NumEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
#puts "$FiberLocBot , $FiberLocTop";
#recorder Node -file WrapperData/MidSpanU2.out -time -node $MidNodeTag -dof 2 disp;
recorder Node -file rib_beam/EndnodeU2.out -time -node $EndNodeTag -dof 2 disp;

recorder Element -file  rib_beam/Rblayer1.out -time -ele $MidEle   section 1 fiber -0.0225 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain1.out -time -ele $MidEle  section 1 fiber -0.0225 0.1 strain;
recorder Element -file  rib_beam/Rblayer2.out -time -ele $MidEle   section 1 fiber -0.0175 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain2.out -time -ele $MidEle  section 1 fiber -0.0175 0.1 strain;
recorder Element -file  rib_beam/Rblayer3.out -time -ele $MidEle   section 1 fiber -0.0125 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain3.out -time -ele $MidEle  section 1 fiber -0.0125 0.1 strain;
recorder Element -file  rib_beam/Rblayer4.out -time -ele $MidEle   section 1 fiber -0.0075 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain4.out -time -ele $MidEle  section 1 fiber -0.0075 0.1 strain;
recorder Element -file  rib_beam/Rblayer5.out -time -ele $MidEle   section 1 fiber -0.0025 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain5.out -time -ele $MidEle  section 1 fiber -0.0025 0.1 strain;
recorder Element -file  rib_beam/Rblayer6.out -time -ele $MidEle   section 1 fiber 0.0025 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain6.out -time -ele $MidEle  section 1 fiber 0.0025 0.1 strain;
recorder Element -file  rib_beam/Rblayer7.out -time -ele $MidEle   section 1 fiber 0.0075 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain7.out -time -ele $MidEle  section 1 fiber 0.0075 0.1 strain;
recorder Element -file  rib_beam/Rblayer8.out -time -ele $MidEle   section 1 fiber 0.0125 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain8.out -time -ele $MidEle  section 1 fiber 0.0125 0.1 strain;
recorder Element -file  rib_beam/Rblayer9.out -time -ele $MidEle   section 1 fiber 0.0175 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain9.out -time -ele $MidEle  section 1 fiber 0.0175 0.1 strain;
recorder Element -file  rib_beam/Rblayer10.out -time -ele $MidEle  section 1 fiber 0.0225 0.1 TempElong;
recorder Element -file  rib_beam/Rbstrain10.out -time -ele $MidEle section 1 fiber 0.0225 0.1 strain;


# recorder Element -file WrapperData/Element1Secm50T.out -time -eleRange 1  $MidSpanEle section 1  fiber -50 10 TempElong;
# recorder Element -file WrapperData/Element1Secm15T.out -time -eleRange 1  $MidSpanEle section 1  fiber -42.5 0 TempElong;
# recorder Element -file WrapperData/Element1Secm5T.out -time -eleRange 1  $MidSpanEle section 1  fiber -37.5 10 TempElong;
# recorder Element -file WrapperData/Element1Sec5T.out -time -eleRange 1  $MidSpanEle section 1  fiber 37.5 10 TempElong;
# recorder Element -file WrapperData/Element1Sec15T.out -time -eleRange 1  $MidSpanEle section 1  fiber 42.5 10 TempElong;
# recorder Element -file WrapperData/Element1Sec25T.out -time -eleRange 1  $MidSpanEle section 1  fiber 50 10 TempElong;



# recorder Element -file WrapperData/EleSecSSm50.out -time -eleRange 1  $MidSpanEle section 1  fiber -50 10 stressStrainTangent;
# recorder Element -file WrapperData/EleSecSSm25.out -time -eleRange 1  $MidSpanEle section 1  fiber -42.5 10 stressStrainTangent;
# recorder Element -file WrapperData/EleSecSSm5.out -time -eleRange 1  $MidSpanEle section 1  fiber -5 10 stressStrainTangent;
# recorder Element -file WrapperData/EleSecSS5.out -time -eleRange 1  $MidSpanEle section 1  fiber 5 10 stressStrainTangent;
# recorder Element -file WrapperData/EleSecSS25.out -time -eleRange 1  $MidSpanEle section 1  fiber 42.5 10 stressStrainTangent;
# recorder Element -file WrapperData/EleSecSS50.out -time -eleRange 1  $MidSpanEle section 1  fiber 50 10 stressStrainTangent;

#recorder Element -file WrapperData/SteelReinL.out -time -ele  $MidSpanEle section 1  fiber -25 -50 stressStrainTangent;
#recorder Element -file WrapperData/SteelReinU.out -time -eleRange 1  1 section 1  fiber 30 -50 stressStrainTangent;

# recorder Element -file WrapperData/EleForceSec1.out -time -eleRange 1  $MidSpanEle section 1 Force;
# recorder Element -file WrapperData/EleForceSec2.out -time -eleRange 1  $MidSpanEle section 2  basic;
# recorder Element -file WrapperData/EleForceSec3.out -time -eleRange 1  $MidSpanEle section 3  forces;
# recorder Element -file WrapperData/EleForceSec4.out -time -eleRange 1  $MidSpanEle section 4  forces;
# recorder Element -file WrapperData/EleForceSec5.out -time -eleRange 1  $MidSpanEle section 5  forces;


#display 2D deformation shape
set ViewScale 0.00001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel2D DeformedShape $ViewScale;

if {$MechANALYSIS == "HasPoint"} {

pattern Plain 1 Linear {
    set Load 100;
	# for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
		# eleLoad -ele $eleID -type -beamUniform $UDL 0 0;
	# }
	load  $EndNodeTag 0  $Load 0;

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

set minusHalfD [expr -$d/2]
set HalfD [expr $d/2]

pattern Plain 2 Linear {

#load 1 -nodalThermal 700 $minusHalfD 400 $HalfD;
#load 2 -nodalThermal 525 $minusHalfD 290 $HalfD;
#load 3 -nodalThermal 380 $minusHalfD 200 $HalfD
#load 4 -nodalThermal 265 $minusHalfD 130 $HalfD
#load 5 -nodalThermal 180 $minusHalfD 80 $HalfD
#load 6 -nodalThermal 125 $minusHalfD 50 $HalfD
#load 7 -nodalThermal 100 $minusHalfD 40 $HalfD
# load 1 -nodalThermal 800 $minusHalfD 400 $HalfD;

# load $MidNodeTag -nodalThermal 200 $minusHalfD 100 $HalfD

# load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD
#load $EndNodeTag -nodalThermal 800 $minusHalfD 400 $HalfD
#eleLoad -range 1 $NumEles -type -beamThermal -source -node
#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $MidNodeTag 0.5 $EndNodeTag 1;

eleLoad -range 1 $NumEles -type -beamThermal 500 -$HalfD 100 $HalfD

}

set stepSize 0.01
set tFinal 1

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1e-2 1000 1
algorithm KrylovNewton
#algorithm Newton
integrator LoadControl $stepSize
analysis Static
analyze 100;
}
wipe;