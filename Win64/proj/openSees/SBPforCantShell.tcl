#Korean Fire Tests on CHS
#Experimental Study on Limiting Temperatures of Circular Hollow Sections

# Units mm, N/mm²
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
set ANALYSIS "Has0Thermo"; 

wipe;
file mkdir WrapperData;
model BasicBuilder -ndm 2 -ndf 3;
#source DisplayPlane.tcl;		# procedure for displaying a plane in model
#source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model

#define node
set NumEles 24;
set BeamLen 3000.0;


set EleLen [expr $BeamLen/$NumEles]
for {set NodeID 0} {$NodeID <= $NumEles} {incr NodeID} {
set locX [expr $NodeID*$EleLen];
set NodeTag [expr $NodeID+1];
node $NodeTag $locX 0;
}

set EndNodeTag [expr $NumEles+1]
set MidNodeTag [expr $NumEles/2+1]
#define boundary condition;
#fix 1 1 1 1;
fix 1 1 1 1;
#fix $EndNodeTag 0 1 0;


#define an elastic material with Tag=1 and E=2e11.
#uniaxialMaterial ElasticThermal 2 19200 1.4e-5;
#uniaxialMaterial ElasticThermal 3 206000 0.5e-5;
uniaxialMaterial SteelECThermal 3 EC2NH 345 206000;
#uniaxialMaterial Steel02Thermal 3 280 200000 0.00001;
#define fibred section; Two fibres: fiber $yLoc $zLoc $A $matTag 
set fpc -30
set epsc0 -0.0025
set fpcu [expr $fpc*0.2];
set epsU -0.02
set lambda 0.1
set ft 3.0
set Ets [expr $ft/0.00838];

uniaxialMaterial ConcreteECThermal 2 $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets


set NumFibers 32;
set b 400.0;
set d 100.0;
set secTag 1;
set NumInts 5;
section fiberSecThermal $secTag {
	for {set fiberID 1} {$fiberID<=$NumFibers} {incr fiberID} {
		set fiberLocy [expr ($d/$NumFibers)*($fiberID-0.5)-$d/2];
		set fiberArea [expr ($d/$NumFibers)*$b*0.5];
		fiber $fiberLocy -10 $fiberArea 2;
		fiber $fiberLocy 10 $fiberArea 2;
		puts "$fiberLocy, $fiberArea";
		}
	fiber -30 -50 78.54 3;
	fiber -30 50 78.54 3;	
	fiber 30 -50 78.54 3;
	fiber 30 50 78.54 3;
	}


#define coordinate transforamtion: geomTransf $type $TransfTag;
#three transformation types can be chosen: Linear, PDelta, Corotational)
geomTransf Corotational 1 ; 
#geomTransf PDelta 1 ;

#define beam element: dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $TransfTag;
#"numIntgrPts" is the number of integration points along the element;
#"TransfTag" is pre-defined coordinate-transformation;
for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
	#element forceBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
}

#define output
set MidSpanNode [expr 1+$NumEles/2];
set MidSpanEle [expr $NumEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
puts "$FiberLocBot , $FiberLocTop";
recorder Node -file WrapperData/N_MidSpan_DOF.out -time -nodeRange 1 1 -dof 1 2 3 disp;
recorder Node -file WrapperData/N_NodeDeflect_DOF.out -time -node $EndNodeTag -dof 2 disp;
recorder Element -file WrapperData/N_ElementGF.out -time -ele 1   globalForce;
recorder Element -file WrapperData/N_ElementDef.out -time -ele1   basicDeformation;


recorder Element -file WrapperData/Element1Secm50T.out -time -eleRange 1  $NumEles section 1  fiber -50 10 TempElong;
recorder Element -file WrapperData/Element1Secm15T.out -time -eleRange 1  $NumEles section 1  fiber -25 0 TempElong;
recorder Element -file WrapperData/Element1Secm5T.out -time -eleRange 1  $NumEles section 1  fiber -5 10 TempElong;
recorder Element -file WrapperData/Element1Sec5T.out -time -eleRange 1  $NumEles section 1  fiber 5 10 TempElong;
recorder Element -file WrapperData/Element1Sec15T.out -time -eleRange 1  $NumEles section 1  fiber 25 10 TempElong;
recorder Element -file WrapperData/Element1Sec25T.out -time -eleRange 1  $NumEles section 1  fiber 50 10 TempElong;



recorder Element -file WrapperData/EleSecSSm50.out -time -eleRange 1  $NumEles section 1  fiber -50 10 stressStrainTangent;
recorder Element -file WrapperData/EleSecSSm25.out -time -eleRange 1  $NumEles section 1  fiber -25 10 stressStrainTangent;
recorder Element -file WrapperData/EleSecSSm5.out -time -eleRange 1  $NumEles section 1  fiber -5 10 stressStrainTangent;
recorder Element -file WrapperData/EleSecSS5.out -time -eleRange 1  $NumEles section 1  fiber 5 10 stressStrainTangent;
recorder Element -file WrapperData/EleSecSS25.out -time -eleRange 1  $NumEles section 1  fiber 25 10 stressStrainTangent;
recorder Element -file WrapperData/EleSecSS50.out -time -eleRange 1  $NumEles section 1  fiber 50 10 stressStrainTangent;

recorder Element -file WrapperData/SteelReinL.out -time -eleRange 1  1 section 1  fiber -30 -50 stressStrainTangent;
recorder Element -file WrapperData/SteelReinU.out -time -eleRange 1  1 section 1  fiber 30 -50 stressStrainTangent;

recorder Element -file WrapperData/EleForceSec1.out -time -eleRange 1  $NumEles section 1 Force;
recorder Element -file WrapperData/EleForceSec2.out -time -eleRange 1  $NumEles section 2  basic;
recorder Element -file WrapperData/EleForceSec3.out -time -eleRange 1  $NumEles section 3  forces;
recorder Element -file WrapperData/EleForceSec4.out -time -eleRange 1  $NumEles section 4  forces;
recorder Element -file WrapperData/EleForceSec5.out -time -eleRange 1  $NumEles section 5  forces;	


#display 2D deformation shape
set ViewScale 0.0001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel2D DeformedShape $ViewScale;

if {$MechANALYSIS == "HasPoint"} {
set Load 100;
#define Uniform load
puts "Now applying uniform load"

pattern Plain 1 Linear {
    set UDL -1;
	#for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
		#eleLoad -ele $eleID -type -beamUniform $UDL 0 0;
	#}

   load	$EndNodeTag 0  $Load 0;

};



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-4 100 2 ;
algorithm Newton;
integrator LoadControl 0.2;
analysis Static;
analyze 100;
}

if {$ANALYSIS == "HasThermo"} {
loadConst -time 0.0
puts "temperature loading"
# Define Thermal Load
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
load 1 -nodalThermal 800 $minusHalfD 400 $HalfD;

load $MidNodeTag -nodalThermal 200 $minusHalfD 100 $HalfD

load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD
#load $EndNodeTag -nodalThermal 800 $minusHalfD 400 $HalfD
#eleLoad -range 1 $NumEles -type -beamThermal -source -node
#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $MidNodeTag 0.5 $EndNodeTag 1; 

eleLoad -range 1 $NumEles -type -beamThermal 1000 -$HalfD 0 $HalfD

}


constraints Plain;
numberer Plain;
system BandGeneral;
test NormUnbalance 1.0e-4 50 4;
algorithm Newton;
integrator LoadControl 0.01;
analysis Static;
analyze 100;

}
wipe;