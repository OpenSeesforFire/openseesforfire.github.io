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

set ANALYSIS "HasLoad"; 
set TANALYSIS "Has0Thermo"

wipe;
file mkdir WrapperData;
model BasicBuilder -ndm 3 -ndf 6;
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl;		# procedure for displaying 2D perspective of model

#define node
set NumEles 2;
set BeamLen 3000.0;
set EleLen [expr $BeamLen/$NumEles]
for {set NodeID 0} {$NodeID <= $NumEles} {incr NodeID} {
set locX [expr $NodeID*$EleLen];
set NodeTag [expr $NodeID+1];
node $NodeTag $locX 0 0;
}

set EndNodeTag [expr $NumEles+1]
set MidNodeTag [expr $NumEles/2+1]
#define boundary condition;
fix 1 1 1 1 1 1 1 ; 
#fix $EndNodeTag 1 1 1 1 1 1 ;


#define an elastic material with Tag=1 and E=2e11.
uniaxialMaterial ElasticThermal 1 200000 1.2e-5;
uniaxialMaterial StainlessECThermal 2  Grade14571 240 2e5  520;
uniaxialMaterial SteelECThermal 3 300 200000;
#uniaxialMaterial Steel01Thermal 1 300 200000 0.00001;
#define fibred section; Two fibres: fiber $yLoc $zLoc $A $matTag 

set NumFibers 4;
set b 40.0;
set d 80.0;
set secTag 1;
set NumInts 5;
set GJ 3e10;
section fiberSecThermal $secTag -GJ $GJ {
	for {set fiberID 1} {$fiberID<=$NumFibers} {incr fiberID} {
		set fiberLocy [expr ($d/$NumFibers)*($fiberID-0.5)-$d/2];
		set fiberArea [expr ($d/$NumFibers)*$b*0.5];
		fiber $fiberLocy -10 $fiberArea 2;
		fiber $fiberLocy 10 $fiberArea 2;
		puts "$fiberLocy, $fiberArea";
		}
	}


#define coordinate transforamtion: geomTransf $type $TransfTag;
#three transformation types can be chosen: Linear, PDelta, Corotational)
geomTransf Corotational 1 0 0 1 ;
geomTransf Linear 2 0 0 1 ;
geomTransf PDelta 3 0 0 1 ;
set GeomTransf 2;

#define beam element: dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $TransfTag;
#"numIntgrPts" is the number of integration points along the element;
#"TransfTag" is pre-defined coordinate-transformation;
for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 $GeomTransf;
}

#define output
set MidSpanNode [expr 1+$NumEles/2];
set MidSpanEle [expr $NumEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
puts "$FiberLocBot , $FiberLocTop";
recorder Node -file WrapperData/N_EndNode_DOF.out -time -node $EndNodeTag -dof 1 2 3 disp;
recorder Node -file WrapperData/N_NodeDeflect_DOF.out -time -nodeRange 1 $EndNodeTag -dof 2 disp;
recorder Element -file WrapperData/N_ElementGF.out -time -ele 1   globalForce;
recorder Element -file WrapperData/N_ElementDef.out -time -ele 1   basicDeformation;

if {$TANALYSIS == "HasThermo"} {
recorder Element -file WrapperData/Element1Secm35.out -time -eleRange 1  $NumEles section 1  fiber -35 -5 TempElong;
recorder Element -file WrapperData/Element1Secm25.out -time -eleRange 1  $NumEles section 1  fiber -25 0 TempElong;
recorder Element -file WrapperData/Element1Secm15.out -time -eleRange 1  $NumEles section 1  fiber -15 0 TempElong;
recorder Element -file WrapperData/Element1Secm5.out -time -eleRange 1  $NumEles section 1  fiber -5 0 TempElong;
recorder Element -file WrapperData/Element1Secp5.out -time -eleRange 1  $NumEles section 1  fiber 5 0 TempElong;
recorder Element -file WrapperData/Element1Secp15.out -time -eleRange 1  $NumEles section 1  fiber 15 0 TempElong;
recorder Element -file WrapperData/Element1Secp25.out -time -eleRange 1  $NumEles section 1  fiber 25 0 TempElong;
recorder Element -file WrapperData/Element1Secp35.out -time -eleRange 1  $NumEles section 1  fiber 35 0 TempElong;

recorder Element -file WrapperData/Element1Sec2.out -time -eleRange 1  $NumEles section 2  fiber -35 0 TempElong;
recorder Element -file WrapperData/Element1Sec3.out -time -eleRange 1  $NumEles section 3  fiber -35 0 TempElong;
recorder Element -file WrapperData/Element1Sec4.out -time -eleRange 1  $NumEles section 4  fiber -35 0 TempElong;
recorder Element -file WrapperData/Element1Sec5.out -time -eleRange 1  $NumEles section 5  fiber -35 0 TempElong;
}

recorder Element -file WrapperData/EleSec1SS.out -time -eleRange 1  $NumEles section 1  fiber -35 -5 stressStrainTangent;
recorder Element -file WrapperData/EleSec2SS.out -time -eleRange 1  $NumEles section 1  fiber -25 -5 stressStrainTangent;
recorder Element -file WrapperData/EleSec3SS.out -time -eleRange 1  $NumEles section 1  fiber -15 -5 stressStrainTangent;
recorder Element -file WrapperData/EleSec4SS.out -time -eleRange 1  $NumEles section 1  fiber 5 -5 stressStrainTangent;
recorder Element -file WrapperData/EleSec5SS.out -time -eleRange 1  $NumEles section 1  fiber 35 -5 stressStrainTangent;

recorder Element -file WrapperData/EleForceSec1.out -time -eleRange 1  $NumEles section 1 forces;
recorder Element -file WrapperData/EleForceSec2.out -time -eleRange 1  $NumEles section 2  forces;
recorder Element -file WrapperData/EleForceSec3.out -time -eleRange 1  $NumEles section 3  forces;
recorder Element -file WrapperData/EleForceSec4.out -time -eleRange 1  $NumEles section 4  forces;
recorder Element -file WrapperData/EleForceSec5.out -time -eleRange 1  $NumEles section 5  forces;	


#display 2D deformation shape
set ViewScale 0.000001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D DeformedShape $ViewScale;

if {$ANALYSIS == "HasLoad"} {
puts "Point load";

pattern Plain 1 Linear {
load $EndNodeTag -300000 0 0 0 0 0;
#eleLoad -ele 1 2 -type -beamUniform -0 0 0;
}
constraints Plain;
numberer Plain;
system BandGeneral;
test NormUnbalance 1.0e-4 10 1;
algorithm Newton;
integrator LoadControl 0.1;
analysis Static;
analyze 10;
loadConst -time 0.0
}


if {$TANALYSIS == "HasThermo"} {


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
#load 1 -nodalThermal 800 $minusHalfD 400 $HalfD;

#load $MidNodeTag -nodalThermal 200 $minusHalfD 100 $HalfD

#load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD

#load $EndNodeTag -nodalThermal 800 $minusHalfD 400 $HalfD;
#eleLoad -range 1 $NumEles -type -beamThermal -source -node
#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $MidNodeTag 0.5 $EndNodeTag 1; 

eleLoad -range 1 1 -type -beamThermal 1000 -$HalfD 1000 $HalfD
#eleLoad -ele 1 2 -type -beamUniform -10 0 0;
}


constraints Plain;
numberer Plain;
system BandGeneral;
test NormUnbalance 1.0e-4 10 1;
algorithm Newton;
integrator LoadControl 0.01;
analysis Static;
analyze 100;
}; #end of thermo-mechanical analysis



