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
set NumEles 12;
set BeamLen 1140.0;


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
fix 1 1 1 0;
fix $EndNodeTag 0 1 0;


#define an elastic material with Tag=1 and E=2e11.
#uniaxialMaterial ElasticThermal 2 19200 1.4e-5;
#uniaxialMaterial ElasticThermal 3 206000 0.5e-5;
uniaxialMaterial SteelECThermal 1 EC3 235 206000;
#uniaxialMaterial Steel02Thermal 3 280 200000 0.00001;
#define fibred section; Two fibres: fiber $yLoc $zLoc $A $matTag 


###########################################################################
proc Wsection { secID matID d bf tf tw nfdw nftw nfbf nftf} {
	# ###################################################################
	# Wsection  $secID $matID $d $bf $tf $tw $nfdw $nftw $nfbf $nftf
	# ###################################################################
	# create a standard W section given the nominal section properties
	# written: Remo M. de Souza
	# date: 06/99
	# modified: 08/99  (according to the new general modelbuilder)
	# input parameters
	# secID - section ID number
	# matID - material ID number 
	# d  = nominal depth
	# tw = web thickness
	# bf = flange width
	# tf = flange thickness
	# nfdw = number of fibers along web depth 
	# nftw = number of fibers along web thickness
	# nfbf = number of fibers along flange width
	# nftf = number of fibers along flange thickness
  	
	set dw [expr $d - 2 * $tf]
	set y1 [expr -$d/2]
	set y2 [expr -$dw/2]
	set y3 [expr  $dw/2]
	set y4 [expr  $d/2]
  
	set z1 [expr -$bf/2]
	set z2 [expr -$tw/2]
	set z3 [expr  $tw/2]
	set z4 [expr  $bf/2]
  
	section fiberSecThermal  $secID  {
   		#                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
   		patch quadr  $matID  $nfbf $nftf   $y1 $z4   $y1 $z1   $y2 $z1   $y2 $z4
   		patch quadr  $matID  $nftw $nfdw   $y2 $z3   $y2 $z2   $y3 $z2   $y3 $z3
   		patch quadr  $matID  $nfbf $nftf   $y3 $z4   $y3 $z1   $y4 $z1   $y4 $z4
	}
}

#Wsection { 1 1 80 46 5.2 3.8 8 2 8 2}
Wsection { 1 1 320 160 10 10 8  2 8 2}
#define coordinate transforamtion: geomTransf $type $TransfTag; 
#three transformation types can be chosen: Linear, PDelta, Corotational)
geomTransf Corot 1 ; 

#define beam element: dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $TransfTag;
#"numIntgrPts" is the number of integration points along the element;
#"TransfTag" is pre-defined coordinate-transformation;
for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 5 1 1;
}

#define output

#define output
set MidSpanNode [expr 1+$NumEles/2];
set MidSpanEle [expr $NumEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
puts "$FiberLocBot , $FiberLocTop";
recorder Node -file WrapperData/N_MidSpan_DOF.out -time -node $MidSpanNode -dof 1 2 3 disp;


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


#display 2D deformation shape
set ViewScale 0.0001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
#DisplayModel2D DeformedShape $ViewScale;

if {$MechANALYSIS == "HasPoint"} {
set Load -100;
#define Uniform load
puts "Now applying uniform load"

pattern Plain 1 Linear {
    set UDL -1;
	#for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
		#eleLoad -ele $eleID -type -beamUniform $UDL 0 0;
	#}

   load	$MidSpanNode $Load 0 0;

};



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-3 300 1;
algorithm Newton;
#integrator LoadControl 0.5;
integrator DisplacementControl $EndNodeTag  1  1;	
analysis Static;
analyze 1000;
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