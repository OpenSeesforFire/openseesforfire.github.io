#Korean Fire Tests on CHS
#Experimental Study on Limiting Temperatures of Circular Hollow Sections

# Units mm, N/mm�
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
set TANALYSIS "Has0Thermo"; 

wipe;
file mkdir WrapperData;
model BasicBuilder -ndm 2 -ndf 3;
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model

#define node
set NumEles 24;
set BeamLen 3.0;
set EleLen [expr $BeamLen/$NumEles]
for {set NodeID 0} {$NodeID <= $NumEles} {incr NodeID} {
set locX [expr $NodeID*$EleLen];
set NodeTag [expr $NodeID+1];
node $NodeTag $locX 0;
}

set EndNodeTag [expr $NumEles+1];
set MidNodeTag [expr $NumEles/2+1];
#define boundary condition;
puts "Now uniform load";
fix 1 1 1 1;
#fix $EndNodeTag 0 1 0;


#define an elastic material with Tag=1 and E=2e11.
#uniaxialMaterial ElasticThermal 1 2e11 1.2e-5 -sSoft;
#uniaxialMaterial SteelECThermal 1 EC3 3e8 2e11;
uniaxialMaterial Steel01Thermal 1 3e8 2e11 0;
#define fibred section; Two fibres: fiber $yLoc $zLoc $A $matTag 
set NumFibers 16;
set b 0.2;
set d 0.05;
set secTag 1;
set NumInts 5;
section fiberSecThermal $secTag {
	for {set fiberID 1} {$fiberID<=$NumFibers} {incr fiberID} {
		set fiberLocy [expr ($d/$NumFibers)*($fiberID-0.5)-$d/2];
		set fiberArea [expr ($d/$NumFibers)*$b*0.5];
		fiber $fiberLocy -0.05 $fiberArea 1;
		fiber $fiberLocy 0.05 $fiberArea 1;
		#puts "$fiberLocy, $fiberArea";
		}
	}


#define coordinate transforamtion: geomTransf $type $TransfTag;
#three transformation types can be chosen: Linear, PDelta, Corotational)
#geomTransf Linear 1;
geomTransf Corotational 1 ; 
#geomTransf PDelta 1 ;

#define beam element: dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $TransfTag;
#"numIntgrPts" is the number of integration points along the element;
#"TransfTag" is pre-defined coordinate-transformation;
for {set eleID 1} {$eleID<= $NumEles} {incr eleID} {
	set NodeTag0 $eleID;
	set NodeTag1 [expr $eleID+1];
	#element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
	element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1 -mass 1000;
}

#define output

set MidSpanEle [expr $NumEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
puts "$FiberLocBot , $FiberLocTop";


#display 2D deformation shape
set ViewScale 0.0001;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel2D DeformedShape $ViewScale;

recorder Node -file WrapperData/EndNode_DOF.out -time -node $EndNodeTag -dof 1 2 disp;
recorder Element -file WrapperData/Element1Secm50T.out -time -ele 5 section 1  fiber -0.005 0.005 TempElong;

if {$ANALYSIS == "HasLoad"} {
#define Uniform load
puts "Now applying uniform load";

pattern Plain 1 Linear {

   load	$EndNodeTag  0  3000 0;

};


# constraints Plain;
# numberer Plain;
# system BandGeneral;
# test NormDispIncr 1e-5 100 2 ;
# algorithm Newton;
# integrator LoadControl 0.1;
# #integrator Explicitdifference��
# analysis Static;
# analyze 10;
# loadConst -time 0.0

constraints Transformation;
numberer RCM;
system Diagonal
algorithm Linear
integrator ExplicitDifference
analysis Transient
analyze 20000 0.00005



}

if {$TANALYSIS == "HasThermo"} {

puts "Thermal action"
# Define Thermal Load
set minusHalfD [expr -$d/2]
set HalfD [expr $d/2]

pattern Plain 2 Linear {
#Here we define nodal thermal action

#Here we offer different approaches to appliy the thermal action to element, 1: elemental thermal action, 2:using end nodal thermal action(each node), 3: using thermal action wrapper

#load $EndNodeTag -nodalThermal 800 $minusHalfD 400 $HalfD

#eleLoad -range 1 $NumEles -type -beamThermal -source -node

#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $EndNodeTag 1; 
eleLoad -range 1 $NumEles -type -beamThermal 1000 -$HalfD 1000 $HalfD

}


# constraints Plain;
# numberer Plain;
# system BandGeneral;
# test NormDispIncr 1.0e-3 100 1;
# algorithm Newton;
# integrator LoadControl 0.01;
# analysis Static;
# analyze 100;

# constraints Transformation;
# numberer RCM;
# system SparseSYM
# test NormDispIncr 1e-3 200 2
# algorithm NewtonLineSearch 0.75
# integrator Newmark 0.5 0.25

# analysis Transient
# puts "ok"
# analyze 2000 0.01

constraints Transformation;
numberer RCM;
system Diagonal
algorithm Linear
integrator ExplicitDifference
analysis Transient
analyze 20000 0.00005
wipe;
}
