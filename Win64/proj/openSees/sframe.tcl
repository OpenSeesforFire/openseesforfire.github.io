#Korean Fire Tests on CHS
#Experimental Study on Limiting Temperatures of Circular Hollow Sections

# Units mm, N/mm²
#		-------
#       |      |
#       |      |
#       |o-----|
#       |      |
#       |      |
#       |      |
#
#					y
#					|
#					|______x
#				   /
#				  /z
set ANALYSIS "HasLoad";
set TANALYSIS "HasThermo"; 

wipe;
file mkdir JointData;
model BasicBuilder -ndm 2 -ndf 3;
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model

set storeyH 3.0;
set bayL 6.0;
set numStorey 2;
set numBay 2;
set nx 6;
set ny 4;

#------------------Node Definition------------------------------
#For beam nodes
#For column nodes, first column NodeTag: 1000-1008, second column: 1100-1108
set NumColEles [expr $numStorey*$ny];  
set EleColLen [expr $storeyH/$ny];

for {set BayID 0} {$BayID <= $numBay} {incr BayID} {
for {set NodeID 0} {$NodeID <= $NumColEles} {incr NodeID} {
set locX [expr $BayID*$bayL];
set locY [expr $NodeID*$EleColLen];
set NodeTag [expr $NodeID+1000+100*$BayID];
node $NodeTag $locX $locY;
}
}

#For beam nodes, first floor beam NodeTag: 100-106, 110-116, second floor: 200-206, 210-216
set NumBeamEles $nx;  set BeamLen $bayL;
set EleBeamLen [expr $BeamLen/$NumBeamEles]

for {set storeyID 1} {$storeyID<=$numStorey} {incr storeyID} {
#loop over each storey
for {set BayID 0} {$BayID < $numBay} {incr BayID} {
for {set NodeID 0} {$NodeID <= $NumBeamEles} {incr NodeID} {
set NodeTag [expr int($NodeID+100*$storeyID+10*$BayID)];
set locX [expr $NodeID*$EleBeamLen+$BayID*$bayL];
set locY [expr $storeyH*$storeyID]
node $NodeTag $locX $locY;
}
}
}




#------------------Section Definition------------------------------
#define an elastic material with Tag=1 and E=2e11.
#uniaxialMaterial SteelECThermal 1 EC3 3e8 2e11;
uniaxialMaterial ElasticThermal 1 2e11 1.2e-5 -sSoft ;
uniaxialMaterial ElasticThermal 2 2e16 1.2e-5 -sSoft;
uniaxialMaterial ElasticThermal 3 2e3 1.2e-5 -sSoft;
uniaxialMaterial JointEPThermal 4 2e3 0.001 -0.002 -sSoft;  # materialTag, Elastic Modulus, yielding strain for positive deform, yielding strain for negative deform, softening as EC steel or concrete


#uniaxialMaterial SteelECThermal 1 EC3 3e8 2e11;
#uniaxialMaterial Steel01Thermal 1 3e8 2e11 0.00001;
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

#-----------------------Element Definition-----------------------------
#define beam element: dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts $secTag $TransfTag;
#"numIntgrPts" is the number of integration points along the element;
#"TransfTag" is pre-defined coordinate-transformation;

#Column ele: 1001-1008, 2001-2008; 
for {set BayID 0} {$BayID <= $numBay} {incr BayID} {
	for {set elei 1} {$elei<= $NumColEles} {incr elei} {
		#element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
		set NodeTag0 [expr $elei+1000+$BayID*100-1]; set NodeTag1 [expr $elei+1000+$BayID*100];
		set eleID [expr $elei+1000+$BayID*100];
		element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
	}
}
#Beam ele: 101-106; 201-206
for {set storeyID 1} {$storeyID<=$numStorey} {incr storeyID} {
	for {set BayID 0} {$BayID < $numBay} {incr BayID} {
		for {set elei 1} {$elei<= $NumBeamEles} {incr elei} {
		#element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
			set NodeTag0 [expr $elei+100*$storeyID+10*$BayID-1]; 
			set NodeTag1 [expr $elei+100*$storeyID+10*$BayID];
			set eleID [expr int($elei+100*$storeyID+10*$BayID)];
			element dispBeamColumnThermal $eleID $NodeTag0 $NodeTag1 $NumInts 1 1;
}	
}
}


#Generate boundary conditions---------------------------------------------//
#set EndNodeTag [expr 10];
#set MidNodeTag [expr 10+$NumBeamEles/2];
#define boundary condition;
for {set BayID 0} {$BayID <= $numBay} {incr BayID} {
set NodeTagB [expr 1000+100*$BayID]; 
fix $NodeTagB 1 1 1;
}
#fix $EndNodeTag 0 1 0;
 
#equalDOF 
for {set storeyID 1} {$storeyID <= $numStorey} {incr storeyID} {
for {set BayID 0} {$BayID < $numBay} {incr BayID} {
if {$BayID==0} {
	set leftBeamNodeID [expr $storeyID*100+$BayID*10];
	set leftColumnNodeID [expr 1000+$storeyID*$ny+$BayID*100];
	set JointEleId [expr 10000+ $storeyID*100+$BayID*10];
	element beamColumnJointThermal $JointEleId $leftColumnNodeID $leftBeamNodeID 2 2 2 -limit 1.0 1.0 200;
	
	set rightBeamNodeID [expr $storeyID*100+$NumBeamEles+$BayID*10];
	set rightColumnNodeID [expr 1100+$storeyID*$ny+$BayID*100];
    set JointEleId [expr 10000+ $storeyID*100+$BayID*10+$NumBeamEles];
	#equalDOF $rightBeamNodeID $rightColumnNodeID 1 2 3;
	element beamColumnJointThermal $JointEleId $rightBeamNodeID $rightColumnNodeID 2 2 3 -limit 1.0 1.0 0.01;   #limit  deform limits at each direction
} else {
	set leftBeamNodeID [expr $storeyID*100+$BayID*10];
	set leftColumnNodeID [expr 1000+$storeyID*$ny+$BayID*100];
	equalDOF $leftBeamNodeID $leftColumnNodeID 1 2 3;
	
	set rightBeamNodeID [expr $storeyID*100+$NumBeamEles+$BayID*10];
	set rightColumnNodeID [expr 1100+$storeyID*$ny+$BayID*100];
    equalDOF $rightBeamNodeID $rightColumnNodeID 1 2 3;
}
}
}

#element beamColumnJointThermal 1 1004 100 2 2 2;
#element beamColumnJointThermal 2 1104 106 2 2 2;
#equalDOF 1004 100 1 2 3;

#define output
puts "Now uniform load";
set MidSpanEle [expr 100+$NumBeamEles/2];
set FiberLocBot [expr ($d/$NumFibers)/2-($d/2)];
set FiberLocTop [expr ($d/$NumFibers)*($NumFibers-0.5)-$d/2];
puts "$FiberLocBot , $FiberLocTop";


#display 2D deformation shape
set ViewScale 5;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel2D DeformedShape $ViewScale;

recorder Node -file JointData/Node1004_DOF.out -time -node 1004 -dof 1 2 3 disp;
recorder Node -file JointData/Node100_DOF.out -time -node 100 -dof 1 2 3 disp;
recorder Node -file JointData/Node2004_DOF.out -time -node 1104 -dof 1 2 3 disp;
recorder Node -file JointData/Node106_DOF.out -time -node 106 -dof 1 2 3 disp;

recorder Element -file JointData/Element101GF.out -time -ele 101   globalForce;
recorder Element -file JointData/Element1004GF.out -time -ele 1004   globalForce;

recorder Element -file JointData/Element1D.out -time -ele 10100 deform;

if {$ANALYSIS == "HasLoad"} {
#define Uniform load
puts "Now applying uniform load";

pattern Plain 1 Linear {
#CREATE UNIFORM LOADS FOR BEAMS  ////$numStorey
set  UDL -10;
for {set level 1} {$level <= $NumBeamEles} {incr level 1} {
	for {set storeyID 1} {$storeyID<=1} {incr storeyID} {
	set eleID [expr $level+100*$storeyID];
	eleLoad -ele $eleID -type -beamUniform $UDL 0
	}
}
}


constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-4 100 2 ;
algorithm Newton;
integrator LoadControl 0.1;
analysis Static;
analyze 10;
loadConst -time 0.0
}

if {$TANALYSIS == "HasThermo"} {

puts "Thermal action"
# Define Thermal Load
set minusHalfD [expr -$d/2]
set HalfD [expr $d/2]

pattern Plain 2 Linear {
#Here we define nodal thermal action
#load 1 -nodalThermal 800 $minusHalfD 400 $HalfD;

#load $MidNodeTag -nodalThermal 200 $minusHalfD 100 $HalfD

#load $EndNodeTag -nodalThermal 0 $minusHalfD 0 $HalfD

#Here we offer different approaches to appliy the thermal action to element, 1: elemental thermal action, 2:using end nodal thermal action(each node), 3: using thermal action wrapper

#load $EndNodeTag -nodalThermal 800 $minusHalfD 400 $HalfD

#eleLoad -range 1 $NumEles -type -beamThermal -source -node

#eleLoad -range 1 $NumEles -type -ThermalWrapper -nodeLoc 1 0 $EndNodeTag 1; 
eleLoad -range 10100 10100 -type -beamThermal 1000 -$HalfD 1000 $HalfD
eleLoad -range 101 106 -type -beamThermal 1000 -$HalfD 1000 $HalfD

}


constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1.0e-5 100 1;
algorithm Newton;
integrator LoadControl 0.01;
analysis Static;
analyze 100;

}
#wipe;