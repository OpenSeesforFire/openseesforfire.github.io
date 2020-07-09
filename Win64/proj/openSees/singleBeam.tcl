wipe
model BasicBuilder -ndm 3 -ndf 6;

set h_rib 70.0
set h_UB305 303.4
set h_slab 70.0

for {set i 1} {$i <= 49} {incr i 1} {
	set x [expr ($i-1)*-300]
	set z [expr -0.5*$h_UB305 - $h_rib -0.5*$h_slab]
	
	puts "node number is $i";
	puts "x = $x";
	puts "z = $z"; 
	node $i $x 3000 $z;
};

uniaxialMaterial Steel01Thermal 2 300.0 2.0e5 0.01; # S275 - range of 291 - 318 MPa

set MatTag 2
set SecTag 2
set UB305 $SecTag
section fiberSecThermal $SecTag -GJ 2e11 { 
#fiber       $yLoc   $zLoc          $A  			   $matTag
fiber		150.00000000	100.0	561		$MatTag
fiber		146.60000000	-100.0	561		$MatTag
fiber		143.20000000	0	561		$MatTag
fiber		113.20000000	0	339.6	$MatTag
fiber		 56.6			0	339.6	$MatTag
fiber		-56.60000000	0	339.6	$MatTag
fiber		-113.20000000	0	339.6	$MatTag
fiber		-143.20000000	0	561		$MatTag
fiber		-146.60000000	100.0	561		$MatTag
fiber		-150.00000000	-100.0	561		$MatTag
}

geomTransf Corotational 2 0 1 0;


for {set i 2} {$i <= 49} {incr i 1} {
	set node1 [expr int($i-1)]
	set node2 $i
	set eleID $node1
	puts "node 1 is number: $node1"
	puts "node 2 is number: $node2"
	element dispBeamColumnThermal $eleID $node1 $node2 5 $UB305 2;
};


fix 1 1 1 1 1 0 0;
fix 49 1 1 1 0 0 0;
fix 20 1 1 1 0 0 0;
fix 35 1 1 1 0 0 0;

pattern Plain 1 Linear {
	for {set i 2} {$i < 49} {incr i 1} {
	load $i 	0.0 	0.0 	-50	 0.0 	0.0 	0.0;
	}
}
print nodes.dat -node
print elements.dat -ele
recorder Node -file RibU3.out -time -nodeRange 1 49 -dof 3 disp

record

constraints Transformation;     		
numberer RCM;			
system BandGeneral;		
test NormDispIncr 0.1 100 1;	
algorithm Newton;					
integrator LoadControl 0.1;				
analysis Static;
analyze 10;					
loadConst

wipe