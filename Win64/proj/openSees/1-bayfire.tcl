wipe
puts "System"
model basic -ndm 3 -ndf 6
# Load display engine
 source DisplayPlane.tcl
 source DisplayModel2D.tcl
 source DisplayModel3D.tcl
 
file mkdir oneBayData;
puts "restraint"
node 1 3.999E+003 0.000E+000 0.000E+000
node 2 3.999E+003 0.000E+000 4.500E+003
node 3 3.999E+003 3.999E+003 0.000E+000
node 4 3.999E+003 3.999E+003 4.500E+003
node 5 0.000E+000 3.999E+003 0.000E+000
node 6 0.000E+000 3.999E+003 4.500E+003
node 7 0.000E+000 0.000E+000 0.000E+000
node 8 0.000E+000 0.000E+000 4.500E+003
node 9 1.333E+003 3.999E+003 4.500E+003
node 10 2.666E+003 3.999E+003 4.500E+003
node 11 1.333E+003 0.000E+000 4.500E+003
node 12 2.666E+003 0.000E+000 4.500E+003
node 13 3.999E+003 1.333E+003 4.500E+003
node 14 3.999E+003 2.666E+003 4.500E+003
node 15 0.000E+000 2.666E+003 4.500E+003
node 16 0.000E+000 1.333E+003 4.500E+003
puts "rigidDiaphragm"
puts "mass"
mass 1 4.861E-001 4.861E-001 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 2 4.218E+000 4.218E+000 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 3 4.861E-001 4.861E-001 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 4 4.218E+000 4.218E+000 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 5 4.861E-001 4.861E-001 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 6 4.218E+000 4.218E+000 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 7 4.861E-001 4.861E-001 0.000E+000 0.000E+000 0.000E+000 0.000E+000
mass 8 4.218E+000 4.218E+000 0.000E+000 0.000E+000 0.000E+000 0.000E+000
puts "node"
fix 1 1 1 1 1 1 1;
fix 3 1 1 1 1 1 1;
fix 5 1 1 1 1 1 1;
fix 7 1 1 1 1 1 1;
puts "material"
set gt [expr 0/34000*0*2]
set gc [expr 35/34000*35*6]
# uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
uniaxialMaterial Steel01Thermal 1 308 2.1e5 0.01; 
# uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
uniaxialMaterial ConcreteECThermal 2 -35 -0.003 -700 -0.02 0.1 0 0;
#nDMaterial CDPPlaneStressThermal $matTag $E0 $Poisson_ratio $ft $fc $gt $gc
#nDMaterial  CDPPlaneStressThermal 4 34000 0.2 0 35 $gt $gc
#nDMaterial PlateRebarThermal matTag uniaxialMatTag orientation(degrees)
#nDMaterial PlateRebarThermal 5 1 0
#nDMaterial PlateRebarThermal 6 1 90
#nDMaterial   PlateFromPlaneStressThermal    matTag   	nDMatTag    forStability
#nDMaterial PlateFromPlaneStressThermal 7 4 1e9
puts "section"
#set SectionType Elastic ;		# options: Elastic FiberSection
set SectionType FiberSection ;		# options: Elastic FiberSection
if {$SectionType == "Elastic"} {
#BEAM250X600 
#section Elastic $secTag $E $A $Iz $Iy $G $J <$alphaY $alphaZ>
section Elastic 1 34000 150000 781250000 4500000000 14166.67 2306747907.473717 
#Colmn 300*300
section Elastic 2 34000 90000 675000000 675000000 14166.67 1140750000 
} elseif {$SectionType == "FiberSection"} {
##BEAM250X600 
section FiberThermal 1 -GJ 3.27E+13 {
fiber	-1.25E+02	-3.00E+02	6.00E+03	2
fiber	-6.25E+01	-3.00E+02	6.00E+03	2
fiber	0.00E+00	-3.00E+02	6.00E+03	2
fiber	6.25E+01	-3.00E+02	6.00E+03	2
fiber	1.25E+02	-3.00E+02	6.00E+03	2
fiber	-1.25E+02	-1.50E+02	6.00E+03	2
fiber	-6.25E+01	-1.50E+02	6.00E+03	2
fiber	0.00E+00	-1.50E+02	6.00E+03	2
fiber	6.25E+01	-1.50E+02	6.00E+03	2
fiber	1.25E+02	-1.50E+02	6.00E+03	2
fiber	-1.25E+02	0.00E+00	6.00E+03	2
fiber	-6.25E+01	0.00E+00	6.00E+03	2
fiber	0.00E+00	0.00E+00	6.00E+03	2
fiber	6.25E+01	0.00E+00	6.00E+03	2
fiber	1.25E+02	0.00E+00	6.00E+03	2
fiber	-1.25E+02	1.50E+02	6.00E+03	2
fiber	-6.25E+01	1.50E+02	6.00E+03	2
fiber	0.00E+00	1.50E+02	6.00E+03	2
fiber	6.25E+01	1.50E+02	6.00E+03	2
fiber	1.25E+02	1.50E+02	6.00E+03	2
fiber	-1.25E+02	3.00E+02	6.00E+03	2
fiber	-6.25E+01	3.00E+02	6.00E+03	2
fiber	0.00E+00	3.00E+02	6.00E+03	2
fiber	6.25E+01	3.00E+02	6.00E+03	2
fiber	1.25E+02	3.00E+02	6.00E+03	2
fiber	-1.00E+02	2.70E+02	2.54E+02	1
fiber	1.00E+02	2.70E+02	2.54E+02	1
fiber	-1.00E+02	-2.70E+02	2.54E+02	1
fiber	-4.00E+01	-2.70E+02	2.54E+02	1
fiber	4.00E+01	-2.70E+02	2.54E+02	1
fiber	1.00E+02	-2.70E+02	2.54E+02	1


}
##COL300X300 
section FiberThermal 2 -GJ 1.62E+13 {
fiber	-1.50E+02	-1.50E+02	3.60E+03	2
fiber	-9.00E+01	-1.50E+02	3.60E+03	2
fiber	0.00E+00	-1.50E+02	3.60E+03	2
fiber	9.00E+01	-1.50E+02	3.60E+03	2
fiber	1.50E+02	-1.50E+02	3.60E+03	2
fiber	-1.50E+02	-9.00E+01	3.60E+03	2
fiber	-9.00E+01	-9.00E+01	3.60E+03	2
fiber	0.00E+00	-9.00E+01	3.60E+03	2
fiber	9.00E+01	-9.00E+01	3.60E+03	2
fiber	1.50E+02	-9.00E+01	3.60E+03	2
fiber	-1.50E+02	0.00E+00	3.60E+03	2
fiber	-9.00E+01	0.00E+00	3.60E+03	2
fiber	0.00E+00	0.00E+00	3.60E+03	2
fiber	9.00E+01	0.00E+00	3.60E+03	2
fiber	1.50E+02	0.00E+00	3.60E+03	2
fiber	-1.50E+02	9.00E+01	3.60E+03	2
fiber	-9.00E+01	9.00E+01	3.60E+03	2
fiber	0.00E+00	9.00E+01	3.60E+03	2
fiber	9.00E+01	9.00E+01	3.60E+03	2
fiber	1.50E+02	9.00E+01	3.60E+03	2
fiber	-1.50E+02	1.50E+02	3.60E+03	2
fiber	-9.00E+01	1.50E+02	3.60E+03	2
fiber	0.00E+00	1.50E+02	3.60E+03	2
fiber	9.00E+01	1.50E+02	3.60E+03	2
fiber	1.50E+02	1.50E+02	3.60E+03	2
fiber	-1.25E+02	-1.25E+02	3.14E+02	1
fiber	1.25E+02	-1.25E+02	3.14E+02	1
fiber	-1.25E+02	1.25E+02	3.14E+02	1
fiber	1.25E+02	1.25E+02	3.14E+02	1
}
};
#section LayeredShellThermal   secTag  numoflayers mat1Tag  mat1Thickness mat2Tag  mat2Thickness ... matnTag  matnThickness
#section  LayeredShellThermal   3 10 7 25 7 25 7 25 5 25 6 25 7 25 7 25 7 25 7 25 7 25 ;

puts "transformation"
geomTransf Linear 1 1.000 0.000 0.000 
geomTransf Linear 2 1.000 0.000 0.000 
geomTransf Linear 3 1.000 0.000 0.000 
geomTransf Linear 4 1.000 0.000 0.000 
geomTransf Linear 5 0.000 0.000 1.000 
geomTransf Linear 6 0.000 0.000 1.000 
geomTransf Linear 7 0.000 0.000 1.000 
geomTransf Linear 8 0.000 0.000 1.000 
geomTransf Linear 9 0.000 0.000 1.000 
geomTransf Linear 10 0.000 0.000 1.000 
geomTransf Linear 11 0.000 0.000 1.000 
geomTransf Linear 12 0.000 0.000 1.000 
geomTransf Linear 13 0.000 0.000 1.000 
geomTransf Linear 14 0.000 0.000 1.000 
geomTransf Linear 15 0.000 0.000 1.000 
geomTransf Linear 16 0.000 0.000 1.000 
puts "Weight"
set GammaConcrete 2.452e-5;   			#Reinforced-Concrete weight density (weight per volume) 
set HCol 300;
set BCol 300;
set HBeam 600;
set BBeam 250;
set DLfactor 1.35;				# scale dead load up a little
set QdlCol [expr $GammaConcrete*$HCol*$BCol*$DLfactor];	# self weight of Column, weight per length
set QBeam [expr $GammaConcrete*$HBeam*$BBeam*$DLfactor];	# self weight of Beam, weight per length
set Tslab 250;			# 6-inch slab
set Lslab [expr 4000/2]; 			    # slab extends a distance of $LGird/2 in/out of plane
set Qslab [expr $GammaConcrete*$Tslab*$Lslab*$DLfactor]; 
set QdlBeam [expr $Qslab + $QBeam]; 	# dead load distributed along beam (one-way slab)

puts "element"
element dispBeamColumnThermal 1 1 2 3 2 1;
element dispBeamColumnThermal 2 3 4 3 2 2;
element dispBeamColumnThermal 3 5 6 3 2 3;
element dispBeamColumnThermal 4 7 8 3 2 4;
element dispBeamColumnThermal 5 2 13 3 1 5;
element dispBeamColumnThermal 6 13 14 3 1 6;
element dispBeamColumnThermal 7 14 4 3 1 7;
element dispBeamColumnThermal 8 6 9 3 1 8;
element dispBeamColumnThermal 9 9 10 3 1 9;
element dispBeamColumnThermal 10 10 4 3 1 10;
element dispBeamColumnThermal 11 6 15 3 1 11;
element dispBeamColumnThermal 12 15 16 3 1 12;
element dispBeamColumnThermal 13 16 8 3 1 13;
element dispBeamColumnThermal 14 8 11 3 1 14;
element dispBeamColumnThermal 15 11 12 3 1 15;
element dispBeamColumnThermal 16 12 2 3 1 16;
#nDMaterial PlateFiberThermal 601 4	
#section PlateFiberThermal 3 601 0.1
#element ShellNLDKGQThermal 17 2 4 6 8 3;
puts "recorder"
recorder Node -file oneBayData/node1reaction.out -time -node 1 -dof 1 2 3 reaction;
#recorder Node -file oneBayData/node3reaction.out -time -node 3 -dof 3 reaction;
#recorder Node -file oneBayData/node5reaction.out -time -node 5 -dof 3 reaction;
#recorder Node -file oneBayData/node7reaction.out -time -node 7 -dof 3 reaction;
recorder Node -file oneBayData/node13disp.out -time -node 13 -dof 1 2 3 4 5 6 disp;
recorder Node -file oneBayData/node14disp.out -time -node 14 -dof 1 2 3 4 5 6 disp;
recorder Node -file oneBayData/node2TopColumndisp.out -time -node 2 -dof 1 2 3 4 5 6 disp;
#recorder Element -file ShellData/EleForceSec1Temp.out -time -ele $MidEle   material 1 fiber 1 TempAndElong;
#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 100;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels


puts "loading"
# Load Case = DEAD
pattern Plain 1 Linear {
eleLoad -range 5 16 -type -beamUniform 0 -47.25 0;
eleLoad -range 5 7 -type -beamUniform 0 -$QdlBeam 0
eleLoad -range 8 10 -type -beamUniform 0 -$QBeam 0
eleLoad -range 11 13 -type -beamUniform 0 -$QdlBeam 0
eleLoad -range 14 16 -type -beamUniform 0 -$QBeam 0
eleLoad -range 1 4 -type -beamUniform 0 0 -$QdlCol	
};
	constraints Plain;					# how it handles boundary conditions
	##numberer Plain;	
	numberer RCM;					# renumber dof's to minimize band-width (optimization)
	#system BandGeneral;					# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	system UmfPack
	test NormUnbalance 1.0e-5 20 1;		# tolerance, max iterations
	algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
	
	integrator LoadControl 0.1;				# determine the next time step for an analysis, # apply gravity in 10 steps
	analysis Static;					# define type of analysis static or transient
	analyze 10;
	
puts "Gravity Analysis Completed!!!!"
loadConst -time 0.0
puts "Thermal Load"

loadConst -time 0.0	;
set T 1000; 	#Max Temp - DEG CELCIUS 
set Y1 300;	    #top fiber of beam
set Y2 -300;	#Bottom fiber of beam


#pattern Plain $PatternTag Linear { eleLoad -ele $eleTag -type -beamThermal $MaxTemp $ExtremeFiberLoc1 $MaxTemp $ExtremeFiberLoc2 };
pattern Plain 2 Linear { 
#load 14 -nodalThermal $T $Y1 $T $Y2;
eleLoad -range 5 16 -type -beamThermal $T $Y2 $T $Y1;

};

puts "analysis"
	#constraints Transformation
	#constraints Penalty 1e10 1e10
	constraints Plain;					# how it handles boundary conditions
	numberer Plain;						# renumber dof's to minimize band-width (optimization)
	#system BandGeneral;					# how to store and solve the system of equations in the analysis (large model: try UmfPack)
	system UmfPack
	test NormUnbalance 1.0e-4 80 1;		# tolerance, max iterations
	#test NormDispIncr 1e-6 10
	#test EnergyIncr 1.0e-8 50
	algorithm Newton;					# use Newton's solution algorithm: updates tangent stiffness at every iteration
	integrator LoadControl 0.01;
    analysis Static;
    analyze 100; #100steps to apply thermal action
puts "Fire done"
