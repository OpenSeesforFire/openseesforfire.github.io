source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl
wipe;
file mkdir HTData;

SIFBuilder display 2;

#AssignSection
SIFXBay 9 9 9;
SIFZBay 6 9 6;
SIFStorey 3 3 3 3;

AddMaterial steel 1 -type EC3 2.40e8  2e11;
AddMaterial stainlessSteel 3 -type Grade14571 2.40e8 5.20e8 2e11;
AddMaterial concrete 2 -type EC2 0 30e6;


AddSection ISection 1 1 0.355 0.1715 0.0074 0.0115;  # $d $bf $tw $tf  UB356x171x51
AddSection ISection 11 1 0.355 0.1715 0.0074 0.0115 -protected 0.02;  # $d $bf $tw $tf  UB356x171x51
AddSection ISection 2 1 0.6026 0.2276 0.0105 0.0148;  # $d $bf $tw $tf   UB610*229*101
AddSection ISection 3 1 0.3034 0.165 0.006 0.0102;  # $d $bf $tw $tf   UC305*165*40

AddSection ISection 13 1 0.254 0.254 0.0010 0.0173;  # $d $bf $tw $tf   UC254*254*89
AddSection ISection 14 1 0.254 0.254 0.0010 0.0173 -protected 0.02;;  # $d $bf $tw $tf   UC254*254*89
AddSection ISection 15 3 0.355 0.1715 0.0074 0.0115;  # $d $bf $tw $tf  UB356x171x51

AddSection  SlabSection 5 2 0.15;


#AddSection Rect 2 1 0.1 0.1;  # $d $bf $tw $tf   UC203*203*46

AssignSection 1 Xbeam -all;
AssignSection 11 Xbeam -ZBay 4 ;
AssignSection 1 Zbeam -all;
AssignSection 2 Zbeam -ZBay 2;
AssignSection 11 Zbeam -XBay 4;
AssignSection 14 column -all;
AssignSection 5 slab -all;
AssignSection 15 Zbeam 3302;
#AssignSection 6 column 2301;
AddSecBeam  XBeam 3 Compartment 1102 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 1202 -numBeams 2;
AddSecBeam  XBeam 3 Compartment 1302 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 2102 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 2202 -numBeams 2;
AddSecBeam  XBeam 3 Compartment 2302 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 3102 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 3202 -numBeams 2;
AddSecBeam  XBeam 3 Compartment 3302 -numBeams 1;
SetBC fixedJoint -Locy 0;
SetBC pinnedJoint -Locx 0;
 #set boundary condition
#AddLoad -joint 222 -load 0 600000 0;
AddLoad -member allslabs -load 0 -5000 0;
AddLoad -member allbeams -load 0 -4000 0;
#AddLoad -liveLoad 2.0;
AddFire -compartment 2302 -type standard;
AddFire -compartment 3302 -type standard;
#AddFire -compartment 1201 1301 2201 2301 -type EC1Local -origin 10 0 12 -HRR 3e6 -Dia 1.0;
#BuildModel elastic;
BuildModel -MeshCtrl 10 6 12;

# Define DISPLAY -------------------------------------------------------------
	set  xPixels 1600;	# height of graphical window in pixels
	set  yPixels 1600;	# height of graphical window in pixels
	set  xLoc1 -800;	    # horizontal location of graphical window (0=upper left-most corner)
	set  yLoc1 -800;	    # vertical location of graphical window (0=upper left-most corner)
	set ViewScale 10;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model

	######################################################################################
	##DisplayModel2D $ShapeType $dAmp $xLoc $yLoc $xPixels $yPixels $nEigen

	## display Node Numbers, Deformed or Mode Shape in 2D problem
	##			Silvia Mazzoni & Frank McKenna, 2006
	##
	## 	ShapeType : 	type of shape to display. # options: ModeShape , NodeNumbers , DeformedShape 
	## 	dAmp : 		relative amplification factor for deformations
	## 	xLoc,yLoc  : 	horizontal & vertical location in pixels of graphical window (0,0=upper left-most corner)
	## 	xPixels,yPixels :	width & height of graphical window in pixels
	## 	nEigen : 		if nEigen not=0, show mode shape for nEigen eigenvalue
	##
	#######################################################################################
	
DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1 $xPixels $yPixels 0

#SIFApply gravity
#SIFApply 
#SIFApply Fire
#SIFApply load 
SIFRecorder Joint -file Joint3302.out -joint 3302 disp;
SIFRecorder Joint -file Joint2302.out -joint 2302 disp;
SIFRecorder Joint -file Joint3402.out -joint 3402 disp;
SIFRecorder Joint -file Joint4302.out -joint 4302 disp;

SIFRecorder Member -file ZBeam3302.out -zBeam 3302 Mideflect;
#SIFRecorder Member -file Column2301.out -column 2301 Mideflect;
#SIFRecorder Member -file SecXBeam111.out -SecXBeam 111 Mideflect;
SIFRecorder Member -file Slab232.out -slab 2302 Mideflect;
#SIFRecorder Member -file Slab111all.out -slab 111 deflect;


recorder Element -file ElementColSecY.out -time -eleRange 877 888 section 3 fiber -0.127 0 stressStrainTangent;
#recorder Element -file ElementColSecYu.out -time -eleRange 1057 1059 section 3  fiber 0.127 0 stressStrainTangent;
recorder Element -file ElementColSecTEY.out -time -eleRange 877 888 section 3  fiber -0.127 0 TempElong;
#recorder Element -file ElementColSecTEYu.out -time -eleRange 321 330 section 3  fiber 0.127 0 TempElong;
recorder Element -file ElementColForces.out -time -eleRange 877 888 section 3  forces;
#recorder Element -file ElementColTopForces.out -time -ele 330 section 3  forces;
#recorder Element -file ElementColTopGForces.out -time -ele 330 globalForce;
#SIFAnalyze SelfWeight -dt 0.5 Load -dt 0.5 Fire -dt 30 -duration 360 -output HTData;
SIFAnalyze Load -dt 0.1 Fire -dt 15 -duration 3600 -output HTData;
#SIFAnalyze Load -dt 0.05;
#SIFAnalyze SelfWeight -dt 0.2 
#print domain.out
