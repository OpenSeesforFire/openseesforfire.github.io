source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl
wipe;
file mkdir HTData;

SIFBuilder display 2;

#AssignSection
SIFXBay 10 ;
SIFStorey 4 4 4 4 4 4 4 4 4 4 ;

AddMaterial steel 1 -type EC3 300e6  2.1e11;
AddMaterial concrete 2 -type EC2 0 35e6;

AddSection ISection 2 1 0.8349 0.2917 0.014 0.0188;  # $d $bf $tw $tf   UB838*292*176
AddSection ISection 3 1 0.5331 0.2093 0.0101 0.0156;  # $d $bf $tw $tf   UB533*210*92
AddSection ISection 4 1 0.4366 0.4122 0.0358 0.058 -protected 0.02;  # $d $bf $tw $tf   UC356*406*467
AddSection ISection 5 1 0.381 0.3948 0.0184 0.0302 -protected 0.02; # $d $bf $tw $tf   UC356*406*235
#AddSection Rect 2 1 0.1 0.1;  # $d $bf $tw $tf   UC203*203*46
AddSection  SlabSection 6 2 0.1 6.0;

AssignSection 3 beams ;
AssignSection 5 column -all;
AssignSection 6 slab -all;
SetBC fixedJoint -Locy 0;
SetBC pinnedJoint -Locx 10;
AddLoad -member allbeams -load 0 -30000;
#AddFire -compartment 1101 -type Standard -start 0;
AddFire -compartment 1104 -type Parametric -start 0;
AddFire -compartment 1105 -type Parametric -start 600;
AddFire -compartment 1106 -type Parametric -start 1200;
BuildModel -MeshCtrl 10 4 4 -pinned;

# Define DISPLAY -------------------------------------------------------------
	set  xPixels 1600;	# height of graphical window in pixels
	set  yPixels 2000;	# height of graphical window in pixels
	set  xLoc1 -800;	    # horizontal location of graphical window (0=upper left-most corner)
	set  yLoc1 0;	    # vertical location of graphical window (0=upper left-most corner)
	set ViewScale 1;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model

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
SIFRecorder Joint -file Joint1104.out -joint 1104 disp;
SIFRecorder Joint -file Joint1105.out -joint 1105 disp;
SIFRecorder Joint -file Joint1106.out -joint 1106 disp;
SIFRecorder Member -file XBeam1104.out -xBeam 1104 Mideflect;
SIFRecorder Member -file XBeam1105.out -xBeam 1105 Mideflect;
SIFRecorder Member -file XBeam1106.out -xBeam 1106 Mideflect;
#SIFRecorder Joint -file Joint221.out -joint 111 disp;
#SIFRecorder Member -file XBeam113.out -xBeam 1103 Mideflect;
#SIFRecorder Member -file XBeam114.out -xBeam 1104 Mideflect;
#SIFRecorder Member -file XBeam115.out -xBeam 1105 Mideflect;
#SIFRecorder Member -file col113.out -Column 1103 Mideflect;
#SIFRecorder Member -file col114.out -Column 1104 Mideflect;
#SIFRecorder Member -file col115.out -Column 1105 Mideflect;
#SIFRecorder Member -file XBeam121.out -xBeam 121 Mideflect;
#SIFRecorder Member -file SecXBeam111.out -SecXBeam 111 Mideflect;
#SIFRecorder Member -file Slab111.out -slab 111 Mideflect;
#SIFRecorder Member -file Slab111all.out -slab 111 deflect;
#SIFAnalyze Fire -dt 10 -duration 7200 -output HTData;
SIFAnalyze Load -dt 0.05 Fire -dt 10 -duration 7200 -output HTData;

print domain.output