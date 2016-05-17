source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl
wipe;
file mkdir HTData;

SIFBuilder;

#AssignSection
SIFXBay 6 9 ;
SIFZBay 6 9;
SIFStorey 5 4  ;

AddMaterial steel 1 -type EC3 2.75e8  2e11;
AddMaterial concrete 2 -type EC2 0 30;

AddSection ISection 1 1 0.355 0.1715 0.0074 0.0115;  # $d $bf $tw $tf  UB356x171x51
AddSection ISection 2 1 0.6026 0.2276 0.0105 0.0148;  # $d $bf $tw $tf   UB610*229*101
AddSection ISection 3 1 0.3034 0.165 0.006 0.0102;  # $d $bf $tw $tf   UB305*165*40

AddSection ISection 4 1 0.254 0.254 0.0086 0.0142 -protected 0.02;;  # $d $bf $tw $tf   UC254*254*73

#AddSection Rect 2 1 0.1 0.1;  # $d $bf $tw $tf   UC203*203*46
AddSection  SlabSection 5 2 0.1;
AssignSection Xbeams 3 ;
#AssignSection Zbeams 1 ;
AssignSection Zbeams 1 -ZBay 1;
AssignSection Zbeams 2 -ZBay 2;
AssignSection columns 4;
AssignSection slabs 5;

AddSecBeam  XBeam 3 Compartment 111 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 211 -numBeams 1;
AddSecBeam  XBeam 3 Compartment 121 -numBeams 2;
AddSecBeam  XBeam 3 Compartment 221 -numBeams 2;
SetBC fixedJoint -locy 0;
 #set boundary condition
#AddLoad -joint 111 121 -load 5000 0 0;
AddLoad -member allslabs -load 0 -4950 0;
#AddLoad -liveLoad 2.0;
AddFire -compartment 111 -type standard;
#AddFire -compartment 111 -type EC1Local -origin 6 0 6 -HRR 3e6 -Dia 0.5;
#BuildModel elastic;
BuildModel -MeshCtrl 6 6 6;

# Define DISPLAY -------------------------------------------------------------
	set  xPixels 800;	# height of graphical window in pixels
	set  yPixels 800;	# height of graphical window in pixels
	set  xLoc1 100;	    # horizontal location of graphical window (0=upper left-most corner)
	set  yLoc1 60;	    # vertical location of graphical window (0=upper left-most corner)
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
SIFRecorder Joint -file Joint111.out -joint 111 disp;
SIFRecorder Member -file XBeam111.out -xBeam 111 Mideflect;
SIFRecorder Member -file XBeam121.out -xBeam 121 Mideflect;
SIFRecorder Member -file Slab111.out -slab 111 Mideflect;
time {
#SIFAnalyze SelfWeight -dt 0.5 Load -dt 0.5 Fire -dt 30 -duration 360 -output HTData;
SIFAnalyze  Fire -dt 15 -duration 3600 -output HTData;
#SIFAnalyze Load -dt 0.02;
}
#SIFAnalyze SelfWeight -dt 0.2 
print domain.out
