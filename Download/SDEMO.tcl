source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl

wipe;
SIFBuilder; 
SIFXBay 6 9 ; 
SIFZBay 6 9;
SIFStorey 5 4;
AddMaterial steel 1 -type EC3 3e8 2e11;
AddMaterial concrete 2 -type EC2 0 30e6;

AddSection ISection 1 1 0.203 0.102 0.0054 0.009;
AddSection ISection 3 1 0.254 0.254 0.0086 0.0142;
AddSection  SlabSection 5 2 0.1;
AssignSection 1 beams;
AssignSection 3 column -all;
AssignSection 5 slab -all;
SetBC fixedJoint -Locy 0; 
AddLoad -member allslabs -load 0 -1000 0;
AddFire -compartment 1101 -type standard;
BuildModel -MeshCtrl 6 6 6;
SIFRecorder Joint -file Joint111.out -joint 1101 disp; 

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
SIFAnalyze Load -dt 0.2 Fire -dt 30 -duration 1800 -output HTData; 
print domain.out
