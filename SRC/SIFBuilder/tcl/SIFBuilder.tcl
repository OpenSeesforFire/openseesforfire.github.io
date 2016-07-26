source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
source DisplayModel3D.tcl

file mkdir HTData;     #define the directory for storing data
SIFBuilder;  #initialise SIFBuilder, (SIFBuilder frame) is accepted for defining frame only without slab

#[BUILDING INFO]
SIFXBay 6 9 ;  #XBAY SPAN|<----6m------->|<----------9m------->| along global x direction
SIFZBay 6 9;   #ZBAY SPAN|<----6m------->|<----------9m------->| along global Z direction
SIFStorey 5 4; #Storey Height|<----5m---->|<---4m--->| along global y direction

AddSecBeam xBeam -compartment 111 -num 2;

#[DEFINE MATERIAL AND SECTION]
AddMaterial steel 1 -type EC3 3e8  2e11;   #E0 : 3e8 , fy: 2e11, EN-1993-1-2 Steel Material
AddMaterial concrete 2 -type EC2 0 30;  #moisture ratio:0 , fc :30, EN-1992-1-2 Concrete material

AddSection ISection 1 1 0.203 0.102 0.0054 0.009;  # $d $bf $tw $tf  UB203x102x23
AddSection ISection 2 1 0.203 0.203 0.007 0.011;  # $d $bf $tw $tf   UC203*203*46
#AddSection Rect 2 1 0.1 0.1;  # $d $bf $tw $tf
AddSection  SlabSection 3 2 0.1 -sRatx 1 200 -sRaty 1 200;  #thickness Reinforcement ratio $steelPar $steelArea

#[ASSIGN SECTION]
AssignSection beams 1;
AssignSection columns 2;
AssignSection slabs 3;

#[DEFINE BC AND LOAD]
SetBC fixedJoint -locy 0;     #set boundary condition
#AddLoad -joint 111 121 -load 5000 0 0;
AddLoad -member allslabs -load 0 -1000 0;
AddFire -compartment 111 -type standard;

#[BUILD MODEL]
BuildModel -MeshCtrl 6 6 6;        #Number of Eles meshed for each member (along global X,Y,Z)

#[Define DISPLAY]
set  xPixels 800;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 100;	    # horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 60;	    # vertical location of graphical window (0=upper left-most corner)
set ViewScale 1;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D DeformedShape $ViewScale $xLoc1 $yLoc1 $xPixels $yPixels 0

#[Define SIFRECORDER]
SIFRecorder Joint -file Joint111.out -joint 111 disp;                  #Define recorders for obtaining displacement of the specific joint
SIFRecorder Member -file XBeam111.out -xBeam 111 Mideflect;			   #Recorder for mid-span deflection of XBeams
SIFRecorder Member -file Slab111.out -slab 111 Mideflect;			   #Recorder for centre deflection of slab

#[RUN ANALYSES]
#SIFAnalyze SelfWeight -dt 0.5 Load -dt 0.5 Fire -dt 30 -duration 360 -output HTData;
SIFAnalyze Load -dt 0.2 Fire -dt 30 -duration 1800 -output HTData;       #First apply mechanical load and then fire action
#SIFAnalyze Load -dt 0.02;
print domain.out