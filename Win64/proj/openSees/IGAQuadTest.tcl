#########################################################
#					                #
# Coarse-mesh cantilever beam analysis.  The beam is    #
# modeled with only 4 elements and uses anti-symmetry.  #
#							#
# ---> Basic units used are kN and meters		#
#							#
#########################################################

wipe

model BasicBuilder -ndm 2 -ndf 2;       #for IGA quad element:each node has 3 DoF, to include rotation.
source DisplayPlane.tcl;		# procedure for displaying a plane in model
source DisplayModel2D.tcl;		# procedure for displaying 2D perspective of model
# beam dimensions
set L1 1.0
set L2 2.0
set D1 0.25
set D2 0.5

# define number and size of elements 
set nElemX 1
set nElemY 1
set nElemT [expr $nElemX*$nElemY]
set sElemX [expr $L2/$nElemX]
set sElemY [expr $D2/$nElemY]

set nNodeX [expr $nElemX + 1]
set nNodeY [expr $nElemY + 1]
set nNodeT [expr $nNodeX*$nNodeY]

# create the nodes
set nid   1
set count 0.0

#Defining nodes
node 1 0 0;
node 2 $L1 0;
node 3 $L2 0;
node 4 0 $D1;
node 5 $L1 $D1;
node 6 $L2 $D1;
node 7 0 $D2;
node 8 $L1 $D2;
node 9 $L2 $D2;
# "Node Tag: $nid crdx: $crdx crdy:  $crdy    "


# boundary conditions
# fix 1   1 1 
# fix 3   1 1 
fixX 0.0   1 1


# define material
set matID 1
set E     20000
set nu    0.25
nDMaterial ElasticIsotropic $matID $E $nu

# create elements
set thick 1.0

element IGAquad 1 1 1 1 2 2 1 2 3 4 5 6 7 8 9 6 0. 0. 0. 1. 1. 1. 6 0. 0. 0. 1. 1. 1. 1 1 $thick "PlaneStress" $matID;
# element IGAquad 1 1 1 1 1 1 1 3 7 9 4 0. 0. 1. 1. 4 0. 0. 1. 1. 1 1 $thick "PlaneStress" $matID;
#element IGANURBS $ex $nex $ey $ney $order_x/y $controlPoints &patchKnotVect_x/y $nMult_x/y $thick "PlaneStress" $matID
#currently the overlapped patch cp is ignored

puts "element okay"
# create recorders
set step 0.1

recorder Node -time -file IGAd1p1m1.out -dT $step -nodeRange 1 $nNodeT -dof 1 2 disp
recorder Element -eleRange 1 $nElemT -time -file IGA1p1m1.out  -dT $step  stress
recorder Element -eleRange 1 $nElemT -time -file IGAe1p1m1.out  -dT $step  strain

#display 2D deformation shape
set ViewScale 0.2;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel2D DeformedShape $ViewScale;
# create loading
set P -1000.0;

pattern Plain 3 {Series -time {0 10 15} -values {0 1 1} -factor 1} { 
	load 3   0.0 [expr 0.25*$P]
	load 6   0.0 [expr 0.5*$P] 
	load 9   0.0 [expr 0.25*$P] 

}

# create analysis
puts "load okay"

integrator LoadControl 0.02
numberer RCM
system SparseGeneral
constraints Transformation
test NormDispIncr 1e-5 40 1
algorithm Newton
analysis Static

analyze 50

wipe