# units: m,N
wipe;	
				
				
set ANALYSIS "HasPoint";
set TANALYSIS "HasThermo"; 

####################################################################################################################

# choose analysis way 

#set analysis "M3";
set analysis "M4";

set material "elastic";
#set material "plastic";
####################################################################################################################

model BasicBuilder -ndm 3 -ndf 6;

file mkdir ShellData;

# source DisplayPlane.tcl
# source DisplayModel2D.tcl
# source DisplayModel3D.tcl

#8.72
#Concrete model
#these should both be even, number of elements per edge
set nx 20;
set ny 8
set slabT 0.1;
set slabF 0.05;
set slabB 0.4;
set slabL 3.0;
set elemx [expr $slabL/$nx]
set elemy [expr $slabB/$ny]
set offset [expr ($slabT*0.5-$slabF*0.5)]
set UDL 0.5E3;

#loaded nodes
set midNode [expr (($ny/2+1)*100+($nx/2+1))];

#################################################################################################################

# creat nodes
for {set i 0} {$i <= $nx} {incr i 1} {
   set x [expr $i*$elemx]
   for {set j 0} {$j <= $ny} {incr j 1} {
       set y [expr $j*$elemy]
       set nodeID [expr ($j+1)*100+$i+1]
	   node $nodeID $x $y 0
	}
}	
#################################################################################################################

# creat materials

# concrete

if {$material == "elastic"} {
    nDMaterial ElasticIsotropic3DThermal 2 2e10 0.2 0 1.4e-5;
    nDMaterial PlateFiberThermal 4 2;
} elseif {$material == "plastic"} {
    set gt [expr 3.0e6/1.79e10*3.0e6*2];
    set gc [expr 30e6/1.79e10*30e6*6];
    nDMaterial  CDPPlaneStressThermal 100  1.79e10 0.2  3.0e6  30e6 $gt $gc;
    nDMaterial   PlateFromPlaneStressThermal    4   100  10e9;
}

# Rebar

#Rebar mesh
uniaxialMaterial SteelECThermal 1 EC2NH 3.45e8 2e11;
#uniaxialMaterial ElasticThermal 1 2e11 1.2e-5;
nDMaterial PlateRebarThermal 3 1 0;
nDMaterial PlateRebarThermal 5 1 90;

#Rebar J2 steel layer
set Es 2e11; set v 0.3;
set K [expr $Es/2.0/(1-$v)];
set G [expr $Es/2.0/(1+$v)]; 
#nDMaterial  J2PlasticityThermal 23 22 $K $G 3.45e8 4.00e8 0.1 0;
#nDMaterial   PlateFiberThermal    24   23;
nDMaterial  J2PlaneStressThermal 23 21 2e11 0.3 3.45e8 4.00e8 0 0;
nDMaterial   PlateFromPlaneStressThermal    24   23   20e10;

#################################################################################################################

# Creat section

#section LayeredShellThermal  2  14  4  0.01  4 0.009607301 3 0.000392699 5 0.000392699 4 0.009607301 4  0.01  4  0.01 4  0.01  4  0.01 4 0.009607301 3 0.000392699 5 0.000392699 4 0.009607301 4  0.01 ;
section LayeredShellThermal  2  5  4  0.01  4 0.01 4 0.01 4  0.01  4  0.01;
section LayeredShellThermal  3  10  4  0.01  4 0.01 4 0.01 4  0.01  4  0.01 4  0.01  4  0.01  4 0.01 4 0.01  4 0.01;
section LayeredShellThermal  4 -offset $offset 10  4  0.01  4 0.01 4 0.01 4  0.01  4  0.01 4  0.01  4  0.01  4 0.01 4 0.01  4 0.01;
section CompositeShellThermal 10 2 0.5 4 0.5; 

#section LayeredShellThermal  3  4  14  0.02 14 0.02  14 0.02  24 0.02 24 0.02;

#################################################################################################################
puts "here0";

# Creat elements
if {$analysis == "M3"} {
    for {set i 1} {$i <= $nx} {incr i 1} {
        for {set j 1} {$j <= $ny} {incr j 1} {
		   if {$j == 3 || $j == 4 || $j == 5 || $j ==6} {
			      set sec 4
			} else {
                  set sec 2
            } 
		  set node1 [expr $j*100+$i+1]
		  set node2 [expr ($j+1)*100+$i+1]
		  set node3 [expr $node2-1]
		  set node4 [expr $node1-1]
		  set elementID $node4
		  element ShellNLDKGQThermal $elementID $node1 $node2 $node3 $node4 $sec
		}	
	}
}  elseif {$analysis == "M4"} {
    for {set i 1} {$i <= $nx} {incr i 1} {
        for {set j 1} {$j <= $ny} {incr j 1} {
		  set node1 [expr $j*100+$i+1]
		  set node2 [expr ($j+1)*100+$i+1]
		  set node3 [expr $node2-1]
		  set node4 [expr $node1-1]
		  set elementID $node4
		  element ShellNLDKGQThermal $elementID $node1 $node2 $node3 $node4 10
		}
	}
}  

#section LayeredShellThermal  2  10  4  0.005 4 0.005  4 0.005  4  0.005  4 0.005 4  0.005  4  0.005 4  0.005  4  0.005 4  0.005 ;
#block2D $nx $ny 1 1 ShellNLDKGQThermal 2  ShellMITC4Thermal ShellMITC4GNLThermal
#################################################################################################################

# Creat constraints

#fully simply supported
fixX 0  1 0 1 1 1 1 ;
#fixX 3  0 0 1 1 0 1 ;
#fixX 6.0   0 1 0 1 0 1 ;
#fixX 3.0   0 1 0 1 0 1 ;
#fixY 0  0 1 0 1 0 1 ;
#fixY 0.2   0 1 0 1 0 1 ;

#fix 6   1 1 1 1 1 1 ;
#fix 11   1 0 1 1 1 1 ;
#fixX 3.0 0 0 1 0 0 0 
#fix 1 1 1 1 0 0 0;
#fix [expr 1+($nx+1)*$ny/2]  0 1 0 0 0 0 ;
#fix [expr ($nx+1)] 0 0 1 0 0 0 ;
#fix [expr ($nx+1)*($ny+1)] 0 0 1 0 0 0;
#fixY 0. 0 1 0 0 0 0;
#fixY 1 0 1 0 0 0 0;

#################################################################################################################

set midSlabsb [expr $nx/2+($nx+1)*$ny/2];
#display 3D deformation shape
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.01;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
# DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

# set MidEle [expr 2+($nx)*$ny/2];
# set EndEle [expr 1+$nx];
set MidEle [expr ($ny/2)*100+($nx/2)];
set EndEle [expr $ny*100+$nx]
recorder Node -file ShellData/DFreeSlabxT.out -time -nodeRange [expr ($ny/2+1)*100+1] [expr ($ny/2+1)*100+($nx+1)] -dof 1  disp; #nodes in Y=0.4;
recorder Node -file ShellData/DFreeSlabyT.out -time -nodeRange [expr ($ny/2+1)*100+1] [expr ($ny/2+1)*100+($nx+1)] -dof 2 disp; #nodes in Y=0.4;
recorder Node -file ShellData/DFreeSlabzT.out -time -node [expr ($ny/2+1)*100+($nx+1)] -dof 3 disp; #nodes in x=3.0, y=0.4

recorder Node -file ShellData/DFreeSlabzTUP.out -time -nodeRange [expr 1*100+1] [expr 1*100+($nx+1)] -dof 3 disp; #nodes in Y=0
recorder Node -file ShellData/DFreeSlabzTDo.out -time -nodeRange [expr ($ny+1)*100+1] [expr ($ny+1)*100+($nx+1)] -dof 3 disp; #nodes in Y=0.8

recorder Element -file ShellData/EleForce.out -time -ele $MidEle  forces;	

recorder Element -file ShellData/EleForceSec1sigma.out -time -ele $MidEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/EleForceSec1Eps.out -time -ele $MidEle  material 1 fiber 1 strain;
recorder Element -file ShellData/EleForceSec1Temp.out -time -ele $MidEle   material 1 fiber 1 TempAndElong;

recorder Element -file ShellData/EleForceSec2sigma.out -time -ele $MidEle  material 1 fiber 2 stress;	
recorder Element -file ShellData/EleForceSec2Eps.out -time -ele $MidEle  material 1 fiber 2 strain;
recorder Element -file ShellData/EleForceSec2Temp.out -time -ele $MidEle   material 1 fiber 2 TempAndElong;

recorder Element -file ShellData/EleForceSec3sigma.out -time -ele $MidEle  material 1 fiber 3 stress;	
recorder Element -file ShellData/EleForceSec3Eps.out -time -ele $MidEle  material 1 fiber 3 strain;
recorder Element -file ShellData/EleForceSec3Temp.out -time -ele $MidEle   material 1 fiber 3 TempAndElong;

recorder Element -file ShellData/EleForceSec4sigma.out -time -ele $MidEle  material 1 fiber 4 stress;	
recorder Element -file ShellData/EleForceSec4Eps.out -time -ele $MidEle  material 1 fiber 4 strain;
recorder Element -file ShellData/EleForceSec4Temp.out -time -ele $MidEle   material 1 fiber 4 TempAndElong;

recorder Element -file ShellData/EleForceSec5sigma.out -time -ele $MidEle  material 1 fiber 5 stress;	
recorder Element -file ShellData/EleForceSec5Eps.out -time -ele $MidEle  material 1 fiber 5 strain;
recorder Element -file ShellData/EleForceSec5Temp.out -time -ele $MidEle   material 1 fiber 5 TempAndElong;

recorder Element -file ShellData/EleForceSec6sigma.out -time -ele $MidEle  material 1 fiber 6 stress;	
recorder Element -file ShellData/EleForceSec6Eps.out -time -ele $MidEle  material 1 fiber 6 strain;
recorder Element -file ShellData/EleForceSec6Temp.out -time -ele $MidEle   material 1 fiber 6 TempAndElong;


recorder Element -file ShellData/EleForceSec7sigma.out -time -ele $MidEle  material 1 fiber 7 stress;	
recorder Element -file ShellData/EleForceSec7Eps.out -time -ele $MidEle  material 1 fiber 7 strain;
recorder Element -file ShellData/EleForceSec7Temp.out -time -ele $MidEle   material 1 fiber 7 TempAndElong;

recorder Element -file ShellData/EleForceSec8sigma.out -time -ele $MidEle  material 1 fiber 8 stress;	
recorder Element -file ShellData/EleForceSec8Eps.out -time -ele $MidEle  material 1 fiber 8 strain;
recorder Element -file ShellData/EleForceSec8Temp.out -time -ele $MidEle   material 1 fiber 8 TempAndElong;

recorder Element -file ShellData/EleForceSec9sigma.out -time -ele $MidEle  material 1 fiber 9 stress;	
recorder Element -file ShellData/EleForceSec9Eps.out -time -ele $MidEle  material 1 fiber 9 strain;
recorder Element -file ShellData/EleForceSec9Temp.out -time -ele $MidEle   material 1 fiber 9 TempAndElong;

recorder Element -file ShellData/EleForceSec10sigma.out -time -ele $MidEle  material 1 fiber 10 stress;	
recorder Element -file ShellData/EleForceSec10Eps.out -time -ele $MidEle  material 1 fiber 10 strain;
recorder Element -file ShellData/EleForceSec10Temp.out -time -ele $MidEle   material 1 fiber 10 TempAndElong;

recorder Element -file ShellData/EleForceSec11sigma.out -time -ele $MidEle  material 1 fiber 11 stress;	
recorder Element -file ShellData/EleForceSec11Eps.out -time -ele $MidEle  material 1 fiber 11 strain;
recorder Element -file ShellData/EleForceSec11Temp.out -time -ele $MidEle   material 1 fiber 11 TempAndElong;

recorder Element -file ShellData/EleForceSec12sigma.out -time -ele $MidEle  material 1 fiber 12 stress;	
recorder Element -file ShellData/EleForceSec12Eps.out -time -ele $MidEle  material 1 fiber 12 strain;
recorder Element -file ShellData/EleForceSec12Temp.out -time -ele $MidEle   material 1 fiber 12 TempAndElong;

recorder Element -file ShellData/EleForceSec13sigma.out -time -ele $MidEle  material 1 fiber 13 stress;	
recorder Element -file ShellData/EleForceSec13Eps.out -time -ele $MidEle  material 1 fiber 13 strain;
recorder Element -file ShellData/EleForceSec13Temp.out -time -ele $MidEle   material 1 fiber 13 TempAndElong;

recorder Element -file ShellData/EleForceSec14sigma.out -time -ele $MidEle  material 1 fiber 14 stress;	
recorder Element -file ShellData/EleForceSec14Eps.out -time -ele $MidEle  material 1 fiber 14 strain;
recorder Element -file ShellData/EleForceSec14Temp.out -time -ele $MidEle   material 1 fiber 14 TempAndElong;

recorder Element -file ShellData/EleForceSec1Lefsigma.out -time -ele $EndEle  material 1 fiber 1 stress;	
recorder Element -file ShellData/EleForceSec1LeftEps.out -time -ele $EndEle  material 1 fiber 1 strain;
recorder Element -file ShellData/EleForceSec1LeftTemp.out -time -ele $EndEle   material 1 fiber 1 TempAndElong;
#print domain.out


if {$ANALYSIS == "HasPoint"} {
#set UDLP [expr -$UDL*$slabB*$slabL/$nx/$ny];

pattern Plain 1 Linear {
   # set NumNodes [expr ($nx+1)*($ny+1)]
	#for {set nodeID 1} {$nodeID<=$NumNodes} {incr nodeID} {
		#load $nodeID 0 0 $UDLP 0 0 0 ;
		#}
    set Load 200;	
	for {set ID 0} {$ID<=$ny} {incr ID} {
		set nodeID [expr ($ID+1)*100+($nx+1)];
		load $nodeID  0 0 [expr $Load/($ny+1)] 0 0 0 ;
    }
		
}
puts "Point";



constraints Plain;
numberer Plain;
system BandGeneral;
test NormDispIncr 1e-3 600 2;
algorithm Newton;
integrator LoadControl 0.5;	
analysis Static;			
analyze 10;
loadConst -time 0.0
}


if {$TANALYSIS == "HasThermo"} {

puts "Thermal action to slab"

#set NumEles [expr $nx*$ny];
#set MidNodeTag [expr (($ny/2+1)*100+($nx/2+1))];
#set EndNodeTag [expr ($ny+1)*100+($nx+1)]

#set minusHalfD [expr -$slabT/2];
#set HalfD [expr $slabT/2];

pattern Plain 3 Linear {
	for {set i 1} {$i <= $nx} {incr i 1} {
	   for {set j 1} {$j <= $ny} {incr j 1} {
	       if {$j == 3 || $j == 4 || $j == 5 || $j ==6} {
	          set sec "rib"
	        } else {
              set sec "flat"
	        }
	        if {$sec=="rib"} { 	  
		       eleLoad  -ele [expr $j*100+$i]  -type -shellThermal  1000    [expr -0.5*$slabT-$offset]  100  [expr 0.5*$slabT-$offset]
	        }  elseif {$sec=="flat"} { 
      	       eleLoad  -ele [expr $j*100+$i]  -type -shellThermal  1000    [expr -0.5*$slabF]  100  [expr 0.5*$slabF]
            }
        }
	}
}		

constraints Plain;
numberer Plain;
system BandGeneral;
#test NormUnbalance 1.0e-4 10 1;
test NormDispIncr 1e-3  600 1;
algorithm Newton;
integrator LoadControl 0.01;	
analysis Static;			
analyze 60;

}

wipe;
