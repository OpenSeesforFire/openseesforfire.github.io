# Recorders
# recorder Node -file output/nodeU3.out -time -nodeRange 100001 [expr 100001 + $num_slab_nodes] -dof 3 disp
recorder Node -file output/midBeamDE.out -time -node 7000121 -dof 3 disp
# recorder Node -file output/midBeam200U3.out -time -node 266 -dof 3 disp
# recorder Node -file output/midBeam300U3.out -time -node 366 -dof 3 disp
# recorder Node -file output/midBeam1300U3.out -time -node 1311 -dof 3 disp
# recorder Node -file output/midBeam1200U3.out -time -node 1211 -dof 3 disp

# recorder Node -file output/Beam100U3.out -time -nodeRange 150 181 -dof 3 disp
# recorder Node -file output/Beam200U3.out -time -nodeRange 250 281 -dof 3 disp
# recorder Node -file output/Beam300U3.out -time -nodeRange 350 381 -dof 3 disp
# recorder Node -file output/Beam1300U3.out -time -nodeRange 1301 1321 -dof 3 disp
# recorder Node -file output/Beam1200U3.out -time -nodeRange 1201 1221 -dof 3 disp

# recorder Node -file output/columnA1x.out -time -node 3309 -dof 1 disp
# recorder Node -file output/columnA1y.out -time -node 3305 -dof 2 disp

# recorder Node -file output/columnB1x.out -time -node 2713 -dof 1 disp
# recorder Node -file output/columnB1y.out -time -node 2709 -dof 2 disp

# recorder Node -file output/columnA2.out -time -node 3409 -dof 1 2 disp

# recorder Node -file output/columnB2.out -time -node 2821 -dof 1 2 disp

# Display
set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 10;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

# apply UDL Load of 5810 N/m2
set UDL -5810
set f3r [expr $UDL*$s_elemx1*$s_elemy]
set f3f [expr $UDL*$s_elemx2*$s_elemy]
#return 
pattern Plain 1 Linear {
	# for {set j 0} {$j <= $s_ny} {incr j 1} {
	   # for {set i 0} {$i <= $s_nx} {incr i 1} {
	     # if {$i < 60 && $j > 80 && $j < 105} {
			# # opening
		 # } elseif {$i > 247 && $j > 98 && $j < 150} {
		    # # openging
		 # } elseif {$i > 49 && $i < 60 && $j > 105 && $j < 150} {
			# # opening
		 # } else {
           # if {$i==0 || $i==$s_nx} {
	         # set f3 [expr $f3r*0.5]
	       # } elseif {$i%4==1 || $i%4==3}  {
             # set f3 [expr 0.5*($f3r+$f3f)]
           # } elseif {$i%4==2 } {
	         # set f3  $f3f
	       # } else {
	         # set f3 $f3r
	       # }		 
	      # set nodeID [expr int(($j+1)*1000+$i+1)] 
	      # #load nodeTag	f1 	f2 	f3 	 f11  f22  f33
		  # load  $nodeID 	0   0   $f3  0 	  0    0
        # }
    # }
# }
	load	181121	0.0	0.0	-2305.0	0.0	0.0	0.0	
	# load	2117	0.0	0.0	-235305.0	0.0	0.0	0.0
	# load	2217	0.0	0.0	-627480.0	0.0	0.0	0.0
	# load	2317	0.0	0.0	-941220.0	0.0	0.0	0.0
	# load	2417	0.0	0.0	-470610.0	0.0	0.0	0.0
	# load	2517	0.0	0.0	-862785.0	0.0	0.0	0.0
	# load	2617	0.0	0.0	-627480.0	0.0	0.0	0.0
	# load	2717	0.0	0.0	-627480.0	0.0	0.0	0.0
	# load	2825	0.0	0.0	-1568700.0	0.0	0.0	0.0
	# load	2925	0.0	0.0	-1568700.0	0.0	0.0	0.0
	# load	3017	0.0	0.0	-627480.0	0.0	0.0	0.0
	# load	3117	0.0	0.0	-392175.0	0.0	0.0	0.0
	# load	3217	0.0	0.0	-522900.0	0.0	0.0	0.0
	# load	3317	0.0	0.0	-313740.0	0.0	0.0	0.0
	# load	3417	0.0	0.0	-174300.0	0.0	0.0	0.0
	# load	3517	0.0	0.0	-104580.0	0.0	0.0	0.0
	# load	3617	0.0	0.0	-122010.0	0.0	0.0	0.0
	# load	3717	0.0	0.0	-313740.0	0.0	0.0	0.0
}
set ok 0
set stepSize 0.1
set currentTime 0.0
set tFinal 1

constraints Transformation     		
numberer RCM			
system SparseSYM	
test NormDispIncr 1e-3 500 1	
#algorithm Newton
algorithm Newton;				
integrator LoadControl $stepSize				
analysis Static

while {$currentTime < $tFinal && $ok == 0} {
	puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	set ok [analyze 1];
	if {$ok == 0} {
		set currentTime [getTime]
		puts "Analysis step successful. Current time is $currentTime"
	} else {
		puts "Analysis failed at time = $currentTime"
		return
	}
}
loadConst -time 0.0

########################
# Ambient Solver END
########################


# set tFinal 144;			# full time of the applied heating

# pattern Plain 2 Linear {
		
		# eleLoad -range 3301 3306  -type 	   -beamThermal -source "temps/ColF1B.dat" -127.2668 127.2668 -64.0751 64.0751
		# eleLoad -range 2701 2706  -type 	   -beamThermal -source "temps/ColE1B.dat" -157.0 157.0 -78.0 78.0
		# eleLoad -range 3401 3406  -type 	   -beamThermal -source "temps/ColF2B.dat" -157.0 157.0 -78.0 78.0
		# eleLoad -range 2801 2806  -type 	   -beamThermal -source "temps/ColE2B.dat" -157.0 157.0 -78.0 78.0
		
		# eleLoad -range 3307 3310  -type 	   -beamThermal -source "temps/ColF1A.dat" -128.0 128.0 -65.0 65.0
		# eleLoad -range 2707 2710  -type 	   -beamThermal -source "temps/ColE1A.dat" -157.0 157.0 -78.0 78.0
		# eleLoad -range 3407 3410  -type 	   -beamThermal -source "temps/ColF2A.dat" -157.0 157.0 -78.0 78.0
		# eleLoad -range 2807 2818  -type 	   -beamThermal -source "temps/ColE2A.dat" -157.0 157.0 -78.0 78.0
	
	
	
		# eleLoad -range 145 148 -type 	   -beamThermal -source "temps/Beam1.dat" -175.5831 175.5831 -42.8751 42.8751
		# eleLoad -range 150 180 -type 	   -beamThermal -source "temps/Beam1.dat" -175.5831 175.5831 -42.8751 42.8751
		
		# eleLoad -range 245 248 -type 	   -beamThermal -source "temps/Beam100.dat" -150.1 150.1 -41.251 41.251
		# eleLoad -range 250 280 -type 	   -beamThermal -source "temps/Beam100.dat" -150.1 150.1 -41.251 41.251
		
		# eleLoad -range 345 348 -type 	   -beamThermal -source "temps/Beam200.dat" -150.1 150.1 -41.251 41.251
		# eleLoad -range 350 380 -type 	   -beamThermal -source "temps/Beam200.dat" -150.1 150.1 -41.251 41.251
		
		
		# eleLoad -range 1201 1220 -type 	   -beamThermal -source "temps/Beam400.dat" -175.5831 175.5831 -42.8751 42.8751
		# eleLoad -range 1222 1225 -type 	   -beamThermal -source "temps/Beam400-II.dat" -298.831 298.831 -56.91 56.91
		
		# eleLoad -range 1301 1320 -type 	   -beamThermal -source "temps/Beam300.dat" -175.5831 175.5831 -42.8751 42.8751
		# eleLoad -range 1322 1325 -type 	   -beamThermal -source "temps/Beam300.dat" -175.5831 175.5831 -42.8751 42.8751
		
		# eleLoad -range 1501 1507 -type 	   -beamThermal -source "temps/Beam500.dat" -175.5831 175.5831 -42.8751 42.8751
		
		# for {set i 47} {$i <= 81} {incr i} {
			# if {$i <= 50} {
				# if {!fmod($i,2)} {; #if even
					# set sec "rib"
				# } else {; #if odd
					# set sec "slab"
				# }
			# } elseif {$i >= 51} {
				# if {!fmod($i,2)} {  #if even
					# set sec "slab"
				# } else { #if odd
					# set sec "rib"
				# }
			# }
			# if {$sec == "rib"} {
				# eleLoad -range [expr $i*1000 + 1] [expr $i*1000 + 28] -type 	 -shellThermal -source "temps/Ribs.dat" [expr -$halfT2-27.5] [expr $halfT2-27.5]
			# } elseif {$sec == "slab"} {
				# eleLoad -range [expr $i*1000 + 1] [expr $i*1000 + 28] -type 	 -shellThermal -source "temps/Slab.dat" [expr -$halfT] [expr $halfT]
			# }
		# }
# }


# #################################
# # Thermo-mechanical Solver
# #################################

# set ok 0
# set stepSize 10
# set currentTime 0.0
# set tFinal 7200;

# constraints Transformation     		
# numberer RCM			
# system SparseSYM			
# test NormDispIncr 1e-2 800 1
# algorithm Newton	
# #algorithm KrylovNewton				
# integrator LoadControl $stepSize				
# analysis Static
# # analyze 450
# while {$currentTime < $tFinal && $ok == 0} {
	# puts "Attempting analysis for time [expr $currentTime + $stepSize], with step size = $stepSize"
	# set ok [analyze 1];
	# if {$ok == 0} {
		# set currentTime [getTime]
		# puts "Analysis step successful. Current time is $currentTime"
	# } else {
		# puts "Analysis failed at time = $currentTime"
		
		# return
	# }
# }
# wipe;
		
