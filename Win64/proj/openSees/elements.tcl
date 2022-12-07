geomTransf Corotational   11     1     0     0; #Vertical
geomTransf Corotational   22     0    -1     0; #Horizontal
geomTransf Corotational   33    -1     0     0; #Column

#columns
for {set i 1} {$i <= [expr int(($column_top - $column_base)/$column_H)]} {incr i 1} {
	#element dispBeamColumnThermal $eleTag $iNode $jNode $numIntgrPts  $secTag1  $transfTag
	element dispBeamColumnThermal [expr $i + 1500] [expr $i + 1500]  [expr ($i+1) + 1500]   5 $Sec_198 33
	element dispBeamColumnThermal [expr $i + 2500] [expr $i + 2500]  [expr ($i+1) + 2500]   5 $Sec_198 33
	element dispBeamColumnThermal [expr $i + 3500] [expr $i + 3500]  [expr ($i+1) + 3500]   5 $Sec_198 33
	element dispBeamColumnThermal [expr $i + 4500] [expr $i + 4500]  [expr ($i+1) + 4500]   5 $Sec_198 33
	element dispBeamColumnThermal [expr $i + 5500] [expr $i + 5500]  [expr ($i+1) + 5500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 6500] [expr $i + 6500]  [expr ($i+1) + 6500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 7500] [expr $i + 7500]  [expr ($i+1) + 7500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 8500] [expr $i + 8500]  [expr ($i+1) + 8500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 9500] [expr $i + 9500]  [expr ($i+1) + 9500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 10500] [expr $i + 10500]  [expr ($i+1) + 10500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 11500] [expr $i + 11500]  [expr ($i+1) + 11500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 12500] [expr $i + 12500]  [expr ($i+1) + 12500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 13500] [expr $i + 13500]  [expr ($i+1) + 13500]   5 $Sec_89 33
    element dispBeamColumnThermal [expr $i + 14500] [expr $i + 14500]  [expr ($i+1) + 14500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 15500] [expr $i + 15500]  [expr ($i+1) + 15500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 16500] [expr $i + 16500]  [expr ($i+1) + 16500]   5 $Sec_198 33
    element dispBeamColumnThermal [expr $i + 17500] [expr $i + 17500]  [expr ($i+1) + 17500]   5 $Sec_89 33		
	#puts "$nodeID + $x + $y+ $z"
}
# horizontal beams
# CD 
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 1} {$i <= [expr (0.5*$c_x1)/0.3*4]} {incr i 1} {
	if {$j==1 || $j==2 || $j==6} {
	   set beamSec $Sec_305
	   } else {
	   set beamSec $Sec_356
	   }
    set nodeID_H11 [expr int(($j+1)*1000000+$i)]  
	set nodeID_H12 [expr int($nodeID_H11+1)]
	set elemID_H11 $nodeID_H11
	element dispBeamColumnThermal $elemID_H11 $nodeID_H11 $nodeID_H12  5  $beamSec  22
	  #puts "$elemID_H11 + $nodeID_H11 + $nodeID_H12+ $beamSec"
    }
}
#return
# DE 
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 62} {$i <= [expr (1.5*$c_x1)/0.3*4+1]} {incr i 1} {
        if {$j==0 || $j==7} {
	   set beamSec $Sec_356
	   } else {
	   set beamSec $Sec_305
	   } 
    set nodeID_H21 [expr int(($j+1)*1000000+$i)]  
	set nodeID_H22 [expr int($nodeID_H21+1)]
	set elemID_H21 $nodeID_H21
	element dispBeamColumnThermal $elemID_H21 $nodeID_H21 $nodeID_H22  5  $beamSec  22
	#puts "$elemID_H21 + $nodeID_H21 + $nodeID_H22+ $beamSec"
    }
}
# EF
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 183} {$i <= [expr (2.5*$c_x1)/0.3*4+2]} {incr i 1} {
	    if {$j==3} {
		set j 6
		}
        if {$j==0 || $j== 2 || $j==7} {
	    set beamSec $Sec_356
	    } else {
	    set beamSec $Sec_305
	    } 
    set nodeID_H31 [expr int(($j+1)*1000000+$i)]  
	set nodeID_H32 [expr int($nodeID_H31+1)]
	set elemID_H31 $nodeID_H31
	element dispBeamColumnThermal $elemID_H31 $nodeID_H31 $nodeID_H32  5  $beamSec  22
	#puts "$elemID_H31 + $nodeID_H31 + $nodeID_H32+ $beamSec"
    }
}
#EF(2-3-left)
for {set j 3} {$j <= 5} {incr j 1} {
    for {set i 183} {$i <= [expr (1.5*$c_x1+4.8)/0.3*4+2+3]} {incr i 1} {
    set beamSec $Sec_305
    set nodeID_H41 [expr int(($j+1)*1000000+$i)]  
	set nodeID_H42 [expr int($nodeID_H41+1)]
	set elemID_H41 $nodeID_H41
	element dispBeamColumnThermal $elemID_H41 $nodeID_H41 $nodeID_H42  5  $beamSec  22
	#puts "$elemID_H41 + $nodeID_H41 + $nodeID_H42+ $beamSec"
   }
}
#EF(2-3-right)
for {set j 3} {$j <= 5} {incr j 1} {
    for {set i 251} {$i <= [expr (2.5*$c_x1)/0.3*4+3]} {incr i 1} {
	set beamSec $Sec_356
    set nodeID_H51 [expr int(($j+1)*1000000+$i)]  
	set nodeID_H52 [expr int($nodeID_H51+1)]
	set elemID_H51 $nodeID_H51
	element dispBeamColumnThermal $elemID_H51 $nodeID_H51 $nodeID_H52  5  $beamSec  22
	#puts "$elemID_H51 + $nodeID_H51 + $nodeID_H52+ $beamSec"
   }
}
#return
# Vertical beams
# D&E&F（1-2）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j 1} {$j <= [expr int($c_y1/$s_elemy)]} {incr j 1} {
	if {$i==2} {
    set i 3
	} 
	set beamSec $Sec_356
	set nodeID_V11 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V12 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V11 $nodeID_V11
	element dispBeamColumnThermal $elemID_V11 $nodeID_V11 $nodeID_V12  5  $beamSec  11
	#puts "$elemID_V11 + $nodeID_V11 + $nodeID_V12+ $beamSec"
    }
}
#return
#D&F（2-2/3）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1/$s_elemy)+2)]} {$j <= [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+1]} {incr j 1} {
    if {$i==1} {
    set i 3
	} 
	set beamSec $Sec_356
	set nodeID_V21 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V22 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V21 $nodeID_V21
	element dispBeamColumnThermal $elemID_V21 $nodeID_V21 $nodeID_V22  5  $beamSec  11
	#puts "$elemID_V21 + $nodeID_V21 + $nodeID_V22+ $beamSec"
    }
}
#return
#D&F（2/3-3）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+3]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
    if {$i==1} {
    set i 3
	} 
	set beamSec $Sec_356
	set nodeID_V31 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V32 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V31 $nodeID_V31
	element dispBeamColumnThermal $elemID_V31 $nodeID_V31 $nodeID_V32  5  $beamSec  11
	#puts "$elemID_V31 + $nodeID_V31 + $nodeID_V32+ $beamSec"
    }
}
#return	
#E(2-3)
for {set j [expr int($c_y1/$s_elemy)+2]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+1]} {incr j 1} {
    set i 1
    set beamSec $Sec_610
	set nodeID_V41 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V42 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V41 $nodeID_V41
	element dispBeamColumnThermal $elemID_V41 $nodeID_V41 $nodeID_V42  5  $beamSec  11
	#puts "$elemID_V41 + $nodeID_V41 + $nodeID_V42+ $beamSec"
    }
#puts "E(2-3)"	
#return
#(E-F)&(2-3)
for {set j [expr int($c_y1/$s_elemy)+2]} {$j <= [expr int(($c_y1+$c_y3)/$s_elemy)+1]} {incr j 1} {
    set i 2
	set beamSec $Sec_254
	set nodeID_V51 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V52 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V51 $nodeID_V51
	element dispBeamColumnThermal $elemID_V51 $nodeID_V51 $nodeID_V52  5  $beamSec  11
	#puts "$elemID_V51 + $nodeID_V51 + $nodeID_V52+ $beamSec"
    }
#return
for {set j [expr int(($c_y1+$c_y3)/$s_elemy)+3]} {$j <= [expr (($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
    set i 2
	set beamSec $Sec_356
	set nodeID_V61 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V62 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V61 $nodeID_V61
	element dispBeamColumnThermal $elemID_V61 $nodeID_V61 $nodeID_V62  5  $beamSec  11
	#puts "$elemID_V61 + $nodeID_V61 + $nodeID_V62+ $beamSec"
    }
#return
#D&E&F(3-4)
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1+$c_y2)/$s_elemy)+4]} {$j <= [expr int((2*$c_y1+$c_y2)/$s_elemy)+3]} {incr j 1} {
	if {$i==2} {
    set i 3
	} 
	set beamSec $Sec_356
	set nodeID_V71 [expr int($j*1000000+5000+$i+1)]  
    set nodeID_V72 [expr int(($j+1)*1000000+5000+$i+1)] 
	set elemID_V71 $nodeID_V71
	element dispBeamColumnThermal $elemID_V71 $nodeID_V71 $nodeID_V72  5  $beamSec  11
	#puts "$elemID_V71 + $nodeID_V71 + $nodeID_V72+ $beamSec"
    }
}
#CD(2/3-3)
for {set j [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+3]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
   set i 0
   set beamSec $Sec_305
   set nodeID_V81 [expr int($j*1000000+5000+$i)]  
   set nodeID_V82 [expr int(($j+1)*1000000+5000+$i)] 
   set elemID_V81 $nodeID_V81 
   element dispBeamColumnThermal $elemID_V81 $nodeID_V81 $nodeID_V82  5  $beamSec  11
} 
#puts "OK"
#return
# shell elements
for {set j 1} {$j <= $s_ny} {incr j 1} {
	for {set i 1} {$i <= $s_nx} {incr i 1} {
	    if {$i < 61 && $j >= 81 && $j < 106} {
			# opening
		} elseif {$i >= 248 && $j >= 99 && $j < 151} {
		    # openging
		} elseif {$i >= 50 && $i < 61 && $j >= 106 && $j < 151} {
			# opening
		} else {	   
	        if {$i%4==0 || $i%4==1} {
	        set sec $ribsec
            }  else {
            set sec $flatsec 
            }
	        set node1 [expr int($j*1000+$i)]	
            set node2 [expr $node1+1]			
            set node3 [expr int(($j+1)*1000+$i+1)]			
            set node4 [expr $node3-1]
            set elemID $node1				
	        element ShellNLDKGQThermal	 $elemID $node1 $node2 $node3 $node4 $sec
	        #puts "$elemID + $node1 + $node2 + $node3 + $node4 + $sec"
        }
    }
}
#return
	