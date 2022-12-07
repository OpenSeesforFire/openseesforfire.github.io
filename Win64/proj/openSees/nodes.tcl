# column nodes
set column_base -4
set column_top  4
set column_H 0.1
set c_x1 9.0
set c_x2 5.006
set c_x3 4.0
set c_x4 3.694
set c_y1 6.0
set c_y2 9.0
set c_y3 2.0
# Columns
for {set i 1} {$i <= [expr 1 + abs($column_top - $column_base)/$column_H]} {incr i 1} {
	set z [expr $column_base + ($i - 1)*$column_H]
	node [expr $i + 1500]          0       [expr $c_y1+0.5*$c_y2]   $z
	node [expr $i + 2500]  [expr 0.5*$c_x1]        0                 $z
	node [expr $i + 3500]  [expr 0.5*$c_x1]      $c_y1               $z
	node [expr $i + 4500]  [expr 0.5*$c_x1] [expr $c_y1+0.5*$c_y2]   $z
	node [expr $i + 5500]  [expr 0.5*$c_x1] [expr $c_y1+$c_y2]       $z
	node [expr $i + 6500]  [expr 0.5*$c_x1] [expr 2*$c_y1+$c_y2]     $z
	node [expr $i + 7500]  [expr 1.5*$c_x1]        0                 $z
	node [expr $i + 8500]  [expr 1.5*$c_x1]      $c_y1               $z
	node [expr $i + 9500]  [expr 1.5*$c_x1] [expr $c_y1+$c_y2]       $z
	node [expr $i + 10500] [expr 1.5*$c_x1] [expr 2*$c_y1+$c_y2]     $z
	node [expr $i + 11500] [expr 1.5*$c_x1+$c_x2] [expr $c_y1+$c_y3] $z
	node [expr $i + 12500] [expr 1.5*$c_x1+$c_x2] [expr $c_y1+$c_y2] $z
	node [expr $i + 13500] [expr 2.5*$c_x1]                0         $z
	node [expr $i + 14500] [expr 2.5*$c_x1]              $c_y1       $z
	node [expr $i + 15500] [expr 2.5*$c_x1] [expr $c_y1+0.5*$c_y2]   $z
	node [expr $i + 16500] [expr 2.5*$c_x1] [expr $c_y1+$c_y2]       $z
	node [expr $i + 17500] [expr 2.5*$c_x1] [expr 2*$c_y1+$c_y2]     $z
	#puts "$i + $z"
}
#return
# slab nodes
set s_H 0
set s_elemx1 [expr 0.5*$uw]
set s_elemx2 [expr 0.5*$fw]
set s_nx [expr (2.5*$c_x1)/($uw+$fw)*4]
set s_elemy 0.1
set s_ny [expr ($c_y1*2+$c_y2)/$s_elemy]
for {set j 0} {$j <= $s_ny} {incr j 1} {       
    set y [expr $j*$s_elemy]
	for {set i 0} {$i <= $s_nx} {incr i 1} {
		 if {$i < 60 && $j > 80 && $j < 105} {
			# opening
		 } elseif {$i > 247 && $j > 98 && $j < 150} {
		    # openging
		 } elseif {$i > 49 && $i < 60 && $j > 104 && $j < 150} {
			# opening
		}  else  {
           if {$i%4 == 1} {
	       set x [expr (($i-1)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	       }   elseif {$i%4 == 2} {
	       set x [expr (($i-2)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	       }   elseif  {$i%4==3} {
	       set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	       }   else {
	       set x [expr ($i/4)*(2*$s_elemx1+2*$s_elemx2)]
	       }
           set nodeID_S [expr int(($j+1)*1000+$i+1)]	
	       node $nodeID_S $x $y $s_H
	       #puts "$nodeID_S + $x + $y + $s_H"
        }
    }
}
#return
# Horizionatal beam nodes
set beamH_610 [expr -0.5*$bh_610-$ribT-0.5*$concT]
set beamH_356 [expr -0.5*$bh_356-$ribT-0.5*$concT]
set beamH_305 [expr -0.5*$bh_305-$ribT-0.5*$concT]
set beamH_254 [expr -0.5*$bh_254-$ribT-0.5*$concT]

# CD 
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 0} {$i <= [expr (0.5*$c_x1)/0.3*4]} {incr i 1} {
        if {$j==1 || $j==2 || $j==6} {
	   set beamH1 $beamH_305
	   } else {
	   set beamH1 $beamH_356
	   } 
	    if {$j<=2} {
		set y [expr $j*(0.5*$c_y1)]
		} elseif {$j==3} {
        set y [expr $c_y1+$c_y3] 
		} elseif {$j==4}  {
		set y [expr $c_y1+0.5*$c_y2]  
		} elseif {$j==5}  {
	    set y [expr $c_y1+$c_y2]  
		} elseif {$j==6}  {
		set y [expr 1.5*$c_y1+$c_y2]
		} else {
		set y [expr 2*$c_y1+$c_y2]
        }
		if {$i%4 == 1} {
	    set x [expr (($i-1)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	}   elseif {$i%4 == 2} {
	    set x [expr (($i-2)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	}   elseif  {$i%4==3} {
	    set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	}   else {
	    set x [expr ($i/4)*(2*$s_elemx1+2*$s_elemx2)]
	}
	   set nodeID_HB1 [expr int(($j+1)*1000000+$i+1)]  
	  node $nodeID_HB1 $x $y $beamH1
	  #puts "$nodeID_HB1 + $x + $y+ $beamH1"
}
}
#return
# DE 
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 61} {$i <= [expr (1.5*$c_x1)/0.3*4+1]} {incr i 1} {
        if {$j==0 || $j==7} {
	   set beamH2 $beamH_356
	   } else {
	   set beamH2 $beamH_305
	   } 
	    if {$j<=2} {
		set y [expr $j*(0.5*$c_y1)]
		} elseif {$j>=3 & $j<=5} {
        set y [expr $c_y1+ ($j-2)*($c_y2/3.0)] 
		} else {
		set y [expr $c_y1+$c_y2+($j-5)*(0.5*$c_y1)]
        }
		if {($i-1)%4 == 1} {
	    set x [expr (($i-2)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	}   elseif {($i-1)%4 == 2} {
	    set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	}   elseif  {($i-1)%4==3} {
	    set x [expr (($i-4)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	}   else {
	    set x [expr (($i-1)/4)*(2*$s_elemx1+2*$s_elemx2)]
	}
	   set nodeID_HB2 [expr int(($j+1)*1000000+$i+1)]  
	  node $nodeID_HB2 $x $y $beamH2
	  #puts "$nodeID_HB2 + $x + $y+ $beamH2"
}
}
# EF
for {set j 0} {$j <= 7} {incr j 1} {
    for {set i 182} {$i <= [expr (2.5*$c_x1)/0.3*4+2]} {incr i 1} {
	    if {$j==3} {
		set j 6
		}
        if {$j==0 || $j== 2 || $j==7} {
	    set beamH3 $beamH_356
	    } else {
	    set beamH3 $beamH_305
	    } 
	    if {$j<=2} {
		set y [expr $j*(0.5*$c_y1)]
		}  else {
		set y [expr $c_y1+$c_y2+($j-5)*(0.5*$c_y1)]
        }
		if {($i-2)%4 == 1} {
	    set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	}   elseif {($i-2)%4 == 2} {
	    set x [expr (($i-4)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	}   elseif  {($i-2)%4==3} {
	    set x [expr (($i-5)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	}   else {
	    set x [expr (($i-2)/4)*(2*$s_elemx1+2*$s_elemx2)]
	}
	   set nodeID_HB3 [expr int(($j+1)*1000000+$i+1)]  
	  node $nodeID_HB3 $x $y $beamH3
	  #puts "$nodeID_HB3 + $x + $y+ $beamH3"
}
}
#EF(2-3-left)
for {set j 3} {$j <= 5} {incr j 1} {
    for {set i 182} {$i <= [expr (1.5*$c_x1+4.8)/0.3*4+2+3]} {incr i 1} {
	   set beamH4 $beamH_305
	   set y [expr $c_y1+($j-2)*($c_y2/3.0)]
		if {($i-2)%4 == 1} {
	    set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	}   elseif {($i-2)%4 == 2} {
	    set x [expr (($i-4)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	}   elseif  {($i-2)%4==3} {
	    set x [expr (($i-5)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	}   else {
	    set x [expr (($i-2)/4)*(2*$s_elemx1+2*$s_elemx2)]
	}
	   set nodeID_HB4 [expr int(($j+1)*1000000+$i+1)]  
	  node $nodeID_HB4 $x $y $beamH4
	  #puts "$nodeID_HB4 + $x + $y+ $beamH4"
}
}
#EF(2-3-right)
set c_y4 1.8
set c_y5 5.2
for {set j 3} {$j <= 5} {incr j 1} {
    for {set i 250} {$i <= [expr (2.5*$c_x1)/0.3*4+3]} {incr i 1} {
	   set beamH5 $beamH_356
	   if {$j==3} {
	   set y [expr $c_y1+$c_y3]
	   } elseif { $j==4} {
	   set y [expr $c_y1+$c_y3+$c_y4]
	   } else {
	   set y [expr $c_y1+$c_y3+$c_y4+$c_y5]
	   }
	   if {($i-3)%4 == 1} {
	   set x [expr (($i-4)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1]
	   } elseif {($i-3)%4 == 2} {
	   set x [expr (($i-5)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+$s_elemx2]
	   } elseif  {($i-3)%4==3} {
	   set x [expr (($i-6)/4)*(2*$s_elemx1+2*$s_elemx2)+$s_elemx1+2*$s_elemx2]
	   } else {
	   set x [expr (($i-3)/4)*(2*$s_elemx1+2*$s_elemx2)]
	   }
	   set nodeID_HB5 [expr int(($j+1)*1000000+$i+1)]  
	  node $nodeID_HB5 $x $y $beamH5
	  #puts "$nodeID_HB5 + $x + $y+ $beamH5"
}
}
#return
# Vertical beam nodes
# D&E&F（1-2）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j 0} {$j <= [expr int($c_y1/$s_elemy)]} {incr j 1} {
    set beamV1 $beamH_356
    set y [expr $j*$s_elemy]
    if {$i==2} {
    set i 3
	} 
	if {$i==3} {
	set x [expr (0.5*$c_x1+($i-1)*$c_x1)]
	} else {
	set x [expr (0.5*$c_x1+$i*$c_x1)]
	}
	set nodeID_VB1 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB1 $x $y $beamV1
	#puts "$nodeID_VB1 + $x + $y+ $beamV1"
    }
}
#return
#D&F（2-2/3）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1/$s_elemy)+1)]} {$j <= [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+1]} {incr j 1} {
    set beamV2 $beamH_356
    set y [expr ($j-1)*$s_elemy]
    if {$i==1} {
    set i 3
	} 
	if {$i==3} {
	set x [expr (0.5*$c_x1+($i-1)*$c_x1)]
	} else {
	set x [expr (0.5*$c_x1+$i*$c_x1)]
	}
	set nodeID_VB2 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB2 $x $y $beamV2
	#puts "$nodeID_VB2 + $x + $y+ $beamV2"
    }
}
#return
#D&F（2/3-3）
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+2]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
    set beamV3 $beamH_356
    set y [expr ($j-2)*$s_elemy]
    if {$i==1} {
    set i 3
	} 
	if {$i==3} {
	set x [expr (0.5*$c_x1+($i-1)*$c_x1)]
	} else {
	set x [expr (0.5*$c_x1+$i*$c_x1)]
	}
	set nodeID_VB3 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB3 $x $y $beamV3
	#puts "$nodeID_VB3 + $x + $y+ $beamV3"
    }
}
#return	
#E(2-3)
for {set j [expr int($c_y1/$s_elemy)+1]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+1]} {incr j 1} {
    set i 1
	set x [expr (0.5*$c_x1+$c_x1)]
	set beamV4 $beamH_610
    set y [expr ($j-1)*$s_elemy]
	set nodeID_VB4 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB4 $x $y $beamV4
	#puts "$nodeID_VB4 + $x + $y+ $beamV4"
}
#return
#(E-F)&(2-3)
for {set j [expr int($c_y1/$s_elemy)+1]} {$j <= [expr int(($c_y1+$c_y3)/$s_elemy)+1]} {incr j 1} {
    set i 2
	set x [expr (0.5*$c_x1+$c_x1+$c_x2)]
	set beamV5 $beamH_254
    set y [expr ($j-1)*$s_elemy]
	set nodeID_VB5 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB5 $x $y $beamV5
	#puts "$nodeID_VB5 + $x + $y+ $beamV5"
}
#return
for {set j [expr int(($c_y1+$c_y3)/$s_elemy)+2]} {$j <= [expr (($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
    set i 2
	set x [expr (0.5*$c_x1+$c_x1+$c_x2)]
	set beamV6 $beamH_356
    set y [expr ($j-2)*$s_elemy]
	set nodeID_VB6 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB6 $x $y $beamV6
	#puts "$nodeID_VB6 + $x + $y+ $beamV6"
}
#return
#D&E&F(3-4)
for {set i 0} {$i <= 3} {incr i 1} {
  for {set j [expr int(($c_y1+$c_y2)/$s_elemy)+3]} {$j <= [expr int((2*$c_y1+$c_y2)/$s_elemy)+3]} {incr j 1} {
    set beamV7 $beamH_356
    set y [expr ($j-3)*$s_elemy]
    if {$i==2} {
    set i 3
	} 
	if {$i==3} {
	set x [expr (0.5*$c_x1+($i-1)*$c_x1)]
	} else {
	set x [expr (0.5*$c_x1+$i*$c_x1)]
	}
	set nodeID_VB7 [expr int(($j+1)*1000000+5000+$i+1)]  
	node $nodeID_VB7 $x $y $beamV7
	#puts "$nodeID_VB7 + $x + $y+ $beamV7"
    }
}
#CD(2/3-3)
for {set j [expr int(($c_y1+0.5*$c_y2)/$s_elemy)+2]} {$j <= [expr int(($c_y1+$c_y2)/$s_elemy)+2]} {incr j 1} {
   set i 0
   set x $c_x4
   set beamV8 $beamH_305
   set y [expr ($j-2)*$s_elemy] 
   set nodeID_VB8 [expr int(($j+1)*1000000+5000+$i)]  
   node $nodeID_VB8 $x $y $beamV8 
   #puts "$nodeID_VB8 + $x + $y+ $beamV8"
}   
#return