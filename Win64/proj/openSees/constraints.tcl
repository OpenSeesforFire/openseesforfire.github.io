# Beam (symmetric boundary conditions)
for {set i 1} {$i <= 8} {incr i 1} {
    if {$i == 5} {
	#column
	} else {
	fix [expr $i*1000000 + 1] 1 1 0 1 1 1
    }
}	

# Column bases and tops
for {set i 1} {$i <= 17} {incr i 1} {
	fix [expr $i*1000 + 501] 1 1 1 1 1 1
    fix [expr $i*1000 + 581] 1 1 0 1 1 1 
}

# Column mastering slab's translation and rotation about z (along the vertial line)
equalDOF	1541	106001	1 2 3 6

equalDOF	2541	1061	1 2 3 6
equalDOF	3541	61061	1 2 3 6
equalDOF	4541	106061	1 2 3 6
equalDOF	5541	151061	1 2 3 6
equalDOF	6541	211061	1 2 3 6

equalDOF	7541	1181	1 2 3 6
equalDOF	8541	61181	1 2 3 6
equalDOF	9541	151181	1 2 3 6
equalDOF	10541	211181	1 2 3 6

equalDOF	11541	81248	1 2 3 6
equalDOF	12541	151248	1 2 3 6

equalDOF	13541	1301	1 2 3 6
equalDOF	14541	61301	1 2 3 6
equalDOF	16541	151301	1 2 3 6
equalDOF	17541	211301	1 2 3 6
# puts "ok"
# return
# Column mastering beams' translation and rotations ( along the vertial line)
equalDOF	1538	5000001		1 2 3 4 5 6

equalDOF	2538	1000061		1 2 3 4 5 6
equalDOF	2538	1000062		1 2 3 4 5 6
equalDOF	2538	1005001		1 2 3 4 5 6

equalDOF	3539	3000061		1 2 3 4 5 6
equalDOF	3539	3000062		1 2 3 4 5 6
equalDOF	3538	62005001	1 2 3 4 5 6
equalDOF	3538	61005001	1 2 3 4 5 6

equalDOF	4538	107005001	1 2 3 4 5 6
equalDOF	4538	108005001	1 2 3 4 5 6
equalDOF	4538	5000061		1 2 3 4 5 6

equalDOF	5538	6000061		1 2 3 4 5 6
equalDOF	5538	153005001	1 2 3 4 5 6
equalDOF	5538	154005001	1 2 3 4 5 6
equalDOF	5539	6000062		1 2 3 4 5 6

equalDOF	6538	8000061	    1 2 3 4 5 6
equalDOF	6538	8000062		1 2 3 4 5 6
equalDOF	6538	214005001	1 2 3 4 5 6

equalDOF	7538	1000182		1 2 3 4 5 6
equalDOF	7538	1000183		1 2 3 4 5 6
equalDOF	7538	1005002		1 2 3 4 5 6

equalDOF	8539	3000182		1 2 3 4 5 6
equalDOF	8538	3000183		1 2 3 4 5 6
equalDOF	8538	61005002	1 2 3 4 5 6
equalDOF	8537	62005002	1 2 3 4 5 6

equalDOF	9539	6000182		1 2 3 4 5 6
equalDOF	9539	6000183		1 2 3 4 5 6
equalDOF	9537	152005002	1 2 3 4 5 6
equalDOF	9538	154005002	1 2 3 4 5 6

equalDOF	10538	8000182		1 2 3 4 5 6
equalDOF	10538	8000183		1 2 3 4 5 6
equalDOF	10538	214005002	1 2 3 4 5 6

equalDOF	11539	82005003	1 2 3 4 5 6
equalDOF	11538	4000251		1 2 3 4 5 6
equalDOF	11538	83005003	1 2 3 4 5 6

equalDOF	12539	6000250  	1 2 3 4 5 6
equalDOF	12538	6000251		1 2 3 4 5 6
equalDOF	12538	153005003	1 2 3 4 5 6

equalDOF	13538	1000303		1 2 3 4 5 6
equalDOF	13538	1005004		1 2 3 4 5 6

equalDOF	14538	3000303		1 2 3 4 5 6
equalDOF	14538	61005004	1 2 3 4 5 6
equalDOF	14538	62005004	1 2 3 4 5 6

equalDOF	15538	107005004	1 2 3 4 5 6
equalDOF	15538	108005004	1 2 3 4 5 6

equalDOF	16538	6000304	    1 2 3 4 5 6
equalDOF	16538	153005004	1 2 3 4 5 6
equalDOF	16538	154005004	1 2 3 4 5 6

equalDOF	17538	8000303 	1 2 3 4 5 6
equalDOF	17538	214005004	1 2 3 4 5 6
#puts "ok"
#return
# Horizontal beams to primary beams ( along the horizontal lines)
equalDOF	31005001	2000061 	1 2 3 4 5 6
equalDOF	31005001	2000062 	1 2 3 4 5 6
equalDOF	31005002	2000182 	1 2 3 4 5 6
equalDOF	31005002	2000183 	1 2 3 4 5 6
equalDOF	31005004	2000303 	1 2 3 4 5 6

equalDOF	3000250 	62005003 	1 2 3 4 5 6

equalDOF	82005001 	4000061 	1 2 3 4 5 6
equalDOF	82005004 	4000304 	1 2 3 4 5 6

equalDOF	92005001 	4000062 	1 2 3 4 5 6
equalDOF	92005002 	4000182 	1 2 3 4 5 6
equalDOF	92005002 	4000183 	1 2 3 4 5 6
equalDOF	93005003 	4000250 	1 2 3 4 5 6

equalDOF	101005003 	5000251 	1 2 3 4 5 6
equalDOF	100005004 	5000304 	1 2 3 4 5 6

equalDOF	5000050 	108005000 	1 2 3 4 5 6
equalDOF	6000050 	153005000 	1 2 3 4 5 6

equalDOF	123005001 	5000062 	1 2 3 4 5 6
equalDOF	122005002 	5000182 	1 2 3 4 5 6
equalDOF	122005002 	5000183 	1 2 3 4 5 6
equalDOF	123005003 	5000250 	1 2 3 4 5 6

equalDOF	184005001	7000061 	1 2 3 4 5 6
equalDOF	184005001	7000062 	1 2 3 4 5 6
equalDOF	184005002	7000182 	1 2 3 4 5 6
equalDOF	184005002	7000183 	1 2 3 4 5 6
equalDOF	184005004	7000303 	1 2 3 4 5 6
# puts "OK"
# return
# Slab to beams
# Horizontal beams
# C-D
set k 1
for {set i 1} {$i < 61} {incr i 1} {
  for {set j 1} {$j <= 8} {incr j 1} {
   if {$i%4 == 1} {
       if {$j != 5} {
	      if {$j == 2}  { 
          set k 31	
	      } elseif {$j == 3} {
	      set k 61
	      } elseif {$j == 4} {
	      set k 81 
	      } elseif {$j == 6} {
	      set k 151
	      } elseif {$j == 7} {
	      set k 181
          } elseif {$j == 8} {
          set k 211
          }
          set beamID [expr $j*1000000+$i]
          set slabID [expr $k*1000+$i]
	      rigidLink beam $beamID $slabID
		  #puts "$beamID + $slabID"	
          } elseif {$i>1 && $i<=50 && $j == 5} {
          set beamID [expr $j*1000000+$i]
          set slabID [expr 106000+$i]
	      rigidLink beam $beamID $slabID
		  #puts "$beamID + $slabID"	
          } else {
	      #No link
	      }
    } else {
	#No link
    }	
} 
} 
#return 
#D-E
for {set i 62} {$i < 181} {incr i 1} {
  for {set j 1} {$j <= 8} {incr j 1} {
     if {$i%4 == 1} {
	    if {$j == 1} {
		set k 1
		} elseif {$j == 2} { 
        set k 31	
	    } elseif {$j == 3} {
	    set k 61
	    } elseif {$j == 4} {
	    set k 91 
	    } elseif {$j == 5} {
	    set k 121
	    } elseif {$j == 6} {
	    set k 151
	    } elseif {$j == 7} {
	    set k 181
        } elseif {$j == 8} {
        set k 211
        } 
        set beamID [expr $j*1000000+($i+1)]
        set slabID [expr $k*1000+$i]
        rigidLink beam $beamID $slabID
		 #puts "$beamID + $slabID"		
     } else {
	#No link
	 }
  }
}
#return  
#EF
for {set i 182} {$i< 301} {incr i 1} {
  for {set j 1} {$j <= 8} {incr j 1} {
     if {$i%4 == 1} {
        if {$j <= 3 || $j >= 7} { 
          if {$j == 1} {
		  set k 1	
	      } elseif {$j == 2} {
	      set k 31
	      } elseif {$j == 3} {
	      set k 61 
	      } elseif {$j == 7} {
	      set k 181
          } elseif {$j == 8} {
          set k 211
		  }
          set beamID [expr $j*1000000+($i+2)]
          set slabID [expr $k*1000+$i]
		  rigidLink beam $beamID $slabID
		 #puts "$beamID + $slabID"	
          } else {
		     if {$i < 248} {
		     if {$j == 4} {
             set k 91 
             } elseif {$j == 5} {
             set k 121
             } elseif {$j == 6} {
             set k 151
             }
             set beamID [expr $j*1000000+($i+2)]
             set slabID [expr $k*1000+$i]
		     rigidLink beam $beamID $slabID	
		   # puts "$beamID + $slabID"				 
             } else {
		     if {$j == 4} {
             set k 81 
             } elseif {$j == 5} {
             set k 99
             } elseif {$j == 6} {
             set k 151
             }
             set beamID [expr $j*1000000+($i+3)]
             set slabID [expr $k*1000+$i]           			
		     rigidLink beam $beamID $slabID	
		 #puts "$beamID + $slabID"				 
            }
        }
    } else {
	# No link 
	}  		
}
}	
#return
# slab to vertical beams
# D&E$F (1-2)
for {set j 2} {$j <= 60} {incr j 1} {
    for {set i 1} {$i <= 4} {incr i 1} {
	if {$i == 3} {
	set i 4
	}
	if {$i == 1} {
	set k 61 
	} elseif {$i == 2} {
	set k 181
	} else {
	set k 301
	}
	set beamID [expr $j*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"	
}
}
#return
# D (2-2/3)
for {set j 62} {$j <= 105} {incr j 1} {
    set i 1
	set k 61
	set beamID [expr ($j+1)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}
#return
# F (2/3-3)
for {set j 62} {$j <= 99} {incr j 1} {
    set i 4
	set k 301
	set beamID [expr ($j+1)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}
#return
for {set j 107} {$j <= 150} {incr j 1} {
	set i 1 
	set k 61 
	set beamID [expr ($j+2)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}
#return
# D&E$F (3-4)
for {set j 152} {$j <= 210} {incr j 1} {
    for {set i 1} {$i <= 4} {incr i 1} {
	if {$i == 3} {
	set i 4
	}
	if {$i == 1} {
	set k 61 
	} elseif {$i == 2} {
	set k 181
	} else {
	set k 301
	}
	set beamID [expr ($j+3)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID	
    #puts "$beamID + $slabID"		
}
}
#return
# E (2-3)
for {set j 62} {$j <= 150} {incr j 1} {
	set i 2
	set k 181
	set beamID [expr ($j+1)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID	
    #puts "$beamID + $slabID"		
}

# EF (2-2/3)
for {set j 62} {$j <= 80} {incr j 1} {
	set i 3
	set k 248
	set beamID [expr ($j+1)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}

# EF (2/3-3)
for {set j 82} {$j <= 150} {incr j 1} {
	set i 3
	set k 248
	set beamID [expr ($j+2)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}

# CD (2/3-3)
for {set j 107} {$j <= 150} {incr j 1} {
	set i 0
	set k 50
	set beamID [expr ($j+2)*1000000+5000+$i]
    set slabID [expr $j*1000+$k]           			
    rigidLink beam $beamID $slabID
    #puts "$beamID + $slabID"		
}
#return