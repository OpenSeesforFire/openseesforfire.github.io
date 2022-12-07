#Slab
set slabT 0.13
set concT 0.075
set ribT 0.055
set deckingT 0.9e-3; #thickness of steel deck
set rebarT 0.142e-3; #rebar layer thickness
set fw 0.112
set lw 0.136;#the lower width of the rib
set uw 0.188;#the upper width of the rib
# concrete
set fc 35.0e6
set epsc0 0.0025
set lambda 0.1
set v_c 0.2
set ft [expr 0.1*$fc]
set E_c [expr 1.5*abs($fc/$epsc0)];#[expr sqrt(abs($fpc))*4700.0]; 
set G_c [expr $E_c/(2.0*(1+$v_c))]; #N/mm2
set gt [expr abs($ft/$E_c*$ft*2)];
set gc [expr abs($fc/$E_c*$fc*6)];
# reinforcement
set rfy 460e6
set rEs 2.1e11
# steel decking
set dfy 350e6
set dfu 480e6
set dEs 2.1e11
set vs 0.3

# Concrete  
  nDMaterial CDPPlaneStressThermal 1 $E_c $v_c $ft $fc $gt $gc;
  nDMaterial PlateFromPlaneStressThermal    2   1   10.0e9;
# Reinforcement
  uniaxialMaterial SteelECThermal 3 EC2NC $rfy $rEs;
  #uniaxialMaterial Steel01Thermal 3    $rfy $rEs 0.01;
  nDMaterial PlateRebarThermal 4 	 3 		0;
  nDMaterial PlateRebarThermal 5 	 3 		90;  
# Steel deck 
 #uniaxialMaterial Steel01Thermal $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
  uniaxialMaterial Steel01Thermal 6      $dfy $dEs 0.02;   # Steel deck in the ribs
  nDMaterial PlateRebarThermal 7 	  6 		0;
  nDMaterial PlateRebarThermal 8 	  6 		90; 
  
  nDMaterial  J2PlaneStressThermal 9 21 $dEs $vs $dfy $dfu 0.02 0;
  nDMaterial  PlateFromPlaneStressThermal    10   9   2e10;     

# Beam steel
 #610*229*101UB(S275)
set bh_610 0.6026
set bf_610 0.2276
set tf_610 0.0148
set tw_610 0.0105
set hw_610 0.573
set bfy_610 275e6
#set bfu_610 e6
set bEs_610 2.1e11
set v_610 0.3
set bG_610 [expr $bEs_610/(2.0*(1+$v_610))]; 
set bJ_610 [expr ($bf_610*pow($tf_610,3)*2+$hw_610*pow($tw_610,3))/3.0*1.25]
set mat_610 11  
 uniaxialMaterial Steel01Thermal $mat_610 $bfy_610 $bEs_610 0.01;
 
 #356*171*51UB (S355)
set bh_356 0.355
set bf_356 0.1715
set tf_356 0.0115
set tw_356 0.0074
set hw_356 0.332
set bfy_356 355e6
#set bfu_356 e6
set bEs_356 2.1e11
set v_356 0.3
set bG_356 [expr $bEs_356/(2.0*(1+$v_356))]; 
set bJ_356 [expr ($bf_356*pow($tf_356,3)*2+$hw_356*pow($tw_356,3))/3.0*1.25] 
set mat_356 12  
 uniaxialMaterial Steel01Thermal $mat_356 $bfy_356 $bEs_356 0.01; 
 
 #305*165*40UB (S275)
set bh_305 0.3034
set bf_305 0.165
set tf_305 0.0102
set tw_305 0.006
set hw_305 0.283
set bfy_305 275e6
#set bfu_305 e6
set bEs_305 2.1e11
set v_305 0.3
set bG_305 [expr $bEs_305/(2.0*(1+$v_305))]; 
set bJ_305 [expr ($bf_305*pow($tf_305,3)*2+$hw_305*pow($tw_305,3))/3.0*1.25]
set mat_305 13  
 uniaxialMaterial Steel01Thermal $mat_305 $bfy_305 $bEs_305 0.01; 

 #254*146*31UB (S275)
set bh_254 0.2514
set bf_254 0.1461
set tf_254 0.0086
set tw_254 0.006
set hw_254 0.2342
set bfy_254 275e6
#set bfu_254 e6
set bEs_254 2.1e11
set v_254 0.3
set bG_254 [expr $bEs_254/(2.0*(1+$v_254))]; 
set bJ_254 [expr ($bf_254*pow($tf_254,3)*2+$hw_254*pow($tw_254,3))/3.0*1.25]
set mat_254 14  
 uniaxialMaterial Steel01Thermal $mat_254 $bfy_254 $bEs_254 0.01; 
 
# Column steel
 #305*305*137 (S355)
set ch_137 0.3205
set cf_137 0.3092
set tf_137 0.0217
set tw_137 0.0138
set hw_137 0.2771
set cfy_137 355e6
#set cfu_137 e6
set cEs_137 2.1e11
set v_137 0.3
set cG_137 [expr $cEs_137/(2.0*(1+$v_137))]; 
set cJ_137 [expr ($cf_137*pow($tf_137,3)*2+$hw_137*pow($tw_137,3))/3.0*1.25] 
set mat_137 15  
 uniaxialMaterial Steel01Thermal $mat_137 $cfy_137 $cEs_137 0.01; 
 
 #305*305*198(S355)
set ch_198 0.3399
set cf_198 0.3145
set tf_198 0.0314
set tw_198 0.0191
set hw_198 0.2771
set cfy_198 355e6
#set cfu_198 e6
set cEs_198 2.1e11
set v_198 0.3
set cG_198 [expr $cEs_198/(2.0*(1+$v_198))]; 
set cJ_198 [expr ($cf_198*pow($tf_198,3)*2+$hw_198*pow($tw_198,3))/3.0*1.25] 
set mat_198 16  
 uniaxialMaterial Steel01Thermal $mat_198 $cfy_198 $cEs_198 0.01; 
 
 #254*254*89UC(S355)
set ch_89 0.2603
set cf_89 0.2563
set tf_89 0.0173
set tw_89 0.0103
set hw_89 0.2257
set cfy_89 355e6
#set cfu_89 e6
set cEs_89 2.1e11
set v_89 0.3
set cG_89 [expr $cEs_89/(2.0*(1+$v_89))]; 
set cJ_89 [expr ($cf_89*pow($tf_89,3)*2+$hw_89*pow($tw_89,3))/3.0*1.25]  
set mat_89 17  
 uniaxialMaterial Steel01Thermal $mat_89 $cfy_89 $cEs_89 0.01; 
 
# Flat part section 
set flatsec 1
section  LayeredShellThermal   $flatsec  11	10 $deckingT  2 [expr 0.010-$deckingT]  2 0.010   4 $rebarT   5 $rebarT   2 [expr 0.01-2*$rebarT]   2 0.01   2 0.01   2 0.01  2 0.01  2 0.005;
# Rib section
set ribsec 2
section  LayeredShellThermal   $ribsec  -Rib  $concT $ribT  90.0  21   10 0.000653 -0.09205  2 0.00684 -0.08705  7 0.0001456 -0.0825  2 0.00799 -0.0775 8 0.0001456 -0.0725  2 0.0085 -0.0675  8 0.0001456 -0.0625  2 0.009 -0.0575  8 0.0001456 -0.0525  2 0.0095 -0.0475  2 0.005 -0.04  2 0.01 -0.0325  2 0.01 -0.0225  4 $rebarT -0.017429  5 $rebarT -0.017287  2 0.01 -0.012216  2 0.01 -0.002216 2 0.01 0.007784  2 0.01 0.017784  2 0.01 0.027784 2 0.005 0.035142 

# Sections 
# UB356x171x51
set Sec_356 3
set JG_356 [expr $bG_356*$bJ_356]
section fiberSecThermal $Sec_356 -GJ $JG_356 { 
#fiber      $yLoc        $zLoc      $A      $matTag
fiber	175.5833333e-3		42.875e-3	328.7083333e-6	$mat_356
fiber	171.75e-3			42.875e-3	328.7083333e-6	$mat_356
fiber	167.9166667e-3		42.875e-3	328.7083333e-6	$mat_356
fiber	175.5833333e-3		-42.875e-3	328.7083333e-6	$mat_356
fiber	171.75e-3			-42.875e-3	328.7083333e-6	$mat_356
fiber	167.9166667e-3		-42.875e-3	328.7083333e-6	$mat_356
fiber	138.3333333e-3		0.0		    409.4666667e-6	$mat_356
fiber	83e-3			    0.0		    409.4666667e-6	$mat_356
fiber	27.66666667e-3		0.0		    409.4666667e-6	$mat_356
fiber	-27.66666667e-3	    0.0		    409.4666667e-6	$mat_356
fiber	-83e-3				0.0		    409.4666667e-6	$mat_356
fiber	-138.3333333e-3	    0.0		    409.4666667e-6	$mat_356
fiber	-167.9166667e-3	    -42.875e-3	328.7083333e-6	$mat_356
fiber	-171.75e-3			-42.875e-3	328.7083333e-6	$mat_356
fiber	-175.5833333e-3	    -42.875e-3	328.7083333e-6	$mat_356
fiber	-167.9166667e-3	    42.875e-3	328.7083333e-6	$mat_356
fiber	-171.75e-3			42.875e-3	328.7083333e-6	$mat_356
fiber	-175.5833333e-3	    42.875e-3	328.7083333e-6	$mat_356
};

# UB305x165x40
set Sec_305 4
set JG_305 [expr $bG_305*$bJ_305]
section fiberSecThermal $Sec_305 -GJ $JG_305 { 
#fiber  $yLoc   $zLoc    $A     $matTag
fiber	150.0e-3	  41.25e-3	280.5e-6	$mat_305
fiber	146.6e-3	  41.25e-3	280.5e-6	$mat_305
fiber	143.2e-3	  41.25e-3	280.5e-6	$mat_305
fiber	150.0e-3	  -41.25e-3	280.5e-6	$mat_305
fiber	146.6e-3	  -41.25e-3	280.5e-6	$mat_305
fiber	143.2e-3	  -41.25e-3	280.5e-6	$mat_305
fiber	117.9167e-3 	0.0		283.0e-6	$mat_305
fiber	70.75e-3	    0.0		283.0e-6	$mat_305
fiber	23.5833e-3	    0.0		283.0e-6	$mat_305
fiber	-23.5833e-3	    0.0		283.0e-6	$mat_305
fiber	-70.75e-3	    0.0		283.0e-6	$mat_305
fiber	-117.9167e-3	0.0		283.0e-6	$mat_305
fiber	-143.2e-3	-41.25e-3	280.5e-6	$mat_305
fiber	-146.6e-3	-41.25e-3	280.5e-6	$mat_305
fiber	-150e-3	    -41.25e-3	280.5e-6	$mat_305
fiber	-143.2e-3	41.25e-3	280.5e-6	$mat_305
fiber	-146.6e-3	41.25e-3	280.5e-6	$mat_305
fiber	-150e-3	    41.25e-3	280.5e-6	$mat_305
}

# UB610x229x101
set Sec_610 5
set JG_610 [expr $bG_610*$bJ_610]
section fiberSecThermal $Sec_610 -GJ $JG_610 { 
#fiber       $yLoc   $zLoc          $A  		   $matTag
fiber	298.83e-3	56.9e-3	    561.41e-6	$mat_610
fiber	293.90e-3	56.9e-3	    561.41e-6	$mat_610
fiber	288.97e-3	56.9e-3	    561.41e-6	$mat_610
fiber	298.83e-3	-56.9e-3	561.41e-6	$mat_610
fiber	293.90e-3	-56.9e-3	561.41e-6	$mat_610
fiber	288.97e-3	-56.9e-3	561.41e-6	$mat_610
fiber	238.75e-3	0.0		    1002.75e-6	$mat_610
fiber	143.25e-3	0.0		    1002.75e-6	$mat_610
fiber	47.75e-3	0.0		    1002.75e-6	$mat_610
fiber	-47.75e-3	0.0		    1002.75e-6	$mat_610
fiber	-143.25e-3	0.0		    1002.75e-6	$mat_610
fiber	-238.75e-3  0.0		    1002.75e-6  $mat_610
fiber	-288.97e-3	-56.9e-3	561.41e-6	$mat_610
fiber	-293.90e-3	-56.9e-3	561.41e-6	$mat_610
fiber	-298.83e-3	-56.9e-3	561.41e-6	$mat_610
fiber	-288.97e-3	56.9e-3	    561.41e-6	$mat_610
fiber	-293.90e-3	56.9e-3	    561.41e-6	$mat_610
fiber	-298.83e-3	56.9e-3	    561.41e-6	$mat_610
}

# UB254x146x31
set Sec_254 6
set JG_254 [expr $bG_254*$bJ_254]
section fiberSecThermal $Sec_254 -GJ $JG_254 {
#fiber       $yLoc   $zLoc          $A  $matTag
 fiber 124.2667e-3     -36.525e-3  209.41e-6		  $mat_254
 fiber 121.4000e-3     -36.525e-3  209.41e-6		  $mat_254
 fiber 118.5333e-3     -36.525e-3  209.41e-6		  $mat_254
 fiber 124.2667e-3     36.525e-3   209.41e-6		  $mat_254
 fiber 121.4000e-3     36.525e-3   209.41e-6		  $mat_254
 fiber 118.5333e-3     36.525e-3   209.41e-6		  $mat_254
 fiber 97.5833e-3      0.0         234.2e-6           $mat_254
 fiber 58.55e-3        0.0         234.2e-6           $mat_254
 fiber 19.5167e-3      0.0         234.2e-6           $mat_254
 fiber -19.5167e-3     0.0         234.2e-6           $mat_254
 fiber -58.55e-3       0.0         234.2e-6           $mat_254
 fiber -97.5833e-3     0.0         234.2e-6           $mat_254
 fiber -118.5333e-3    -36.525e-3  209.41e-6		  $mat_254
 fiber -121.4000e-3    -36.525e-3  209.41e-6		  $mat_254
 fiber -124.2667e-3    -36.525e-3  209.41e-6		  $mat_254
 fiber -118.5333e-3     36.525e-3  209.41e-6		  $mat_254
 fiber -121.4000e-3     36.525e-3  209.41e-6		  $mat_254
 fiber -124.2667e-3     36.525e-3  209.41e-6		  $mat_254
}

# COLUMN SECTION
# UC305x305x137
set Sec_137 7
set JG_137 [expr $cG_137*$cJ_137]
section fiberSecThermal $Sec_137 -GJ $JG_137 { 
#fiber       $yLoc      $zLoc      $A  		$matTag
 fiber 156.6333e-3     -77.3e-3 	1118.273e-6  	$mat_137
 fiber 149.4000e-3     -77.3e-3 	1118.273e-6  	$mat_137
 fiber 142.1667e-3     -77.3e-3 	1118.273e-6 	$mat_137
 fiber 156.6333e-3      77.3e-3 	1118.273e-6  	$mat_137
 fiber 149.4000e-3      77.3e-3 	1118.273e-6  	$mat_137
 fiber 142.1667e-3      77.3e-3 	1118.273e-6 	$mat_137
 fiber 115.4583e-3      0.0 	    637.33e-6  	    $mat_137
 fiber  69.275e-3       0.0 	    637.33e-6     	$mat_137
 fiber  23.0916e-3      0.0	        637.33e-6    	$mat_137
 fiber -23.0916e-3      0.0 	    637.33e-6   	$mat_137 
 fiber -69.275e-3       0.0 	    637.33e-6   	$mat_137
 fiber -115.4583e-3     0.0 	    637.33e-6   	$mat_137
 fiber -142.1667e-3     77.3e-3 	1118.273e-6  	$mat_137
 fiber -149.4000e-3     77.3e-3 	1118.273e-6 	$mat_137
 fiber -156.6333e-3     77.3e-3 	1118.273e-6  	$mat_137
 fiber -142.1667e-3    -77.3e-3 	1118.273e-6  	$mat_137
 fiber -149.4000e-3    -77.3e-3 	1118.273e-6 	$mat_137
 fiber -156.6333e-3    -77.3e-3 	1118.273e-6  	$mat_137
}

# UC305x305x198
set Sec_198 8
set JG_198 [expr $cG_198*$cJ_198]
section fiberSecThermal $Sec_198 -GJ $JG_198 {
#fiber       $yLoc   $zLoc          $A  $matTag
 fiber 164.7167e-3      -78.625e-3  1645.9e-6  $mat_198
 fiber 154.2500e-3      -78.625e-3  1645.9e-6  $mat_198
 fiber 143.7833e-3      -78.625e-3  1645.9e-6  $mat_198
 fiber 164.7167e-3      78.625e-3   1645.9e-6  $mat_198
 fiber 154.2500e-3      78.625e-3   1645.9e-6  $mat_198
 fiber 143.7833e-3      78.625e-3   1645.9e-6  $mat_198
 fiber 115.4583e-3      0.0         882.1017e-6  $mat_198
 fiber 69.275e-3        0.0         882.1017e-6  $mat_198
 fiber 23.0917e-3       0.0         882.1017e-6  $mat_198
 fiber -23.0917e-3      0.0         882.1017e-6  $mat_198 
 fiber -69.275e-3       0.0         882.1017e-6  $mat_198
 fiber -115.4583e-3     0.0         882.1017e-6  $mat_198
 fiber -143.7833e-3     78.625e-3   1645.9e-6  $mat_198
 fiber -154.2500e-3     78.625e-3   1645.9e-6  $mat_198
 fiber -164.7167e-3     78.625e-3   1645.9e-6  $mat_198
 fiber -143.7833e-3     -78.625e-3  1645.9e-6  $mat_198
 fiber -154.2500e-3     -78.625e-3  1645.9e-6  $mat_198
 fiber -164.7167e-3     -78.625e-3  1645.9e-6  $mat_198
}

# UC254x254x89
set Sec_89 9
set JG_89 [expr $cJ_89*$cG_89]
section fiberSecThermal $Sec_89 -GJ $JG_89 { 
#fiber       $yLoc   $zLoc          $A  $matTag
 fiber 127.2667e-3     -64.075e-3  739.0e-6     $mat_89
 fiber 121.5000e-3     -64.075e-3  739.0e-6     $mat_89
 fiber 115.7333e-3     -64.075e-3  739.0e-6     $mat_89
 fiber 127.2667e-3     64.075e-3   739.0e-6     $mat_89
 fiber 121.5000e-3     64.075e-3   739.0e-6     $mat_89
 fiber 115.7333e-3     64.075e-3   739.0e-6     $mat_89
 fiber 94.04167e-3     0.0         387.4517e-6  $mat_89
 fiber 56.425e-3       0.0         387.4517e-6  $mat_89
 fiber 18.80833e-3     0.0         387.4517e-6  $mat_89
 fiber -18.80833e-3    0.0         387.4517e-6  $mat_89
 fiber -56.425e-3      0.0         387.4517e-6  $mat_89
 fiber -94.04167e-3    0.0         387.4517e-6  $mat_89 
 fiber -115.7333e-3    64.075e-3   739.0e-6     $mat_89
 fiber -121.5000e-3    64.075e-3   739.0e-6     $mat_89
 fiber -127.2667e-3    64.075e-3   739.0e-6     $mat_89
 fiber -115.7333e-3    -64.075e-3  739.0e-6     $mat_89
 fiber -121.5000e-3    -64.075e-3  739.0e-6     $mat_89
 fiber -127.2667e-3    -64.075e-3  739.0e-6     $mat_89
}

