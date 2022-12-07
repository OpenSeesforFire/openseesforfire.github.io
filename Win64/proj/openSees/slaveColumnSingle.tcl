


wipe;

# ------------------------------
# Start of model generation
# ------------------------------
# create ModelBuilder (with two-dimensions and 2 DOF/node)
model BasicBuilder -ndm 3 -ndf 6
# Load OpenFresco package
# -----------------------
# (make sure all dlls are in the same folder as openSees.exe)
#loadPackage OpenFresco

# Define geometry for model Nodes of 1 flange



# -------------------------
# node $tag $xCrd $yCrd
# Define materials

source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl

#source LibUnits.tcl;			# define units
#source DisplayPlane.tcl;		# procedure for displaying a plane in model
#source DisplayModel3D.tcl;	


set K 1.3e11;
set G 0.75e11;
set sig0 455e6;
set sigInf 555e6;

set delta 1;
set H 0;

#set delta2 0;
set mDen  7800;


#nDMaterial J2Thermal 4 3 $K $G $sig0 $sigInf $delta $H 0;

# Another option for J2 plasticity
#nDMaterial  J2PlaneStressThermal $tag $typeTag $nu $fy $fu $hardening factor
nDMaterial  J2PlaneStressThermal 4 3 2e11 0.3 4.68e8 6.00e8 0.1 2e8;
nDMaterial   PlateFromPlaneStressThermal    2   4   20e10;


    # // EN 1992&1993
    # //typeTag:3   EC3 Structural Steel
    # //typeTag:21  EC2 Reinforcing Steel EC2 NHotRolled
    # //typeTag:22  EC2 Reinforcing Steel EC2 NCold formed
    # //typeTag:23  EC2 Reinforcing Steel EC2 X


#set sigInf 44960;

#nDMaterial DruckerPragerThermal 3 2.4e10 1.3e10 7.6e6 $rho $rhoBar $Kinf $K0 1 1 $h 1 $mDen1;
#nDMaterial ElasticIsotropic3DThermal $matTag 2.9e7 0.3 0 1.4e-5;
nDMaterial CDPPlaneStressThermal 3 3.15e10 0.15 4.4e6 4.0e7 9.3e5  3.16e8;
#nDMaterial ElasticIsotropic3DThermal 3 2.9e10 0.3 0 1.4e-5;

nDMaterial   PlateFromPlaneStressThermal    1 3   10e8;

#nDMaterial PlateFiberThermal 2 4;
#nDMaterial PlateFiberThermal 1 3;

#uniaxialMaterial Steel01Thermal 8 475 2.1e11 0.01;
# section: $secTag $matTag $thickness
uniaxialMaterial SteelECThermal 8 EC2NH 4.60e8 2e11;
#uniaxialMaterial ElasticThermal 1 2e11 1.2e-5;
nDMaterial PlateRebarThermal 9 8 0;  # Transverse steel.
nDMaterial PlateRebarThermal 10 8 90; #longitudinal steel.

#section LayeredShellThermal  5  8  2 0.001  10  0.00015 9  0.00015   1 0.009 1  0.0097  1  0.02 1  0.01 1  0.01     ;
#section LayeredShellThermal  5  7 1  0.02 10  0.0005 9  0.0005 1  0.009 1  0.01  1  0.01  1  0.006    ;

#section PlateFiberThermal 5 2 0.12;

section LayeredShellThermal  3  4  2  0.004  2  0.004 2  0.004 2  0.004   ;
section LayeredShellThermal  2  4  2  0.003 2   0.003 2  0.003 2  0.003;
section LayeredShellThermal  6  4  2  0.0035 2  0.0035  2 0.0035  2 0.0035;
section LayeredShellThermal  7  5  2  0.00472  2  0.00472  2  0.00472  2  0.00472  2  0.00472   ;
puts "Secions defined"

set colW 0.4;
set colhw 0.4;
set colH 1.0;
set nx 8;
set ny 20;
set elex [expr $colW/nx]
set eley [expr $colW/nx]
#For Lower Flange
for {set i 1} {$i <= $nx} {incr i 1} {
	for {set j 1} {$j <= $ny} {incr j 1} {
			set x [expr $i*$elemx];
			set y [expr $i*$elemx];
	
	}
	
	}
node	1	0	0	0






puts "Node defined"










#fix	101	1	1	1	0	0	0
#fix	125	0	1	1	0	0	0

#fix	201	1	1	1	0	0	0
#fix	225	0	1	1	0	0	0
#fix	176	1	1	1	0	0	0
#fix	200	0	1	1	0	0	0


	
	
fix	1	1	1	1	1 1 1
fix	26	1	1	1	1 1 1
fix	51	1	1	1	1 1 1
fix	76	1	1	1	1 1 1
fix	101	1	1	1	1 1 1
fix	126	1	1	1	1 1 1

fix	151	1	1	1	1	1	1

fix	176	1	1	1	1	1	1

fix	201	1	1	1	1	1	1





#fix	1	1	1	1	0 0 0
#fix	25	1	1	1	0 0 0
#fix	51	1	1	1	0 0 0
#fix	75	1	1	1	0 0 0
#fix	101	1	1	1	0 0 0
#fix	125	1	1	1	0 0 0
#fix	26	1	1	1	0 0 0
#fix	50	1	1	1	0 0 0
#fix	76	1	1	1	0 0 0
#fix	100	1	1	1	0 0 0


#fix	126	1	1	1	0	0	0
#fix	150	1	1	1	0	0	0
#fix	151	1	1	1	0	0	0
#fix	175	1	1	1	0	0	0
#fix	176	1	1	1	0	0	0
#fix	200	1	1	1	0	0	0
#fix	201	1	1	1	0	0	0
#fix	225	1	1	1	0	0	0










#
 
#define element for traslational&rotational spring

element ShellMITC4Thermal	1	1	2	27	26	2
element ShellMITC4Thermal	2	2	3	28	27	2
element ShellMITC4Thermal	3	3	4	29	28	2
element ShellMITC4Thermal	4	4	5	30	29	2
element ShellMITC4Thermal	5	5	6	31	30	2
element ShellMITC4Thermal	6	6	7	32	31	2
element ShellMITC4Thermal	7	7	8	33	32	2
element ShellMITC4Thermal	8	8	9	34	33	2
element ShellMITC4Thermal	9	9	10	35	34	2
element ShellMITC4Thermal	10	10	11	36	35	2
element ShellMITC4Thermal	11	11	12	37	36	2
element ShellMITC4Thermal	12	12	13	38	37	2
element ShellMITC4Thermal	13	13	14	39	38	2
element ShellMITC4Thermal	14	14	15	40	39	2
element ShellMITC4Thermal	15	15	16	41	40	2
element ShellMITC4Thermal	16	16	17	42	41	2
element ShellMITC4Thermal	17	17	18	43	42	2
element ShellMITC4Thermal	18	18	19	44	43	2
element ShellMITC4Thermal	19	19	20	45	44	2
element ShellMITC4Thermal	20	20	21	46	45	2
element ShellMITC4Thermal	21	21	22	47	46	2
element ShellMITC4Thermal	22	22	23	48	47	2
element ShellMITC4Thermal	23	23	24	49	48	2
element ShellMITC4Thermal	24	24	25	50	49	2
element ShellMITC4Thermal	25	26	27	52	51	2
element ShellMITC4Thermal	26	27	28	53	52	2
element ShellMITC4Thermal	27	28	29	54	53	2
element ShellMITC4Thermal	28	29	30	55	54	2
element ShellMITC4Thermal	29	30	31	56	55	2
element ShellMITC4Thermal	30	31	32	57	56	2
element ShellMITC4Thermal	31	32	33	58	57	2
element ShellMITC4Thermal	32	33	34	59	58	2
element ShellMITC4Thermal	33	34	35	60	59	2
element ShellMITC4Thermal	34	35	36	61	60	2
element ShellMITC4Thermal	35	36	37	62	61	2
element ShellMITC4Thermal	36	37	38	63	62	2
element ShellMITC4Thermal	37	38	39	64	63	2
element ShellMITC4Thermal	38	39	40	65	64	2
element ShellMITC4Thermal	39	40	41	66	65	2
element ShellMITC4Thermal	40	41	42	67	66	2
element ShellMITC4Thermal	41	42	43	68	67	2
element ShellMITC4Thermal	42	43	44	69	68	2
element ShellMITC4Thermal	43	44	45	70	69	2
element ShellMITC4Thermal	44	45	46	71	70	2
element ShellMITC4Thermal	45	46	47	72	71	2
element ShellMITC4Thermal	46	47	48	73	72	2
element ShellMITC4Thermal	47	48	49	74	73	2
element ShellMITC4Thermal	48	49	50	75	74	2
element ShellMITC4Thermal	49	76	77	2	1	2
element ShellMITC4Thermal	50	77	78	3	2	2
element ShellMITC4Thermal	51	78	79	4	3	2
element ShellMITC4Thermal	52	79	80	5	4	2
element ShellMITC4Thermal	53	80	81	6	5	2
element ShellMITC4Thermal	54	81	82	7	6	2
element ShellMITC4Thermal	55	82	83	8	7	2
element ShellMITC4Thermal	56	83	84	9	8	2
element ShellMITC4Thermal	57	84	85	10	9	2
element ShellMITC4Thermal	58	85	86	11	10	2
element ShellMITC4Thermal	59	86	87	12	11	2
element ShellMITC4Thermal	60	87	88	13	12	2
element ShellMITC4Thermal	61	88	89	14	13	2
element ShellMITC4Thermal	62	89	90	15	14	2
element ShellMITC4Thermal	63	90	91	16	15	2
element ShellMITC4Thermal	64	91	92	17	16	2
element ShellMITC4Thermal	65	92	93	18	17	2
element ShellMITC4Thermal	66	93	94	19	18	2
element ShellMITC4Thermal	67	94	95	20	19	2
element ShellMITC4Thermal	68	95	96	21	20	2
element ShellMITC4Thermal	69	96	97	22	21	2
element ShellMITC4Thermal	70	97	98	23	22	2
element ShellMITC4Thermal	71	98	99	24	23	2
element ShellMITC4Thermal	72	99	100	25	24	2
element ShellMITC4Thermal	73	101	102	77	76	2
element ShellMITC4Thermal	74	102	103	78	77	2
element ShellMITC4Thermal	75	103	104	79	78	2
element ShellMITC4Thermal	76	104	105	80	79	2
element ShellMITC4Thermal	77	105	106	81	80	2
element ShellMITC4Thermal	78	106	107	82	81	2
element ShellMITC4Thermal	79	107	108	83	82	2
element ShellMITC4Thermal	80	108	109	84	83	2
element ShellMITC4Thermal	81	109	110	85	84	2
element ShellMITC4Thermal	82	110	111	86	85	2
element ShellMITC4Thermal	83	111	112	87	86	2
element ShellMITC4Thermal	84	112	113	88	87	2
element ShellMITC4Thermal	85	113	114	89	88	2
element ShellMITC4Thermal	86	114	115	90	89	2
element ShellMITC4Thermal	87	115	116	91	90	2
element ShellMITC4Thermal	88	116	117	92	91	2
element ShellMITC4Thermal	89	117	118	93	92	2
element ShellMITC4Thermal	90	118	119	94	93	2
element ShellMITC4Thermal	91	119	120	95	94	2
element ShellMITC4Thermal	92	120	121	96	95	2
element ShellMITC4Thermal	93	121	122	97	96	2
element ShellMITC4Thermal	94	122	123	98	97	2
element ShellMITC4Thermal	95	123	124	99	98	2
element ShellMITC4Thermal	96	124	125	100	99	2

element ShellMITC4Thermal	97	51	52	152	151	3
element ShellMITC4Thermal	98	52	53	153	152	3
element ShellMITC4Thermal	99	53	54	154	153	3
element ShellMITC4Thermal	100	54	55	155	154	3
element ShellMITC4Thermal	101	55	56	156	155	3
element ShellMITC4Thermal	102	56	57	157	156	3
element ShellMITC4Thermal	103	57	58	158	157	3
element ShellMITC4Thermal	104	58	59	159	158	3
element ShellMITC4Thermal	105	59	60	160	159	3
element ShellMITC4Thermal	106	60	61	161	160	3
element ShellMITC4Thermal	107	61	62	162	161	3
element ShellMITC4Thermal	108	62	63	163	162	3
element ShellMITC4Thermal	109	63	64	164	163	3
element ShellMITC4Thermal	110	64	65	165	164	3
element ShellMITC4Thermal	111	65	66	166	165	3
element ShellMITC4Thermal	112	66	67	167	166	3
element ShellMITC4Thermal	113	67	68	168	167	3
element ShellMITC4Thermal	114	68	69	169	168	3
element ShellMITC4Thermal	115	69	70	170	169	3
element ShellMITC4Thermal	116	70	71	171	170	3
element ShellMITC4Thermal	117	71	72	172	171	3
element ShellMITC4Thermal	118	72	73	173	172	3
element ShellMITC4Thermal	119	73	74	174	173	3
element ShellMITC4Thermal	120	74	75	175	174	3
element ShellMITC4Thermal	121	126	127	52	51	3
element ShellMITC4Thermal	122	127	128	53	52	3
element ShellMITC4Thermal	123	128	129	54	53	3
element ShellMITC4Thermal	124	129	130	55	54	3
element ShellMITC4Thermal	125	130	131	56	55	3
element ShellMITC4Thermal	126	131	132	57	56	3
element ShellMITC4Thermal	127	132	133	58	57	3
element ShellMITC4Thermal	128	133	134	59	58	3
element ShellMITC4Thermal	129	134	135	60	59	3
element ShellMITC4Thermal	130	135	136	61	60	3
element ShellMITC4Thermal	131	136	137	62	61	3
element ShellMITC4Thermal	132	137	138	63	62	3
element ShellMITC4Thermal	133	138	139	64	63	3
element ShellMITC4Thermal	134	139	140	65	64	3
element ShellMITC4Thermal	135	140	141	66	65	3
element ShellMITC4Thermal	136	141	142	67	66	3
element ShellMITC4Thermal	137	142	143	68	67	3
element ShellMITC4Thermal	138	143	144	69	68	3
element ShellMITC4Thermal	139	144	145	70	69	3
element ShellMITC4Thermal	140	145	146	71	70	3
element ShellMITC4Thermal	141	146	147	72	71	3
element ShellMITC4Thermal	142	147	148	73	72	3
element ShellMITC4Thermal	143	148	149	74	73	3
element ShellMITC4Thermal	144	149	150	75	74	3
element ShellMITC4Thermal	145	101	102	202	201	3
element ShellMITC4Thermal	146	102	103	203	202	3
element ShellMITC4Thermal	147	103	104	204	203	3
element ShellMITC4Thermal	148	104	105	205	204	3
element ShellMITC4Thermal	149	105	106	206	205	3
element ShellMITC4Thermal	150	106	107	207	206	3
element ShellMITC4Thermal	151	107	108	208	207	3
element ShellMITC4Thermal	152	108	109	209	208	3
element ShellMITC4Thermal	153	109	110	210	209	3
element ShellMITC4Thermal	154	110	111	211	210	3
element ShellMITC4Thermal	155	111	112	212	211	3
element ShellMITC4Thermal	156	112	113	213	212	3
element ShellMITC4Thermal	157	113	114	214	213	3
element ShellMITC4Thermal	158	114	115	215	214	3
element ShellMITC4Thermal	159	115	116	216	215	3
element ShellMITC4Thermal	160	116	117	217	216	3
element ShellMITC4Thermal	161	117	118	218	217	3
element ShellMITC4Thermal	162	118	119	219	218	3
element ShellMITC4Thermal	163	119	120	220	219	3
element ShellMITC4Thermal	164	120	121	221	220	3
element ShellMITC4Thermal	165	121	122	222	221	3
element ShellMITC4Thermal	166	122	123	223	222	3
element ShellMITC4Thermal	167	123	124	224	223	3
element ShellMITC4Thermal	168	124	125	225	224	3
element ShellMITC4Thermal	169	176	177	102	101	3
element ShellMITC4Thermal	170	177	178	103	102	3
element ShellMITC4Thermal	171	178	179	104	103	3
element ShellMITC4Thermal	172	179	180	105	104	3
element ShellMITC4Thermal	173	180	181	106	105	3
element ShellMITC4Thermal	174	181	182	107	106	3
element ShellMITC4Thermal	175	182	183	108	107	3
element ShellMITC4Thermal	176	183	184	109	108	3
element ShellMITC4Thermal	177	184	185	110	109	3
element ShellMITC4Thermal	178	185	186	111	110	3
element ShellMITC4Thermal	179	186	187	112	111	3
element ShellMITC4Thermal	180	187	188	113	112	3
element ShellMITC4Thermal	181	188	189	114	113	3
element ShellMITC4Thermal	182	189	190	115	114	3
element ShellMITC4Thermal	183	190	191	116	115	3
element ShellMITC4Thermal	184	191	192	117	116	3
element ShellMITC4Thermal	185	192	193	118	117	3
element ShellMITC4Thermal	186	193	194	119	118	3
element ShellMITC4Thermal	187	194	195	120	119	3
element ShellMITC4Thermal	188	195	196	121	120	3
element ShellMITC4Thermal	189	196	197	122	121	3
element ShellMITC4Thermal	190	197	198	123	122	3
element ShellMITC4Thermal	191	198	199	124	123	3
element ShellMITC4Thermal	192	199	200	125	124	3

						






puts "ELements defined"


#element  adapter 1000  -node   175 175 175 75 75 75  150 150 150 25 25 25     225 225 225  125 125 125  200 200 200  -dof 1  -dof 2  -dof 3   -dof 1  -dof 2  -dof 3 -dof 1  -dof 2  -dof 3  -dof 1  -dof 2  -dof 3 -dof 1  -dof 2  -dof 3 -dof 1  -dof 2  -dof 3 -dof 1  -dof 2  -dof 3   -stif  3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07	0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3E+07  44000

# Define DISPLAY -------------------------------------------------------------

set  xPixels 1200;	# height of graphical window in pixels
set  yPixels 800;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set ViewScale 0.01;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D  DeformedShape $ViewScale $xLoc1 $yLoc1  $xPixels $yPixels

recorder Element -file ShellData/EleForceSec2Tempm.out -time -ele 153   material 2 fiber 2 TempAndElong;
recorder Element -file ShellData/EleForceSec2stressm.out -time -ele 10  material 2 fiber 2 stress;	

recorder Element -file ShellData/EleForceSec2stress1m.out -time -ele 1  material 2 fiber 2 stress;
recorder Element -file ShellData/EleForceSec2stress20m.out -time -ele 20  material 2 fiber 2 stress;
recorder Element -file ShellData/EleForceSec2stress22m.out -time -ele 22  material 2 fiber 2 stress;

recorder Element -file ShellData/EleForceSec2stress81m.out -time -ele 81  material 2 fiber 2 stress;
recorder Element -file ShellData/EleForceSec2stress90m.out -time -ele 90  material 2 fiber 2 stress;
recorder Element -file ShellData/EleForceSec2stress121m.out -time -ele 122  material 2 fiber 2 stress;
recorder Element -file ShellData/EleForceSec2stress130m.out -time -ele 130  material 2 fiber 2 stress;

recorder Element -file ShellData/EleForceSec2strainm.out -time -ele 10  material 2 fiber 2 strain;

recorder Node -file MidNode138.out -time -node 138  -dof 1 2 3 4 5 6 disp
recorder Node -file MidNode25.out -time -node 25  -dof 1 2  disp
recorder Node -file MidNode125.out -time -node 125  -dof 1 2  disp
recorder Node -file DFreeSlabyT16m.out -time -nodeRange 126  150 -dof 2 disp;
recorder Node -file DFreeSlabyT163m.out -time -nodeRange 176  200 -dof 2 disp;

recorder Node -file DFreeSlabyTwebottom.out -time -nodeRange 101 125 -dof 3 disp;
recorder Node -file DFreeSlabyTwebtop.out -time -nodeRange 76 100 -dof 3 disp;
recorder Node -file DFreeSlabyT13m.out -time -nodeRange 1 25 -dof 3 disp;

recorder Node -file DFreeSlabyT14m.out -time -nodeRange 26  50 -dof 3 disp;

recorder Node -file DFreeSlabyT15m.out -time -nodeRange 51 75 -dof 3 disp;


puts "To start loading"

set UDLP 43;

pattern Plain 1 Linear {
   
	

load	23	0	-$UDLP	0	0	0	0	;

  
}
constraints Plain;
numberer Plain;
system SparseSYM;
test NormDispIncr 1e-4  10 1;
algorithm Newton;
integrator LoadControl 0.1;	
analysis Static;			
analyze 10;
loadConst -time 0.0

#pattern Plain 2 Linear {

	#eleLoad -range 1 96 -type -shellThermal 700 [expr -0.0032] 700 [expr 0.0032];
    #eleLoad -range 97  144 -type -shellThermal 650 [expr -0.0044] 650 [expr 0.0044];
	#eleLoad -range 145  192 -type -shellThermal 800 [expr -0.0044] 800 [expr 0.0044];
	#eleLoad -range 193  288 -type -shellThermal 450 [expr -0.03] 20 [expr 0.03];
#};
pattern Plain 3 Linear {

for {set level 1} {$level <=96} {incr level 1} { 
set eleID $level
eleLoad -ele $eleID -type -shellThermal  600 -0.006 600  0.006;

}

for {set level 97} {$level <=192} {incr level 1} { 
set eleID $level
eleLoad -ele $eleID -type -shellThermal 600 -0.008 600  0.008;

}

}


# Start of analysis generation
# ------------------------------
# create the system of equations
system SparseSYM

# create the DOF numberer
numberer Plain

# create the constraint handler
constraints Transformation

# create the convergence test
test NormDispIncr 1e-3 50 1;

#test NormUnbalance 1.0e-12 25
#test EnergyIncr 1.0e-12 25

# create the integration scheme
integrator LoadControl 0.1

# create the solution algorithm
algorithm Newton

# create the analysis object 
analysis Static
# ------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------
analyze 100
# --------------------------------
# End of analysis
# --------------------------------








