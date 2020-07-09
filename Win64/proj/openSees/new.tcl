wipe
# UNITS N mm
model BasicBuilder -ndm 2 -ndf 3;
# create nodes
node	1001	0	0
node	1002	6000	0
node	1003	12000	0
node	1004	18000	0
node	1005	24000	0
node	1006	30000	0
node	2001	0	3000
node	2002	6000	3000
node	2003	12000	3000
node	2004	18000	3000
node	2005	24000	3000
node	2006	30000	3000
node	3001	0	6000
node	3002	6000	6000
node	3003	12000	6000
node	3004	18000	6000
node	3005	24000	6000
node	3006	30000	6000
node	4001	0	9000
node	4002	6000	9000
node	4003	12000	9000
node	4004	18000	9000
node	4005	24000	9000
node	4006	30000	9000
node	5001	0	12000
node	5002	6000	12000
node	5003	12000	12000
node	5004	18000	12000
node	5005	24000	12000
node	5006	30000	12000
node	6001	0	15000
node	6002	6000	15000
node	6003	12000	15000
node	6004	18000	15000
node	6005	24000	15000
node	6006	30000	15000

# inner node 
node	20111	3000	3000
node	20121	9000	3000
node	20131	15000	3000
node	20141	21000	3000
node	20151	27000	3000
			
node	30111	3000	6000
node	30121	9000	6000
node	30131	15000	6000
node	30141	21000	6000
node	30151	27000	6000
			
node	40111	3000	9000
node	40121	9000	9000
node	40131	15000	9000
node	40141	21000	9000
node	40151	27000	9000
			
node	50111	3000	12000
node	50121	9000	12000
node	50131	15000	12000
node	50141	21000	12000
node	50151	27000	12000
			
node	60111	3000	15000
node	60121	9000	15000
node	60131	15000	15000
node	60141	21000	15000
node	60151	27000	15000

# aux. nodes
node	10023	6000	0
node	10034	12000	0
node	10043	18000	0
node	10054	24000	0
node	20011	0	3000
node	20021	6000	3000
node	20022	6000	3000
node	20023	6000	3000
node	20031	12000	3000
node	20032	12000	3000
node	20034	12000	3000
node	20041	18000	3000
node	20042	18000	3000
node	20043	18000	3000
node	20051	24000	3000
node	20052	24000	3000
node	20054	24000	3000
node	20062	30000	3000
node	30011	0	6000
node	30021	6000	6000
node	30022	6000	6000
node	30023	6000	6000
node	30031	12000	6000
node	30032	12000	6000
node	30034	12000	6000
node	30041	18000	6000
node	30042	18000	6000
node	30043	18000	6000
node	30051	24000	6000
node	30052	24000	6000
node	30054	24000	6000
node	30062	30000	6000
node	40011	0	9000
node	40021	6000	9000
node	40022	6000	9000
node	40023	6000	9000
node	40031	12000	9000
node	40032	12000	9000
node	40034	12000	9000
node	40041	18000	9000
node	40042	18000	9000
node	40043	18000	9000
node	40051	24000	9000
node	40052	24000	9000
node	40054	24000	9000
node	40062	30000	9000
node	50011	0	12000
node	50021	6000	12000
node	50022	6000	12000
node	50023	6000	12000
node	50031	12000	12000
node	50032	12000	12000
node	50034	12000	12000
node	50041	18000	12000
node	50042	18000	12000
node	50043	18000	12000
node	50051	24000	12000
node	50052	24000	12000
node	50054	24000	12000
node	50062	30000	12000
node	60011	0	15000
node	60021	6000	15000
node	60022	6000	15000
node	60023	6000	15000
node	60031	12000	15000
node	60032	12000	15000
node	60034	12000	15000
node	60041	18000	15000
node	60042	18000	15000
node	60043	18000	15000
node	60051	24000	15000
node	60052	24000	15000
node	60054	24000	15000
node	60062	30000	15000

# mass 
mass	2001	1.799E+00	1.799E+00	0
mass	2002	3.598E+00	3.598E+00	0
mass	2003	3.598E+00	3.598E+00	0
mass	2004	3.598E+00	3.598E+00	0
mass	2005	3.598E+00	3.598E+00	0
mass	2006	1.799E+00	1.799E+00	0
mass	3001	1.799E+00	1.799E+00	0
mass	3002	3.598E+00	3.598E+00	0
mass	3003	3.598E+00	3.598E+00	0
mass	3004	3.598E+00	3.598E+00	0
mass	3005	3.598E+00	3.598E+00	0
mass	3006	1.799E+00	1.799E+00	0
mass	4001	1.799E+00	1.799E+00	0
mass	4002	3.598E+00	3.598E+00	0
mass	4003	3.598E+00	3.598E+00	0
mass	4004	3.598E+00	3.598E+00	0
mass	4005	3.598E+00	3.598E+00	0
mass	4006	1.799E+00	1.799E+00	0
mass	5001	1.799E+00	1.799E+00	0
mass	5002	3.598E+00	3.598E+00	0
mass	5003	3.598E+00	3.598E+00	0
mass	5004	3.598E+00	3.598E+00	0
mass	5005	3.598E+00	3.598E+00	0
mass	5006	1.799E+00	1.799E+00	0
mass	6001	1.799E+00	1.799E+00	0
mass	6002	3.598E+00	3.598E+00	0
mass	6003	3.598E+00	3.598E+00	0
mass	6004	3.598E+00	3.598E+00	0
mass	6005	3.598E+00	3.598E+00	0
mass	6006	1.799E+00	1.799E+00	0


# Single point constraints -- Boundary Conditions		
fixY 0.0 1 1 0; 

print "done1st"

# equalDOF
equalDOF	1002	10023	1	2
equalDOF	1003	10034	1	2
equalDOF	1004	10043	1	2
equalDOF	1005	10054	1	2
equalDOF	2001	20011	1	2
equalDOF	2002	20021	1	2
equalDOF	2002	20022	1	2
equalDOF	2002	20023	1	2
equalDOF	2003	20031	1	2
equalDOF	2003	20032	1	2
equalDOF	2003	20034	1	2
equalDOF	2004	20041	1	2
equalDOF	2004	20042	1	2
equalDOF	2004	20043	1	2
equalDOF	2005	20051	1	2
equalDOF	2005	20052	1	2
equalDOF	2005	20054	1	2
equalDOF	2006	20062	1	2
equalDOF	3001	30011	1	2
equalDOF	3002	30021	1	2
equalDOF	3002	30022	1	2
equalDOF	3002	30023	1	2
equalDOF	3003	30031	1	2
equalDOF	3003	30032	1	2
equalDOF	3003	30034	1	2
equalDOF	3004	30041	1	2
equalDOF	3004	30042	1	2
equalDOF	3004	30043	1	2
equalDOF	3005	30051	1	2
equalDOF	3005	30052	1	2
equalDOF	3005	30054	1	2
equalDOF	3006	30062	1	2
equalDOF	4001	40011	1	2
equalDOF	4002	40021	1	2
equalDOF	4002	40022	1	2
equalDOF	4002	40023	1	2
equalDOF	4003	40031	1	2
equalDOF	4003	40032	1	2
equalDOF	4003	40034	1	2
equalDOF	4004	40041	1	2
equalDOF	4004	40042	1	2
equalDOF	4004	40043	1	2
equalDOF	4005	40051	1	2
equalDOF	4005	40052	1	2
equalDOF	4005	40054	1	2
equalDOF	4006	40062	1	2
equalDOF	5001	50011	1	2
equalDOF	5002	50021	1	2
equalDOF	5002	50022	1	2
equalDOF	5002	50023	1	2
equalDOF	5003	50031	1	2
equalDOF	5003	50032	1	2
equalDOF	5003	50034	1	2
equalDOF	5004	50041	1	2
equalDOF	5004	50042	1	2
equalDOF	5004	50043	1	2
equalDOF	5005	50051	1	2
equalDOF	5005	50052	1	2
equalDOF	5005	50054	1	2
equalDOF	5006	50062	1	2
equalDOF	6001	60011	1	2
equalDOF	6002	60021	1	2
equalDOF	6002	60022	1	2
equalDOF	6002	60023	1	2
equalDOF	6003	60031	1	2
equalDOF	6003	60032	1	2
equalDOF	6003	60034	1	2
equalDOF	6004	60041	1	2
equalDOF	6004	60042	1	2
equalDOF	6004	60043	1	2
equalDOF	6005	60051	1	2
equalDOF	6005	60052	1	2
equalDOF	6005	60054	1	2
equalDOF	6006	60062	1	2

#equalDOF $MasterNode $SlaveNode $DOFs
equalDOF	2001	2002	1	
equalDOF	2001	2003	1	
equalDOF	2001	2004	1	
equalDOF	2001	2005	1	
equalDOF	2001	2006	1	
				
				
equalDOF	3001	3002	1	
equalDOF	3001	3003	1	
equalDOF	3001	3004	1	
equalDOF	3001	3005	1	
equalDOF	3001	3006	1	
				
				
equalDOF	4001	4002	1	
equalDOF	4001	4003	1	
equalDOF	4001	4004	1	
equalDOF	4001	4005	1	
equalDOF	4001	4006	1	
				
				
equalDOF	5001	5002	1	
equalDOF	5001	5003	1	
equalDOF	5001	5004	1	
equalDOF	5001	5005	1	
equalDOF	5001	5006	1	
				
				
equalDOF	6001	6002	1	
equalDOF	6001	6003	1	
equalDOF	6001	6004	1	
equalDOF	6001	6005	1	
equalDOF	6001	6006	1	

# inner
equalDOF	2001	20111	1
equalDOF	2001	20121	1
equalDOF	2001	20131	1
equalDOF	2001	20141	1
equalDOF	2001	20151	1
			
equalDOF	3001	30111	1
equalDOF	3001	30121	1
equalDOF	3001	30131	1
equalDOF	3001	30141	1
equalDOF	3001	30151	1
			
equalDOF	4001	40111	1
equalDOF	4001	40121	1
equalDOF	4001	40131	1
equalDOF	4001	40141	1
equalDOF	4001	40151	1
			
equalDOF	5001	50111	1
equalDOF	5001	50121	1
equalDOF	5001	50131	1
equalDOF	5001	50141	1
equalDOF	5001	50151	1


# material
set matID 100

uniaxialMaterial Steel01Thermal $matID 235 2.0e5 0.01

# sections 
set W16x26 100
set W16x36 101
set W16x67 200
set U100   300
set U140   400

section fiberSecThermal  $W16x26  {
	patch quad  $matID  5 1   -200 70   -200 -70   -190.7 -70   -190.7 70
	patch quad  $matID  1 5   -190.7 3.2   -190.7 -3.2   190.7 -3.2   190.7 3.2
	patch quad  $matID  5 1   190.7 70   190.7 -70   199 -70   199 70
}

section fiberSecThermal  $W16x36  {
	patch quad  $matID  5 1   -202 88   -202 -89   -190.6 -89   -190.6 88
	patch quad  $matID  1 5   -190.6 3.75   -190.6 -3.75   190.6 -3.75   190.6 3.75
	patch quad  $matID  5 1   190.6 88   190.6 -89   201 -89   201 88
}

section fiberSecThermal  $W16x67  {
	patch quad  $matID  5 1   -208 130   -208 -130   -190.6 -130   -190.6 130
	patch quad  $matID  1 5   -190.6 5   -190.6 -5   190.6 -5   190.6 5
	patch quad  $matID  5 1   190.6 130   190.6 -130   207 -130   207 130
}


section fiberSecThermal  $U140  {
	patch	rect	$matID	5	1	-65	60	-5	70
	patch	rect	$matID	1	5	-65	-60	-58	60
	patch	rect	$matID	5	1	-65	-70	-5	-60
									
	patch	rect	$matID	5	1	5	60	65	70
	patch	rect	$matID	1	5	58	-60	65	60
	patch	rect	$matID	5	1	5	-70	65	-60
}

section fiberSecThermal  $U100  {
	patch	rect	$matID	5	1	-51	42.4	-5	50
	patch	rect	$matID	1	5	-51	-42.4	-46.5	42.4
	patch	rect	$matID	5	1	-51	-50	-5	-42.4
									
	patch	rect	$matID	5	1	5	42.4	51	50
	patch	rect	$matID	1	5	46.5	-42.4	51	42.4
	patch	rect	$matID	5	1	5	-50	51	-42.4

}


# element
set colTransf  1
set baemTransf 2
set brTransf   3
set lcTransf   4

geomTransf Linear $colTransf
geomTransf Linear $baemTransf
geomTransf Linear $brTransf
geomTransf PDelta $lcTransf

print "done2nd"

# column def
set section $W16x67

#element	dispBeamColumnThermal	$eleTag	$iNode	$jNode	$numIntgrPts	$secTag	$transfTag	<-mass	$massDens>

element	dispBeamColumnThermal	10101	1001	2001	5	$section	$colTransf
element	dispBeamColumnThermal	10102	1002	2002	5	$section	$colTransf
element	dispBeamColumnThermal	10103	1003	2003	5	$section	$colTransf
element	dispBeamColumnThermal	10104	1004	2004	5	$section	$colTransf
element	dispBeamColumnThermal	10105	1005	2005	5	$section	$colTransf
element	dispBeamColumnThermal	10106	1006	2006	5	$section	$colTransf
							
element	dispBeamColumnThermal	10201	2001	3001	5	$section	$colTransf
element	dispBeamColumnThermal	10202	2002	3002	5	$section	$colTransf
element	dispBeamColumnThermal	10203	2003	3003	5	$section	$colTransf
element	dispBeamColumnThermal	10204	2004	3004	5	$section	$colTransf
element	dispBeamColumnThermal	10205	2005	3005	5	$section	$colTransf
element	dispBeamColumnThermal	10206	2006	3006	5	$section	$colTransf
							
element	dispBeamColumnThermal	10301	3001	4001	5	$section	$colTransf
element	dispBeamColumnThermal	10302	3002	4002	5	$section	$colTransf
element	dispBeamColumnThermal	10303	3003	4003	5	$section	$colTransf
element	dispBeamColumnThermal	10304	3004	4004	5	$section	$colTransf
element	dispBeamColumnThermal	10305	3005	4005	5	$section	$colTransf
element	dispBeamColumnThermal	10306	3006	4006	5	$section	$colTransf
							
element	dispBeamColumnThermal	10401	4001	5001	5	$section	$colTransf
element	dispBeamColumnThermal	10402	4002	5002	5	$section	$colTransf
element	dispBeamColumnThermal	10403	4003	5003	5	$section	$colTransf
element	dispBeamColumnThermal	10404	4004	5004	5	$section	$colTransf
element	dispBeamColumnThermal	10405	4005	5005	5	$section	$colTransf
element	dispBeamColumnThermal	10406	4006	5006	5	$section	$colTransf
							
element	dispBeamColumnThermal	10501	5001	6001	5	$section	$colTransf
element	dispBeamColumnThermal	10502	5002	6002	5	$section	$colTransf
element	dispBeamColumnThermal	10503	5003	6003	5	$section	$colTransf
element	dispBeamColumnThermal	10504	5004	6004	5	$section	$colTransf
element	dispBeamColumnThermal	10505	5005	6005	5	$section	$colTransf
element	dispBeamColumnThermal	10506	5006	6006	5	$section	$colTransf




# beam def
set section1 $W16x26
set section2 $W16x36

element	dispBeamColumnThermal	20101	20011	20111	5	$section1	$baemTransf
element	dispBeamColumnThermal	20102	20021	20121	5	$section2	$baemTransf
element	dispBeamColumnThermal	20103	20031	20131	5	$section1	$baemTransf
element	dispBeamColumnThermal	20104	20041	20141	5	$section2	$baemTransf
element	dispBeamColumnThermal	20105	20051	20151	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20201	30011	30111	5	$section1	$baemTransf
element	dispBeamColumnThermal	20202	30021	30121	5	$section2	$baemTransf
element	dispBeamColumnThermal	20203	30031	30131	5	$section1	$baemTransf
element	dispBeamColumnThermal	20204	30041	30141	5	$section2	$baemTransf
element	dispBeamColumnThermal	20205	30051	30151	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20301	40011	40111	5	$section1	$baemTransf
element	dispBeamColumnThermal	20302	40021	40121	5	$section2	$baemTransf
element	dispBeamColumnThermal	20303	40031	40131	5	$section1	$baemTransf
element	dispBeamColumnThermal	20304	40041	40141	5	$section2	$baemTransf
element	dispBeamColumnThermal	20305	40051	40151	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20401	50011	50111	5	$section1	$baemTransf
element	dispBeamColumnThermal	20402	50021	50121	5	$section2	$baemTransf
element	dispBeamColumnThermal	20403	50031	50131	5	$section1	$baemTransf
element	dispBeamColumnThermal	20404	50041	50141	5	$section2	$baemTransf
element	dispBeamColumnThermal	20405	50051	50151	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20501	60011	60111	5	$section1	$baemTransf
element	dispBeamColumnThermal	20502	60021	60121	5	$section2	$baemTransf
element	dispBeamColumnThermal	20503	60031	60131	5	$section1	$baemTransf
element	dispBeamColumnThermal	20504	60041	60141	5	$section2	$baemTransf
element	dispBeamColumnThermal	20505	60051	60151	5	$section1	$baemTransf


# right
element	dispBeamColumnThermal	20111	20111	20022	5	$section1	$baemTransf
element	dispBeamColumnThermal	20112	20121	20032	5	$section2	$baemTransf
element	dispBeamColumnThermal	20113	20131	20042	5	$section1	$baemTransf
element	dispBeamColumnThermal	20114	20141	20052	5	$section2	$baemTransf
element	dispBeamColumnThermal	20115	20151	20062	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20211	30111	30022	5	$section1	$baemTransf
element	dispBeamColumnThermal	20212	30121	30032	5	$section2	$baemTransf
element	dispBeamColumnThermal	20213	30131	30042	5	$section1	$baemTransf
element	dispBeamColumnThermal	20214	30141	30052	5	$section2	$baemTransf
element	dispBeamColumnThermal	20215	30151	30062	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20311	40111	40022	5	$section1	$baemTransf
element	dispBeamColumnThermal	20312	40121	40032	5	$section2	$baemTransf
element	dispBeamColumnThermal	20313	40131	40042	5	$section1	$baemTransf
element	dispBeamColumnThermal	20314	40141	40052	5	$section2	$baemTransf
element	dispBeamColumnThermal	20315	40151	40062	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20411	50111	50022	5	$section1	$baemTransf
element	dispBeamColumnThermal	20412	50121	50032	5	$section2	$baemTransf
element	dispBeamColumnThermal	20413	50131	50042	5	$section1	$baemTransf
element	dispBeamColumnThermal	20414	50141	50052	5	$section2	$baemTransf
element	dispBeamColumnThermal	20415	50151	50062	5	$section1	$baemTransf
							
element	dispBeamColumnThermal	20511	60111	60022	5	$section1	$baemTransf
element	dispBeamColumnThermal	20512	60121	60032	5	$section2	$baemTransf
element	dispBeamColumnThermal	20513	60131	60042	5	$section1	$baemTransf
element	dispBeamColumnThermal	20514	60141	60052	5	$section2	$baemTransf
element	dispBeamColumnThermal	20515	60151	60062	5	$section1	$baemTransf



# brace
set sectionBr $U140
element	dispBeamColumnThermal	50101	10023	20034	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50102	10034	20023	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50103	10043	20054	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50104	10054	20043	5	$sectionBr	$brTransf
							
element	dispBeamColumnThermal	50201	20023	30034	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50202	20034	30023	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50203	20043	30054	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50204	20054	30043	5	$sectionBr	$brTransf
							
element	dispBeamColumnThermal	50301	30023	40034	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50302	30034	40023	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50303	30043	40054	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50304	30054	40043	5	$sectionBr	$brTransf
							
element	dispBeamColumnThermal	50401	40023	50034	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50402	40034	50023	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50403	40043	50054	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50404	40054	50043	5	$sectionBr	$brTransf
							
set sectionBr $U100							
element	dispBeamColumnThermal	50501	50023	60034	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50502	50034	60023	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50503	50043	60054	5	$sectionBr	$brTransf
element	dispBeamColumnThermal	50504	50054	60043	5	$sectionBr	$brTransf

# loading

# Gravity Loads
timeSeries Linear 100 -factor  1.0
pattern Plain 110 100 { 
	eleLoad	-ele	20101	-type	-beamUniform	-5.884  0.0 
	eleLoad	-ele	20102	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20103	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20104	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20105	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20201	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20202	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20203	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20204	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20205	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20301	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20302	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20303	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20304	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20305	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20401	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20402	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20403	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20404	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20405	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20501	-type	-beamUniform	-2.256
	eleLoad	-ele	20502	-type	-beamUniform	-2.256
	eleLoad	-ele	20503	-type	-beamUniform	-2.256
	eleLoad	-ele	20504	-type	-beamUniform	-2.256
	eleLoad	-ele	20505	-type	-beamUniform	-2.256
							
	# right
	eleLoad	-ele	20111	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20112	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20113	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20114	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20115	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20211	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20212	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20213	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20214	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20215	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20311	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20312	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20313	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20314	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20315	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20411	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20412	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20413	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20414	-type	-beamUniform	-5.884  0.0
	eleLoad	-ele	20415	-type	-beamUniform	-5.884  0.0
						
	eleLoad	-ele	20511	-type	-beamUniform	-2.256
	eleLoad	-ele	20512	-type	-beamUniform	-2.256
	eleLoad	-ele	20513	-type	-beamUniform	-2.256
	eleLoad	-ele	20514	-type	-beamUniform	-2.256
	eleLoad	-ele	20515	-type	-beamUniform	-2.256
						

}

# leaning col.
# ###############################################################
# node
node	777001	36000	0
node	777002	36000	3000
node	777003	36000	6000
node	777004	36000	9000
node	777005	36000	12000
node	777006	36000	15000


# fix
fix 777001 1 1 0

# rigid links
rigidLink	bar	2006	777002
rigidLink	bar	3006	777003
rigidLink	bar	4006	777004
rigidLink	bar	5006	777005
rigidLink	bar	6006	777006

element	elasticBeamColumn	1110101	777001	777002	537094	2.00E+05	15960000000	$lcTransf
element	elasticBeamColumn	1110102	777002	777003	537094	2.00E+05	15960000000	$lcTransf
element	elasticBeamColumn	1110103	777003	777004	537094	2.00E+05	15960000000	$lcTransf
element	elasticBeamColumn	1110104	777004	777005	537094	2.00E+05	15960000000	$lcTransf
element	elasticBeamColumn	1110105	777005	777006	537094	2.00E+05	15960000000	$lcTransf


pattern Plain 111 100 { 
	load	777002	0.0	-176520.0	0.0
	load	777003	0.0	-176520.0	0.0
	load	777004	0.0	-176520.0	0.0
	load	777005	0.0	-176520.0	0.0
	load	777006	0.0	-67680.0	0.0
}
# ###############################################################


# Rayleigh damping
# set damping based on first eigen mode
set lambda [expr [eigen 1]**0.5]
puts "First Mode Period: T= [expr 6.283185/$lambda]\n"

set dampRatio 0.05
rayleigh 0. 0. 0. [expr 2*$dampRatio/$lambda]

# Analysis options
system UmfPack
numberer RCM
constraints Transformation
integrator LoadControl 0.1
test NormDispIncr 1.0e-12  10 3
algorithm Linear
analysis Static
analyze 10
puts "Gravity Analysis completed SUCCESSFULLY\n"
loadConst -time 0.0


# recorder
# Node recorder
recorder Node -file Node_displacementsY_1.txt   -time -node 	20111	20121	20131	20141	20151	-dof 2  disp
recorder Node -file Node_displacementsY_2.txt   -time -node 	30111	30121	30131	30141	30151	-dof 2  disp
recorder Node -file Node_displacementsY_3.txt   -time -node 	40111	40121	40131	40141	40151	-dof 2  disp
recorder Node -file Node_displacementsY_4.txt   -time -node 	50111	50121	50131	50141	50151	-dof 2  disp
recorder Node -file Node_displacementsY_5.txt   -time -node 	60111	60121	60131	60141	60151	-dof 2  disp

# Drift
recorder Drift -file driftX.txt -time -iNode 1001	2001	3001	4001	5001 -jNode 2001	3001	4001	5001	6001 -dof 1 -perpDirn 2

# Element Recorder
recorder Element -file DispBeamColumn_localForce.txt -time -ele 20101 localForce
recorder Element -file DispBeamColumn_TempElong.txt -time -ele 20101 section 1 fiber 0 0 TempElong
# display displacement shape of the column
recorder display "Displaced shape" 650 20 500 550 -wipe
prp  40 70 1;
vup  0  1  0;
vpn  0  0  1;
display 1 5 200;


# Time History Analysis
set GMdirection 1
set DtAnalysis 0.01
set TmaxAnalysis 36.88
set inFile "NORTHR.AT2"
set dt 0.01

# Define Pattern for Time History analysis
set AccelSeries 200
timeSeries Path $AccelSeries -dt $dt -filePath $inFile -factor 9810
pattern UniformExcitation 250  $GMdirection -accel $AccelSeries

# Analysis options

integrator Newmark 0.5 0.25
test NormUnbalance 1.0e-6 150;
algorithm ModifiedNewton -initial
analysis Transient

set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];

# perform analysis - returns 0 if analysis was successful
# analyze $Nsteps $DtAnalysis;

puts "Time History Analysis completed SUCCESSFULLY"
remove loadPattern 250

#fire loading
puts "Fire Loading"
loadConst -time 0.0	

set fThermalLoad1  thermalTH1.txt
set fThermalLoad2  thermalTH2.txt
set fThermalLoad3  thermalTH3.txt

set HalfD1  [expr 399./2.]
set HalfD2  [expr 403./2.]
set HalfDC1 [expr 415./2.]
set HalfDBr1 70
set HalfDBr2 50



pattern Plain 500 Linear {
	# beam
	eleLoad	-ele	20101	-type	-beamThermal	-source	TH0.txt	-$HalfD1	$HalfD1
	
# 
}

set nAllFire  3600
set nFiStep [expr 1.0/$nAllFire]

# Analysis options
test NormUnbalance 1.0e-3 150 2
system UmfPack
integrator LoadControl 0.01
analysis Static
analyze 100