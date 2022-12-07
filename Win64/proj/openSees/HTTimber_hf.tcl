
#This Tcl file is written for Heat Transfer in a timber section modelled using OpenSees


wipe;
HeatTransfer 2D;         #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 


#HTMaterial definition
HTMaterial Timber 2 1 -file TImberPars.dat 5.0 60.0 30.0 600e3;     #Defining HeatTransfer Material with Material tag:2, type tag:1 (any integer);

#Dimension definition
set blockw 0.1;  #Timber block width
set blockt 0.06;	# Timber block thickness (depth)

#HT Entity
HTEntity Block 2 0.0 [expr $blockt/2] $blockw $blockt;


#Mesh control(element size)
set elex [expr $blockw/10];
set eley [expr $blockt/60];

HTMesh 1 2 2 -phaseChange 1 -MeshCtrls $elex $eley;  

#To activate mesh
HTMeshAll;
puts "Mesh complete"

#Set ambient Temperature
SetInitialT 298.15;

#Define HTconstants tag, convection coefficient,ambient temp.,emissity,,absorptivity
HTConstants 1 0.0 298.15 0.8 0.8; #Unexposed surface (ambient air)
HTConstants 2 25.0 298.15 0.8 0.8; #Exposed surface (fire)

#Define Heat flux Boundary condition (Heat loss to air)
HTPattern AmbientBC 1 {
	HeatFluxBC -HTEntity 2 -face 4 -type ConvecAndRad -HTConstants 1;  #for line
	#HeatFluxBC -HTEntity 2 -face 4 -type ConvecAndRad -HTConstants 1; 
}

#Define fire models
#Standard Fire Curve
FireModel standard 1; #-duration 3600; 
#FireModel hydroCarbon 1;
#User defined fire curve: fire dat(time, gasTemp)
#FireModel UserDefined 3 -file fire.dat; #-duration 3600; 
#EC1 Localised fire
FireModel Localised 2 -origin 0.0 -3.0 0.0 -firePars 1.0 2E6 3.0 2 ;
#Idealised heat flux : uniform heat flux q
FireModel Idealised 5 -q 50e3 -uniform;    #for fixed heat flux
FireModel UserDefined 6 -file hf.dat -type 2; #for heat flux definition using external file (time-heat flux data)


#Define Heat Flux Boundary Condition (Heat received from fire)
HTPattern fire 2 model 6 {
	HeatFluxBC -HTEntity 2 -face 1 -type -Prescribed -HTConstants 2; #Fixed heat flux as cone
	#HeatFluxBC -HTEntity 1 -face 1 -type -ConvecAndRad -HTConstants 2;
	HeatFluxBC -HTEntity 2 -face 1 -type -ConvecAndRad -HTConstants 2;
}


#Record concrete beam temperatures
HTRecorder -file TimberTdist.out -xloc 0.0 -yloc 0.001 0.010 0.018 0.042 0.054 0.060;


#Define node set for recording temperature change;
HTNodeSet 1 -Locx 0.0  ;
HTRecorder -file Timbertemp.out NodeSet 1;

HTEleSet 1 -HTEntity 2 -NodeSet 1 -face 2;
HTRecorder -file TimberChar.out EleSet 1 material 1 -phaseTag;

#HTNodeSet 1 -Entity 1 -Face 1 -Locx 0.0 -Locy -0.05 -Locz 0.0;

#Define HT analysis type (default)
HTAnalysis HeatTransfer

#Define number of analysis steps and step size
#HTAnalyze $numSteps $timeStep;
HTAnalyze 8800 0.5;	

#Wipe the analysis when completed
wipe;

	
	