HeatTransfer 2D;
HTMaterial CarbonSteelEC3 1;
HTMaterial ConcreteEC2 2 0.03;
HTEntity Block 1 0.0 0.000375 0.01 0.00075; #Steel profile
 HTEntity Block 2 0.0 0.06575 0.01 0.13;
#HTEntity Block 2 0.0 0.03825 0.1 0.075;
HTMesh 1 1 1 -phaseChange 0 -MeshCtrls 0.0002 0.0001875;
#HTMesh 2 2 2 -phaseChange 1 -MeshCtrls 0.1575 0.01625;
HTMesh 2 2 2 -phaseChange 1 -MeshCtrls 0.0002 0.004063;
#HTRefineMesh -Entity 2 -SeedTag 1 4 -space 0.02 10 0.01 9 0.005 4 0.01 9 0.02 10;
HTMeshAll;

SetInitialT 293.15;
HTNodeSet 1 -Entity 1 -face 4;
HTNodeSet 2 -Entity 2 -face 1;
HTCoupleT -NodeSet 1 2;
HTConstants 1 4.0 293.15 0.7 0.7;
HTConstants 2 25.0 293.15 0.7 0.7;

HTPattern AmbientBC 1 {
HeatFluxBC -HTEntity 2 -face 4 -type ConvecAndRad -HTConstants 1;
}
FireModel standard 1;
HTNodeSet 3 -Entity 2 -Locx 0.0;
# HTEleSet 1 -Entity 2 -NodeSet 3 -face 1;
# HTNodeSet 4 -Entity 2 -Locx 0.1 0.3;
# HTEleSet 2 -Entity 2 -NodeSet 4 -face 1;
HTPattern fire 2 model 1 {
HeatFluxBC -HTEntity 1 -face 1  -type ConvecAndRad -HTConstants 2;
# HeatFluxBC -HTEleSet 1 -face 1 -type ConvecAndRad -HTConstants 2;
# HeatFluxBC -HTEleSet 2 -face 1 -type ConvecAndRad -HTConstants 2;
}
HTRecorder -file temp0.out -NodeSet 3;
#HTRecorder -file TcomTdist.out -xloc 0.0 -yloc 0.018 0.042 0.054 0.060;

#HTRecorder -file temp1.out -NodeSet 2;
HTAnalysis HeatTransfer
HTAnalyze 360 30;

wipeHT;