wipe;
#units are m, N, MPa, ton, and C

# Model builder
model BasicBuilder -ndm 3 -ndf 6;
file mkdir output
source DisplayPlane.tcl
source DisplayModel2D.tcl
source DisplayModel3D.tcl


# Execute
source sections.tcl
puts "loaded 'sections' correctly"
source nodes.tcl
puts "loaded 'nodes' correctly"
source constraints.tcl
puts "loaded 'constraints' correctly"
source elements.tcl
puts "loaded 'elements' correctly"
source loads.tcl

