For some experiments, we want to be able to provide the mesher with 
"poles" (points of the medial axis) that are given by some external 
software.

** POWER CRUST **
1. To use the power crust algorithm to compute these points, we can 
download code from :
http://code.google.com/p/powercrust/

2. Use the example off_to_pts.cpp to generate an input file for power crust
command line : 
>./off_to_pts.exe myshape.off
Note myshape.off should contain a point cloud describing the shape. No
connectivity is needed here.

3. Give it to powercrust
./powercrust.exe -i myshape.pts

4. axisface.off gives an approximation of the medial axis

3. /*todo : give it to the mesher*/



** ANOTHER ONE **
- where to find the code
- how to use it