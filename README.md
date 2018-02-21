# straightPipe
Straight pipe OpenFoam case for low Reynolds turbulent model

Here it is worth noting that the user should issue the following commands to run the case:
- blockMesh : responsible to generate the mesh;
- simpleFoam : resolve the case;
- simpleFoam -postProcessing -func wallShearSteess : post process the case calculating the shear stress;
- paraFoam : visualize the results;
