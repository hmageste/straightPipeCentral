# straightPipe
Straight pipe OpenFoam case for low Reynolds turbulent model

Here it is worth noting that the user should issue the following commands to run the case:
- blockMesh : responsible to generate the mesh;
- simpleFoam : resolve the case;
- simpleFoam -postProcessing -func wallShearSteess : post process the case calculating the shear stress;
- paraFoam : visualize the results;

Initial values for $$k$$ and $$\epsilon$$ usually are taken as:
- $$k=\frac{3}{2}(lU)^2$$;
- $$C_\mu^\frac{3}{4}\frac{k^\frac{3}{2}}{l}$$;

where $$l=0.16Re^{-\frac{1}{8}}$$.


