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

## Monitoring the data

### For monitoring the residuals do:
foamMonitor -l postProcessing/residuals/0/residuals.dat

### Ploting the graph is done with:
gnuplot
gnuplot> set style data linespoints
gnuplot> plot "postProcessing/singleGraph/<time>/line_U.xy" u 2:1

## Skin factor
For a circular pipe, the following relations are true:
- f = 8*tau/(pho*U^2)
- delta_p = 1/4*f*L/R*rho*u^2


