# cpalab
Repository for my o-mesh generation code.
Several pieces of functionality need to be added, the foremost of which at present is a good form of mesh smoothing. At present the base mesh has terrible orthogonality, which I cant understand, however that is what it is. This needs to be smoothed into a much more useable mesh, which I hope to create a flow solver that will use the mesh created. 

Current functionality is being capable to mesh around a cirle, ellipse, or NACA 4 or 5 digit aerofoil. Planned functionality for taking an input file of coordinates, which allows any aerofoil shape to be produced. Inputs are aerofoil number, the mesh distance, number of points on shape, number of points in trailing edge truncation, number of points to the boundary, and starting cell height or aspect ratio. 

Output is to an ASCII .dat file, ready to be used by Tecplot360.
