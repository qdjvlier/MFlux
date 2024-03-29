# MFlux
This software calculates root water uptake and soil and crop potentials according to a process-based model, De Jong van Lier et al. (2013), http://dx.doi.org/10.2136/vzj2013.02.0039.
Scenario parameters are defined in the file MFluxInfo.dat.
Optionally, time-variable potential transpiration rates can be defined in the file MFluxTpData.dat.
The simulation scenario includes the number of days to simulate, the internal timestep, the number of soil layers, the potential transpiration rate (Tp), root and xylem radii, longitudinal and radial root conductance, and minimum leaf water potential.  
A soil profile with N layers of variable thickness can be defined. For each layer, soil hydraulic properties are described by the Van Genuchten-Mualem parameters, as well as the root length density.
The model output includes, for predefined times during the simulated period, the maximum (Tm) and actual (Ta) transpiration rates (as defined in de Jong van Lier et al., 2013), the root water uptake from the soil layers, the water contents and pressure heads of the soil layers, the pressure head at the surface of the roots in the soil layers, and the xylem and leaf water potentials.
