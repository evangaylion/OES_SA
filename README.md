# OES_SA
Optical Emission Spectroscopy with Self Absorption

This code was developed for the dissertation "Modeling Extreme Nuclear Fireball Vapor Chemistry via Laser Ablation Plasma Plumes"

Kurucz.py and NIST_ASD.py return an array of transition data for a given ion and wavelength range. These parameters along with the temperature and other physical parameters are used to compute the emission spectrum. 

Self absorption parameters have also been included. There is the option to use Gaussian, Lorentzian, or Voigt line profiles, though the default setting is Voigt. Calculations are preformed assuming local thermodynamic equilibrium.

The plume_parsing.py code was developed with 2D CFD simulation images in mind. The CFD program in question was HyBurn developed by R. Houim. In order to take advantage of this script, one would need to have 2D data of the distribution of a particular atomic species, temperature, and density. Using VisIt or other Amrex-compatible visualization programs, a lineout can be "drawn" on the 2D simulation space, collection 1D data for the listed variables. One can also simply establish a set of 1D arrays meant to mimic this.

The parsing code mainly serves to create a line-of-sight from which a spatially-integrated signal is produced.

