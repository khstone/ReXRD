# ReXRD
Python script for calculating the resonant diffraction response of a powder for a given structure. 
The script can parse a simple CIF file to extract a structure and then calculate the response of 
specified Bragg peaks across a given x-ray energy range.  Averages over equivalent or overlapped 
peaks as appropriate for powder diffraction are made.

To run the simluation program, run ReXRD_sim.py from the command line, you will be prompted for
a CIF file to load. The program will parse the CIF file to extract the structural information, 
partial occupancies of sites is considered so that occupational disorder may be considered. 

You will then be prompted for the low and high x-ray energy ranges to simulate and the x-ray 
energy step size (1-5eV is generally sufficient). You will then be prompted for the Miller 
indices of the peaks that you wish to consider.  This program uses the energy dependent atomic 
scattering factors as tabulated by Henke et al. and available through the CXRO website 
(www.cxro.lbl.gov).
