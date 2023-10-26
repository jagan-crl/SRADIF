"""
How to Cite:

J.K. Tavares, V. Gururajan, J. Jayachandran, Effects of radiation heat loss on planar and spherical
hydrofluorocarbon/air flames, Combust. Flame 258 (2023) 113067.

"""

import cantera as ct
import numpy as np
import SRADIF as OTL
import time

######################################################################
#  General Controls
######################################################################                                                                                               

# Polynomial fit for exp data
POLY = False

######################################################################
#  Shell Generation Controls
######################################################################

# Number of layers, default was 1000                                                               
NLAY = 1000 #So that polyfit uses odd number

# Chamber volume in cm^3:
Chamber_radius = 3.0e0
VOLU = 4.0e0/3.0e0*np.pi*Chamber_radius**3.0e0

# Ignition shell volume (LTORC data only, set to 0 to ignore)
ign = 0.0e0 #cm
IGN = (4.0e0/3.0e0*np.pi*ign**3.0e0)/7.0e0

######################################################################
#  Tolerance Controls
######################################################################

# Volume tolerance in cm^3                                                       
ATOL = 1.0E-6

# Gama tolerance                                                                 
GTOL = 1.0E-6

# Cooling time tolerance in s                                                    
TTOL = 1.0E-6

######################################################################
#  Gas Initial Conditions
######################################################################

# Initial temperature                                                            
TEMP = 300.0e0

# Initial pressure in atm                                                          
PRES = 1.0e0

# Reactant composition (moles or mole fraction)                              
reactants = 'CH2F2:1.2, O2:1.0, N2:3.76'

# Chemical mechanism
scheme = 'hfcmech_R32_2021Aug20.cti'

#######################################################################

startAll = time.time()# program timer
OTL.spherical(reactants,scheme,PRES,TEMP,NLAY,VOLU,\
                ATOL,GTOL,TTOL,POLY,IGN)
endAll = time.time() # end time
print('  Elapsed Time {0:4.2f} seconds'.format(endAll-startAll))
