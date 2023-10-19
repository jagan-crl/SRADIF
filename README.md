## How to Cite

J.K. Tavares, V. Gururajan, J. Jayachandran, Effects of radiation heat loss on planar and spherical
hydrofluorocarbon/air flames, Combust. Flame 258 (2023) 113067.

## Required Packages

Python
Cantera
Numpy

## Input File

The input file (Input.py) should be modified to specify:
-Initial temperature
-Initial pressure
-Unburned reactant composition
-Chemical model

Other input file control parameters:
-Number of layers
-Initial unburned gas ball radius
-Ignition burned gas ball radius
-Tolerances (ATOL,GTOL,TTOL)
-Switch for polynomial fit of experimental data

## Main SRADIF Code

The main SRADIF code (SRADIF.py) should also be edited to include the experimental data input file name.
The column numbers for time and flame radius should also be specified considering units (must convert to ms and cm):

Line 169:
expdata = np.loadtxt('<experimental data file name>', encoding="utf8")

Line 177-177:
texp = expdata[:,<time column number>] * 1000					# [s] to [ms]
rexp = expdata[:,<flame radius column number>] * 100 			# [m] to [cm]

## How to Run

Run the input file (Input.py), which will call main function "spherical" from main SRADIF code (SRADIF.py).
All other functions are specified within the main SRADIF code.

# Outputs

1. Output.txt

Column 1: Time 								[ms]
Column 2: Flame Radius 						[cm]
Column 3: Flame Stretch						[s^-1]
Column 4: Flame Propagation Speed 			[cm/s]
Column 5: Inward Flow Velocity				[cm/s]
Column 6: Burned Flame Speed				[cm/s]

2. profiles.dat

Column 1: Radial Distance					[cm]
Column 2: Temperature						[K]
Column 3: Radiation-induced Gas Velocity 	[cm/s]
Column 4: H2O Mole Fraction
Column 5: CO2 Mole Fraction
Column 6: CO Mole Fraction
Column 7: CH4 Mole Fraction
Column 8: HF Mole Fraction
Column 9: CF2O Mole Fraction

Data in profiles.dat is printed in time intervals specified with the "writeDeltat" parameter (Line 225) of the main SRADIF code. 

## Core Algorithm

The core algorithm is detailed in the paper "Effects of radiation heat loss on planar and spherical hydrofluorocarbon/air flames" by Justin Tavares and Jagannath Jayachandran.


## Example

An example input file (Input_Example.py) with corresponding SRADIF code file (SRADIF_Example.py), which loads data from "globaloutput_Example.dat,"
is provided in the "Example" folder along with the resulting output files "Output.txt" and "profiles.dat." 

Note: This requires that the included 'hfcmech_R32_2021Aug20.cti' chemical mechanism be added to Cantera's chemical model library folder.


