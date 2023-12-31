## NOTES : Written by Gregory Mohammed, Masters Thesis,
###		  31 August, 2012

EDITS : Written by Gregory Mohammed, after successful completion of the Masters thesis. As of 3 September, 2013 all bugs have been resolved.  

There are a lot of switches and data in here. Read through carefully so that you understand the comments and what each switch and set of data do and are for. Do some test runs if you are not clear.

**Indicate choices:**
Initial conditions: 0 = smooth shock, 1 = shock   
initShock = 1

**showInitialSlice: Show = 1, Don't show = 0**

showExact: 0 => Does not show an exact solution.   
1 => Plots the exact solution.   
Be aware that there is no correct exact solution for the KKT data.

showInitialSlice = 0   
showExact = 0

## !!! NOTE !!!
showExact has to be 1 for the Upwind solution to display. This is also true for the Godunov solution. The KKT solution does not have this problem, but showExact must equal 0. For KKT as well, method = 1, neutrinos = 1, dataSwitch = 2, and vprofile = 2. You can play with any of the neutrino fluxes and the totalTimeKKT value (in seconds).

**Method: Upwind = 0, Godunov = 1**

## !!! NOTE !!! - For either of these methods ensure that showExact = 1 (above).

method = 1

**Neutrinos: Do Not Include = 0, Include = 1**

neutrinos = 1

**Tunes the value of the neutrino flux and energy term.**

fluxTuner = 1.0e-2 (fluxTuner = 0.0e-0)   
energyTuner = 1.0e-2 (energyTuner = 0.0e-0)

**High Resolution Method: ZeroSigma = 0, MinMod = 1, MC = 2**

highResType = 1

**Constant = 0, Linear = 1, CubicHermite = 2**

reconstructionType = 1

**Resolution, measured in the number of cells.**

resolution = 1600

**Test file with parameters for the program Godunov-v.2.0.0, which is written in Java**

gamma = 2.0 (gamma = 1.67)

**for gamma=2, k=epsilon/rho0**

k = 0.125

**number of evolutions**

totalSteps = 10000 (totalSteps = 1)

**Total time for the Sod tube.**

The Sod values are hard-coded for the Sod shock tube. 
This allows the param.dat file to be flexible for "playing" with data to observe results for different physical situations (or unphysical ones too).

**THE ONLY VARIABLE WHICH NEEDS TO BE EDITED HERE IS THE totalTimeSod (below).**
**Set to 0.25 when dataSwitch = 0. These times are already geometrized.**

totalTimeSod = 0.25 (totalTimeSod = 0.114)

Total time for the Kuroda et al. data. This gives a good evolution which shows good data. vprofile can be what you want to see. 

totalTimeKKT = 3.33e-5 (totalTimeKKT = 2.5e-5)

**dump control: allowed change in evolution variables**

epsdmp = 1.0

**interval between dumps**

dmpinterval = 1000

**number of correctors for Upwind.**

corrector = 1

**artificial viscosity**

artvis.k1 = 0.0
artvis.k2 = 0.0

**courant**
**initial delta = p1*freefall**

courant = 0.4 (p1 small for KKT, and Sod with total time > 0.25.p1 = 0.4)

**dataSwitch: Sod = 0, Experimental Sod = 1, Neutrino Model = 2**

If dataSwitch = 0, then the classic Sod tube will be plotted using the hard-coded data. If dataSwitch = 1, then your data will be plotted. If dataSwitch = 2, then the KKT data will be plotted.

dataSwitch = 2

If plotGeo = 1, then plot geometrized units, otherwise plot the ungeometrized units. This is only for KKT. For the Sod data and experimental Sod data the value of this variable is irrelevant.

plotGeo = 0

**Data for the Sod shock tube. Feel free to experiment here, as the classic Sod data are hard-coded, and unaffected by changes here.**

x_left = -0.3
x_right = 0.3

v_left = 0.0
v_right = 0.0

rho0_left = 1.0
rho0_right = 0.125

p_left = 1.0
p_right = 0.1

This data below is ONLY for the Kuroda et al. data against the data used in the exact Riemann solver.

**Neutrino data from Kuroda et al.**

**Here use cgs units.**

xLeft = 8.0e6
xRight = 1.0e8

vLeft = 1.0e7
vRight = -1.0e7

lrho0nu = 2.0e14
rrho0nu = 1.0e9

lpnu = 1.07e8
rpnu = 0.0

Make different velocity profiles. Other variable profiles are hard-coded as noted below.   
vprofile = 1 => Discontinuous, velocity is zero, and the rho0 and p are steps. This is compared with the original Sod data.   
vprofile = 2 => parabolic, left increasing, right decreasing. Other variables are negative ramps. This is compared with the Sod tube for similar input.

vprofile = 2

_THIS IS NOT WORKING. IT WAS INCLUDED FOR FUTURE WORK._
_Mass of the neutronized core, or black hole._
_This is specified by the user, and may be 0 if the user does not want to consider it. However, doing so zeros out the ad-hoc gravity term, which leads to a numerical explosion, which is non-physical._

**Mass of the Sun = 1.98892e30 kg**
**Geometrized Mass of the Sun = 4.5 km**

massCore = 0.0 (massCore = 4.5e-14)
