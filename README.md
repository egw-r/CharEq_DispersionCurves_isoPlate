# CharEq_DispersionCurves_isoPlate

# 1.0 isotropic_plate_intersectionSolution.py
This program was written for calculating the dispersion curves of an isotropic elastodynamic plate waveguide. It solves 
Horace Lamb's characteristic equation for the frequency-wavenumber pairs that make up the dispersion curves. It takes 
the thickness of the waveguide, its material properties, frequency range of insterest and phase velocity of interest as 
input. It outputs a plot of the dispersion curves of a Lamb wave, separated symmetric and anti-symmetric modes, and a 
plot of 'error' which is the difference in wavenumber between adjacent points when the refinement process of the 
root-finding is done. It also saves the plotted data into comma-separated text files separated into symmetric and anti-symmetric
solutions. The data is formatted into three columns, frequency, wavenumber, and error. The error can be thought of as an
error in the wavenumber root at a given frequency.

For better or worse the output is frequency with units of Hz and the wavenumber is units of rad/m. I would have rather 
it be rad/s instead of Hz, but oh well. The plotting is there, in part, to let the user know how to manipulate the data 
that is output by the solver.

## 1.1 Horace Lamb
The original dispersion relation, or characteristic equation, was derived by Horace Lamb:

Horace Lamb. “On Waves in an Elastic Plate.” Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character 93, no. 648 (March 1, 1917): 114–28. https://doi.org/10.1098/rspa.1917.0008.

The original included hyperbolic functions, but modern representations use trigonometric functions:

Anti-symmetric Lamb Waves:

![images](images/antisym_chareq.svg)

Symmetric Lamb Waves:

![images](images/sym_chareq.svg)

where ![images](images/longWavenumber.svg) and ![images](images/shearWavenumber.svg).

## 1.2 Root-Finding Method

The root-finding method compares three adjacent points. The center point is a root if the trend of the points go from 
positive to negative (transition over zero) and if the slope of the absolute value of the points goes from negative to 
positive (a valley). For details, please see the commented code.

### 1.2.1 Prerequisites

As shown in isotropic_plate_intersectionSolution.py, the required packages are very little: 
```
import numpy as np
import elasticity as es
import time
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools impart partial
import matplotlib as mpl
```
'elasticity' is elasticity.py and the rest are common scientific packages: numpy and matplotlib.

'time' and 'multiprocessing' are to keep track of how long a full calculation takes and allows all the 
threads of your computer processor to be used for calculation purposes. Partial is used for formatting the solver function
so that it can be used by multiprocessing.

## 1.3 Installing

This code was written with readability and ease-of-use in mind. As a researcher coming to python from Matlab, this code 
is similar in structure and does not make any attempt to optimize for python. It does not even try that hard to be efficient.
As such, the user need only: 
1. download elasticity.py and isotropic_plate_intersectionSolution.py, 
1. place the two files in the same directory/folder,
1. edit the labelled global parameters in isotropic_plate_intersectionSolution.py (and maybe some of the num_points variables if it runs too slow),
1. run isotropic_plate_intersectionSolution.py.

For reference, I used pycharm to run and develop the code. The interpreter was python3.

## 1.4 Performance and Tests

The processor I used was a Ryzen 3700x, it has 8 cores and 16 threads. The base clock is 3.6 GHz.

When tested by calculating solutions between phase velocities of 200 m/s and 10000 m/s with 200 frequency points 
between 1 Hz and 20 MHz, the symmetric and anti-symmetric modes were calculated in about 96 sec with a conservative 
estimate of the accuracy being 10^-3 rad/m.

The SH0 mode appears in the dispersion curve solutions because both the LHS and RHS of the characteristic equation 
become zero when the phase velocity is equal to the shear wave speed.

# 2.0 groupvelocity_calculator.py
This program was written for calculating the group velocity of Lamb waves when given only the frequency-wavenumber 
solutions to the characteristic equation of a Lamb wave. Some group velocity calculation methods require the dispersion
curves to be sorted, which can result in slower calculations and introduce errors if sorted incorrectly. These methods also
rely on numerical derivatives, the accuracy of which depends on the density of points along the curve. Having a higher density
of points requires significantly more computational time. To avoid these setbacks, the code uses a method published in the 
following article:

Hakoda, Christopher and Lissenden, Cliff J. "Application of a general expression for the group velocity vector of elastodynamic guided waves," Journal of Sound and Vibration, 469 (2020): 115165.https://doi.org/10.1016/j.jsv.2019.115165.

The method uses the displacement profile of the Lamb wave, as opposed to adjacent points in a dispersion curve, to calculate
the group velocity. For anisotropic plate waveguides these displacement profiles can be calculated using partial waves, but 
for Lamb wave's the process is trivial due to the closed-form Lamb wave solution. 

Included in the code is a method for calculating the group velocity of multiple frequency-wavenumber solutions and another
method for calculating the group velocity of a single frequency-wavenumber solution. After plotting the phase-velocity
-versus-frequency and the group-velocity-versus-frequency plot, the code also saves the group velocity
data in a new comma-separated text file that is separated into symmetric and anti-symmetric solutions. The data is formatted
into three columns, frequency, wavenumber, and group velocity. 

PLEASE NOTE: the SH0 mode's displacement profile is not accurately represented by the displacement profile expressions derived for Lamb waves
so the group velocity values calculated for the SH0 modes are not correct. 

## 2.1 Prerequisites
Just like in isotropic_plate_intersectionSolution.py, the required packages are very little and their purposes are much the same: 
```
import numpy as np
import elasticity as es
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
```

# 3.0 dispersion_plotter.py
The code was written mainly for plotting the dispersion curves and group velocity curves, but was also included as a reference
for anyone that wants to plot the results with different formatting. It includes two simple methods that are an example of how
the data can be read from the saved comma-separated text files and then plotted neatly. 

## 3.1 Prerequisites
The prerequisite packages for this code are even less: 
```
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
```

# Authors

* **Christopher Hakoda** - *Initial work* - [egw-r](https://github.com/egw-r/)

# License

This project is licensed under the GNU v3.0 License - see the [LICENSE](LICENSE) file for details
