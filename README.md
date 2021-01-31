# CharEq_DispersionCurves_isoPlate

This program was made for calculating the dispersion curves of an isotropic elastodynamic plate waveguide. It solves 
Horace Lamb's characteristic equation for the frequency-wavenumber pairs that make up the dispersion curves. It takes 
the thickness of the waveguide, its material properties, frequency range of insterest and phase velocity of interest as 
input. It outputs a plot of the dispersion curves of a Lamb wave, separated symmetric and anti-symmetric modes, and a 
plot of 'error' which is the difference in wavenumber between adjacent points when the refinement process of the 
root-finding is done. It also saves the plotted data into a comma-separated text file in the adjacent directory. 

For better or worse the output is frequency with units of Hz and the wavenumber is units of rad/m. I would have rather 
it be rad/s instead of Hz, but oh well. The plotting is there, in part, to let the user know how to manipulate the data 
that is output by the solver.

## Horace Lamb
The original dispersion relation, or characteristic equation, was derived by Horace Lamb:

Horace Lamb. “On Waves in an Elastic Plate.” Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character 93, no. 648 (March 1, 1917): 114–28. https://doi.org/10.1098/rspa.1917.0008.

The original included hyperbolic functions, but modern representations use trigonometric functions:

Anti-symmetric Lamb Waves:

![images](images/antisym_chareq.svg)

Symmetric Lamb Waves:

![images](images/sym_chareq.svg)

where ![images](images/longWavenumber.svg) and ![images](images/shearWavenumber.svg).

## Root-Finding Method

The root-finding method compares three adjacent points. The center point is a root if the trend of the points go from 
positive to negative (transition over zero) and if the slope of the absolute value of the points goes from negative to 
positive (a valley). For details, please see the commented code.

### Prerequisites

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

### Installing

This code was written with readability and ease-of-use in mind. As a researcher coming to python from Matlab, this code 
is similar in structure and does not make any attempt to optimize for python. It does not even try that hard to be efficient.
As such, the user need only: 
1. download elasticity.py and isotropic_plate_intersectionSolution.py, 
1. place the two files in the same directory/folder,
1. edit the labelled global parameters in isotropic_plate_intersectionSolution.py (and maybe some of the num_points variables if it runs too slow),
1. run isotropic_plate_intersectionSolution.py.

For reference, I used pycharm to run and develop the code. The interpreter was python3.

## Performance and Tests

The processor I used was a Ryzen 3700x, it has 8 cores and 16 threads. The base clock is 3.6 GHz.

When tested by calculating solutions between phase velocities of 200 m/s and 10000 m/s with 200 frequency points 
between 1 Hz and 20 MHz, the symmetric and anti-symmetric modes were calculated in about 96 sec with a conservative 
estimate of the accuracy being 10^-3 rad/m.

The SH0 mode appears in the dispersion curve solutions because both the LHS and RHS of the characteristic equation 
become zero when the phase velocity is equal to the shear wave speed.

## Authors

* **Christopher Hakoda** - *Initial work* - [egw-r](https://github.com/egw-r/)

## License

This project is licensed under the GNU v3.0 License - see the [LICENSE](LICENSE) file for details
