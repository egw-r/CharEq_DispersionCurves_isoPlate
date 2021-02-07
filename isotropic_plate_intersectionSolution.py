"""
Started on Thu September 27 2020

Description: This file is for trying the intersection method of solving the plate waveguide dispersion equations that
             were derived using the scalar and vector potential functions.

Usage:
1) Mainly for error comparison with other dispersion curve calculation methods.

@author: Christopher Hakoda
         christopherhakoda@gmail.com

"""

import numpy as np
import elasticity as es
import time
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
import matplotlib as mpl


# Setting Global Parameters
thickness = 0.001
cL, cT = es.wave_speed_calculator(rho=2700, lamb=5.697675e10, mu=2.5947e10)
# phase velocity range of interest for the dispersion curves
cmin = 200  # m/s; approx. min phase velocity
cmax = 10000  # m/s; approx. max phase velocity
# frequency range of interest for the dispersion curves
fmin = 1  # Hz; linspace used so this isn't the real min
fmax = 20e6  # Hz; maximum frequency

error_threshold = 0.01  # rad/m; max variation


def char_eq(k=1, freq=1e6):
    """
    Parameters
    ----------
    k: Type, complex array
    DESCRIPTION. wavenumber value {default: 1 rad/m}

    freq: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    Returns
    -------
    symmetric characteristic equation and antisymmetric characteristic equation
    """

    omega = 2 * np.pi * freq  # rad/s
    h = thickness / 2
    q = np.sqrt(omega ** 2 / (cT ** 2) - k ** 2)  # rad/m
    p = np.sqrt(omega ** 2 / (cL ** 2) - k ** 2)  # rad/m

    symmetric_left = np.tan(q * h) / np.tan(p * h)
    symmetric_right = -4 * k ** 2 * p * q / ((q ** 2 - k ** 2) ** 2)

    antisymmetric_left = np.tan(q * h) / np.tan(p * h)
    antisymmetric_right = -((q ** 2 - k ** 2) ** 2) / (4 * k ** 2 * p * q)

    return symmetric_left-symmetric_right, antisymmetric_left-antisymmetric_right


def my_diff(data_array):
    # gets the point between two adjacent points and saves to a new array
    new_len = len(data_array)-1
    new_array = np.empty(new_len)
    i = 0
    while (i+1) <= (len(data_array)-1):
        new_array[i] = np.mean(np.array([data_array[i], data_array[i+1]]))
        i = i+1
    return new_array


def my_mag(data_array):
    # adds the real and imaginary part
    return np.real(data_array)+np.imag(data_array)


def my_solver(k_start=1, num_points=500, freq=1e6, k_end=100000, solution_type='symmetric'):
    """
    Parameters
    ----------
    k_start: Type, float
    DESCRIPTION. starting wavenumber used when finding roots. {default: 1 rad/m}

    k_end: Type, float
    DESCRIPTION. starting wavenumber used when finding roots. {default: 100000 rad/m}

    num_points: Type, int
    DESCRIPTION. number of points for every 1000 rad/m from which to search for the root. This is directly linked to
    accuracy of the root {default: 500 points}

    freq: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    solution_type: Type, string
    DESCRIPTION. either symmetric or antisymmetric. {default: symmetric}

    Returns
    -------
    k_refine: Type, complex float
    DESCRIPTION. the wavenumber that will be used as input for refining the root

    k_next_root: Type, complex float
    DESCRIPTION. the wavenumber that can be used as input for finding the next root.

    error: Type, complex float
    DESCRIPTION. the difference between the two outer wavenumber points used when solving for roots
    """
    if isinstance(num_points, int):
        pass
    else:
        print('Please make sure that num_points is an int.')
        return

    omega = 2*np.pi*freq
    kmax = omega/cmin  # cmin = omega/kmax
    kmin = omega/cmax  # cmax = omega/kmin

    if k_start < kmin:
        pass  # just use kmin as defined by cmax
    else:
        kmin = k_start  # as long as k0 is larger than kmin (i.e., cp<cmax) then kmin will be changed to k0

    if k_end > kmax:
        pass  # just use kmax as defined by cmin
    else:
        kmax = k_end  # as long as k0 is larger than kmin (i.e., cp<cmax) then kmin will be changed to k0

    n = np.array([0, 1, 2])
    wave_array = np.linspace(start=kmin, stop=kmax, num=num_points, dtype='complex')  # create a complex wavenumber array based on parameters
    fsym, fasym = char_eq(k=wave_array, freq=freq)  # solve the characteristic equation for the wavenumber array and a freq
    if solution_type == 'symmetric':  # both of these cannot be calculated at once; here we choose which one to use
        fk = fsym
    else:
        fk = fasym

    while n[2] <= num_points-1:  # last indice is less than or equal to the maximum number of points minus 1 because it starts from zero
        fn = fk[n]
        # looks at three adjacent points at a time
        min_boolean = (np.abs(fn[1])-np.abs(fn[0])) < 0 and (np.abs(fn[2])-np.abs(fn[1])) > 0  # if _ is negative and _ is positive
        if my_mag(fn[0])*my_mag(fn[2]) < 0 and min_boolean:  # if _ negative and _ true
            k_refine = wave_array[n[0]]
            k_next_root = wave_array[n[2]]
            k_error = wave_array[n[2]]-wave_array[n[0]]
            return k_refine, k_next_root, k_error  # stop while loop and return two k0 values which can be used for refining
        else:
            pass  # Do Nothing
        n += 1  # increment to the next set of points

    # print('Refinement Terminated!')
    return -1, -1, -1  # used for terminating refinement process


def my_refiner(k_start=1, freq=1e6, num_points=500, solution_type='symmetric'):
    """
    Parameters
    ----------
    k_start: Type, float
    DESCRIPTION. starting wavenumber used when finding roots. {default: 1 rad/m}

    num_points: Type, int
    DESCRIPTION. number of points for every 1000 rad/m from which to search for the root. This is directly linked to
    accuracy of the root {default: 500 points}

    freq: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    solution_type: Type, string
    DESCRIPTION. either symmetric or antisymmetric. {default: symmetric}

    Returns
    -------
    k_final: Type, complex float
    DESCRIPTION. the wavenumber that will be used as input for refining the root

    k_next_root: Type, complex float
    DESCRIPTION. the wavenumber that can be used as input for finding the next root.

    error: Type, complex float
    DESCRIPTION. the difference between the two outer wavenumber points used when solving for roots
    """
    omega = 2*np.pi*freq
    kmax = omega/cmin  # cmin = omega/kmax
    kmin = omega/cmax  # cmax = omega/kmin
    if k_start < kmin:
        pass  # just use kmin as defined by cmax
    else:
        kmin = k_start  # as long as k0 is larger than kmin (i.e., cp<cmax) then kmin will be changed to k0

    wavenumber_start = kmin
    wavenumber_end = kmax
    num_points = int(np.real(np.round((kmax-kmin)*num_points/1000)))  # uses the input num_points to determine the number of points needed for the range between kmin and kmax

    error = 100  # arbitrarily large number that is greater than error_threshold
    while error >= error_threshold:
        wavenumber_start, wavenumber_end, error = my_solver(k_start=wavenumber_start, num_points=num_points, freq=freq, k_end=wavenumber_end, solution_type=solution_type)
        num_points = 200
        if error == -1:  # hit the kmax limit without finding a root
            return -1, -1, -1  # escapes function and returns k_final=-1, k_next_root=-1, error=-1

    k_final = wavenumber_start
    k_next_root = wavenumber_end
    return k_final, k_next_root, error


def main_solver(freq=1e6, num_points=500, solution_type='symmetric'):
    """
    Parameters
    ----------
    num_points: Type, int
    DESCRIPTION. number of points for every 1000 rad/m from which to search for the root. This is directly linked to
    accuracy of the root {default: 500 points}

    freq: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    solution_type: Type, string
    DESCRIPTION. either symmetric or antisymmetric. {default: symmetric}

    Returns
    -------
    solution: Type, complex float
    DESCRIPTION. the roots where each row is a solution formatted: frequency, wavenumber, error

    error: Type, complex float
    DESCRIPTION. the difference between the two outer wavenumber points used when solving for roots
    """
    k_next = 1  # starting wavenumber for the next root
    error = 1  # initializing variable
    solution = np.array([[0, 0, 0]], dtype='complex')  # initializing the solution output variable
    while error != -1:
        k_final, k_next, error = my_refiner(k_start=k_next, freq=freq, num_points=num_points, solution_type=solution_type)
        if error != -1:
            add_on = np.array([[complex(freq, 0), k_final, error]])  # new results to add to 'solution'
            solution = np.append(solution, add_on, 0)

    solution = np.delete(solution, 0, 0)  # deletes the initialization row of zeros
    return solution


def saving_data(filename, data):
    """
    Parameters
    ----------
    filename: Type, string
    DESCRIPTION. the filename that the data will be saved to with suffix '.txt'.

    data: Type, Ndarray and complex
    DESCRIPTION. the data to be saved to file.

    Returns
    -------
    file formatted as [freq, real wavenumber, imag wavenumber, error]
    """
    formatted_data = np.column_stack((np.real(data[:, 0]), np.real(data[:, 1]), np.imag(data[:, 1]), np.real(data[:, 2])))
    np.savetxt(filename+".txt", formatted_data, delimiter=',')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Functions Above   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Solving Below   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Note to self: make sure num_points is high enough otherwise some solutions will not be found. This is easier to find
# when sweeping from low frequency to higher frequencies because, generally, num_points has to be higher for higher
# frequencies.


print('Calculating Roots of Characteristic Equation...')
start_time = time.time()  # starts a timer to know how long the solver runs
pool = mp.Pool(processes=None)  # starts multiprocessing pool

freq_array = np.linspace(start=fmin, stop=fmax, num=200)  # define iterable frequency array

# freezes some argurments since the pool can only handle passing in one argument
sym_func = partial(main_solver, num_points=500, solution_type='symmetric')
asym_func = partial(main_solver, num_points=500, solution_type='antisymmetric')

# performing multiprocessing
raw_sym = pool.map(func=sym_func, iterable=freq_array)
raw_asym = pool.map(func=asym_func, iterable=freq_array)
pool.close()
pool.join()

sol_sym = np.concatenate(raw_sym, axis=0)  # takes the raw output and formats it for plotting
sol_asym = np.concatenate(raw_asym, axis=0)  # takes the raw output and formats it for plotting

saving_data('symmetric_solution', sol_sym)
saving_data('antisymmetric_solution', sol_asym)

print('The solver ran for ' + str(time.time()-start_time) + ' secs.')
print('Calculations Done, Close Figure Windows to Finish Process...')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Plotting    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpl.rcParams['font.size'] = 14

f1 = plt.figure(1)
ax11, = plt.plot(np.real(sol_sym[:, 0])/(10**6), np.real(sol_sym[:, 0]*2*np.pi/sol_sym[:, 1])/1000, '.k', label='symmetric')
ax12, = plt.plot(np.real(sol_asym[:, 0])/(10**6), np.real(sol_asym[:, 0]*2*np.pi/sol_asym[:, 1])/1000, '.r', label='antisymmetric')
plt.axis([0, 20, 0, 10])
plt.legend()
plt.ylabel('Phase Velocity (km/s)')
plt.xlabel('Frequency (MHz)')

f2 = plt.figure(2)
plt.semilogy(np.absolute(sol_sym[:, 0])/(10**6), np.absolute(sol_sym[:, 2]), '.k', label='symmetric')
plt.semilogy(np.absolute(sol_asym[:, 0])/(10**6), np.absolute(sol_asym[:, 2]), '.r', label='antisymmetric')
plt.xlim([0, 20])
plt.legend()
plt.title('Errors for each root')

plt.show()
