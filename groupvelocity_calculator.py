"""
Started on Thu September 27 2020

Description: This file is for using the displacement profile of Lamb waves to calculate the group velocity for a given
frequency-wavenumber pair

Usage:
1) Mainly for error comparison with other dispersion curve calculation methods.

@author: Christopher Hakoda
         christopherhakoda@gmail.com

"""

import numpy as np
import elasticity as es
import Isotropic_DisplacementProfile as disp
import time
import matplotlib.pyplot as plt
import matplotlib as mpl

# Setting Global Parameters
num_points = 101  # number of points to integrate over; must be odd


def main():
    mat_thick = 0.001  # m
    mat_rho = 2700  # kg/m^3
    mat_lamb = 5.697675e10  # Pa
    mat_mu = 2.5947e10  # Pa

    # this file is assumed to be comma-separated and formatted as [freq,real(wavenumber),imag(wavenumber),...whatever else]
    filepath1 = "symmetric_solution.txt"
    filepath2 = "antisymmetric_solution.txt"

    print('Calculating the group velocity...')
    start_time = time.time()  # starts a timer to know how long the solver runs

    data1 = np.loadtxt(filepath1, delimiter=',')
    data2 = np.loadtxt(filepath2, delimiter=',')

    omega1 = data1[:, 0] * 2 * np.pi
    omega1 = omega1.astype(np.complex128)  # change imported data to complex type
    k1 = data1[:, 1] + 1j * data1[:, 2]
    k1 = k1.astype(np.complex128)  # change imported data to complex type

    omega2 = data2[:, 0] * 2 * np.pi
    omega2 = omega2.astype(np.complex128)  # change imported data to complex type
    k2 = data2[:, 1] + 1j * data2[:, 2]
    k2 = k2.astype(np.complex128)  # change imported data to complex type

    sym_disp = disp.symmetric_disp_multiple(omega=omega1, k=k1, thickness=mat_thick, num_points=num_points, rho=mat_rho, lamb=mat_lamb, mu=mat_mu)
    asym_disp = disp.antisymmetric_disp_multiple(omega=omega2, k=k2, thickness=mat_thick, num_points=num_points, rho=mat_rho, lamb=mat_lamb,
                                            mu=mat_mu)

    sym_result = isotropic_multiple(disp=sym_disp, omega=omega1, k=k1, rho=mat_rho, lamb=mat_lamb, mu=mat_mu)
    asym_result = isotropic_multiple(disp=asym_disp, omega=omega2, k=k2, rho=mat_rho, lamb=mat_lamb, mu=mat_mu)

    sym_result[:, 0] = sym_result[:, 0] / (2 * np.pi)  # changing from omega to frequency
    asym_result[:, 0] = asym_result[:, 0] / (2 * np.pi)

    # uses the filepath name as a basis for saving new data; this does not save the previous error data and replaces it
    # with the group velocity in the fourth column.
    saving_data(filepath1[:-4], sym_result)
    saving_data(filepath2[:-4], asym_result)

    print('The solver ran for ' + str(time.time() - start_time) + ' secs.')
    print('Calculations Done, Close Figure Windows to Finish Process...')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Plotting    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mpl.rcParams['font.size'] = 14

    f1 = plt.figure(1)
    plt.plot(np.real(sym_result[:, 0]) / (10 ** 6), np.real(sym_result[:, 0] * 2 * np.pi / sym_result[:, 1]) / 1000,
             '.k', label='symmetric')
    plt.plot(np.real(asym_result[:, 0]) / (10 ** 6), np.real(asym_result[:, 0] * 2 * np.pi / asym_result[:, 1]) / 1000,
             '.r', label='antisymmetric')
    plt.axis([0, 20, 0, 10])
    plt.legend()
    plt.ylabel('Phase Velocity (km/s)')
    plt.xlabel('Frequency (MHz)')

    f2 = plt.figure(2)
    # absolute value applied here because the group velocity is negative due to the original wave propagation assumption
    # used in the derivation of the displacement profile expression being a backward propagating wave. In other words the
    # absolute value is purely aesthetic feel free to remove it and adjust axis for a more accurate representation.
    ax11, = plt.plot(np.real(sym_result[:, 0]) / (10 ** 6), np.abs(np.real(sym_result[:, 2])) / 1000, '.k',
                     label='symmetric')
    ax12, = plt.plot(np.real(asym_result[:, 0]) / (10 ** 6), np.abs(np.real(asym_result[:, 2])) / 1000, '.r',
                     label='antisymmetric')
    plt.axis([0, 20, 0, 6])
    plt.legend()
    plt.ylabel('Group Velocity (km/s)')
    plt.xlabel('Frequency (MHz)')

    plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Functions Below  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def isotropic_multiple(disp, omega, k, rho, lamb, mu):
    """
    Parameters
    ----------
    disp: Type, complex
    DESCRIPTION. complex NxM array that is formatted as (thick_vect, UX, UY) where each thick_vect, UX and UY at a given
    frequency-wavenumber pair are 1xM. These are the displacement profiles for a Lamb wave. thick_vect is in {m}.

    omega: Type, complex
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}, wavenumber {rad/m}

    k: Type, complex
    DESCRIPTION. wavenumber {rad/m}

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    Returns
    -------
    complex array that has original omega and k data array but with an extra column appended that has the group velocity
    """

    omega.shape = (len(omega), 1)  # formatting shape
    k.shape = (len(k), 1)  # formatting shape

    M, N = np.shape(disp)
    N_0 = int(N/3)
    Y = disp[:, 0:N_0]
    UX = disp[:, N_0:(2*N_0)]
    UY = disp[:, (2*N_0):(3*N_0)]

    # note to self: the '::' notation ensures that even the last element is included in the range starting from the
    # first number at a step of the second number.
    y = Y[:, 1::2]  # choosing every other value starting from the point between the first and third point
    dy = np.absolute(y[0, 1]-y[0, 0])  # assumes 'y' is a 1xN array after being indexed
    ux = UX[:, 1::2]  # choosing every other value starting from the point between the first and third point
    Dux = np.diff(UX[:, 0::2], axis=1)/dy  # calculating the numerical derivative, by getting the difference between the first and third point
    uy = UY[:, 1::2]
    Duy = np.diff(UY[:, 0::2], axis=1)/dy  # calculating the numerical derivative

    # note that a positive wavenumber and a positive frequency corresponds with a backwards propagating wave so the
    # expression for group velocity below was adjusted to take this into account. The real and imaginary parts of the
    # expression were integrated separately.
    cgx_num = -2*k*(mu*uy*np.conj(uy)+(lamb+2*mu)*ux*np.conj(ux))+1j*(lamb*Duy*np.conj(ux)-lamb*ux*np.conj(Duy)+mu*Dux*np.conj(uy)-mu*uy*np.conj(Dux))
    cgx_num_r = np.trapz(y=np.real(cgx_num), x=y)
    cgx_num_i = np.trapz(y=np.imag(cgx_num), x=y)
    cgx_den = np.trapz(y=(2*omega*rho*(ux*np.conj(ux)+uy*np.conj(uy))), x=y)  # guaranteed real due to dot product with conjugate
    cgx = (cgx_num_r+1j*cgx_num_i)/cgx_den
    cgx.shape = (len(cgx), 1)
    return np.hstack((omega, k, cgx))


def isotropic_single(omega, k, thickness, rho, lamb, mu, solution_type):
    """
    Parameters
    ----------
    omega: Type, complex
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}

    k: Type, complex
    DESCRIPTION. wavenumber {rad/m}

    thickness: Type, real
    DESCRIPTION. thickness of the plate waveguide {m}

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at.

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at.

    solution_type: Type, string
    DESCRIPTION. either 'symmetric' or 'antisymmetric'.

    Returns
    -------
    complex array that has original omega and k data array but with an extra column appended that has the group velocity
    """

    # ~~~~~~~~~~~~~~~~ displacement profile calculation ~~~~~~~~~~~~
    if solution_type == 'symmetric':
        thick_vect, UX, UY = disp.symmetric_disp_single(omega=omega, k=k, thickness=thickness, rho=rho, lamb=lamb, mu=mu)
    else:
        thick_vect, UX, UY = disp.antisymmetric_disp_single(omega=omega, k=k, thickness=thickness, rho=rho, lamb=lamb, mu=mu)

    # ~~~~~~~~~~~~~~~~ group velocity calculation that uses the displacement profile ~~~~~~~~~~~~
    y = thick_vect[0, 1::2]  # choosing every other value starting from the point between the first and third point
    dy = np.absolute(y[1] - y[0])  # assumes 'y' is a 1xN array after being indexed
    ux = UX[0, 1::2]  # choosing every other value starting from the point between the first and third point
    Dux = np.diff(UX[0, 0::2]) / dy  # calculating the numerical derivative, by getting the difference between the first and third point
    uy = UY[0, 1::2]
    Duy = np.diff(UY[0, 0::2]) / dy  # calculating the numerical derivative

    cgx_num = -2*k*(mu*uy*np.conj(uy)+(lamb+2*mu)*ux*np.conj(ux))+1j*(lamb*Duy*np.conj(ux)-lamb*ux*np.conj(Duy)+mu*Dux*np.conj(uy)-mu*uy*np.conj(Dux))
    cgx_num_r = np.trapz(y=np.real(cgx_num), x=y)
    cgx_num_i = np.trapz(y=np.imag(cgx_num), x=y)
    cgx_den = np.trapz(y=2*omega*rho*(ux*np.conj(ux)+uy*np.conj(uy)), x=y)
    cgx = (cgx_num_r + 1j * cgx_num_i) / cgx_den
    print('For a '+solution_type+' Lamb wave mode at a frequency of '+str(omega / (2 * np.pi))+' Hz, the group velocity was calculated to be: '+str(cgx)+' m/s')
    if cgx < 0:
        print('Note: a positive wavenumber and a positive frequency corresponds with a backwards propagating wave which results in a negative group velocity.')
    # ~~~~~~~~~~~~~~~~ plotting displacement profiles ~~~~~~~~~~~~
    f1 = plt.figure(1)
    plt.plot(np.real(thick_vect[0, :]), np.real(UX)[0, :], '-k', label='UX')
    plt.plot(np.real(thick_vect[0, :]), np.imag(UX[0, :]), '-.k', label='UX')
    plt.plot(np.real(thick_vect[0, :]), np.real(UY[0, :]), '-r', label='UY')
    plt.plot(np.real(thick_vect[0, :]), np.imag(UY[0, :]), '-.r', label='UY')
    plt.legend()

    plt.show()


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
    file formatted as [omega,real wavenumber, imag wavenumber, group velocity]
    """
    formatted_data = np.column_stack((np.real(data[:, 0]), np.real(data[:, 1]), np.imag(data[:, 1]), np.real(data[:, 2])))
    np.savetxt(filename+"_addGroupVelocity.txt", formatted_data, delimiter=',')


if __name__ == "__main__":
    main()
