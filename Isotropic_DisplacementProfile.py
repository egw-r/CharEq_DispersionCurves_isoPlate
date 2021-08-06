"""
Started on Thu Aug 06 2021

Description: A library of common elasticity tools

@author: Christopher Hakoda

"""
import numpy as np
import elasticity as es

def symmetric_disp_multiple(omega, k, thickness, num_points, rho, lamb, mu):
    """
    Parameters
    ----------
    omega: Type, float or complex array
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}

    k: Type, complex array
    DESCRIPTION. wavenumber {rad/m}

    thickness: Type, float
    DESCRIPTION. thickness of the plate waveguide. {m}

    num_points: Type, float
    DESCRIPTION. number of points along thickness.

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    Returns
    -------
    symmetric Lamb wave displacement profiles for a given omega and k are arranged in a row: Y, UX, UY
    """
    omega.shape = (len(omega), 1)  # formatting shape
    k.shape = (len(k), 1)  # formatting shape
    cL, cT = es.wave_speed_calculator(rho=rho, lamb=lamb, mu=mu)
    h = thickness / 2
    q = np.sqrt(omega ** 2 / (cT ** 2) - k ** 2)  # rad/m
    p = np.sqrt(omega ** 2 / (cL ** 2) - k ** 2)  # rad/m
    thick_vect = np.linspace(start=-h, stop=h, num=num_points, dtype=np.complex128)
    thick_vect.shape = (1, len(thick_vect))  # formatting shape

    A2 = (k**2-q**2)*np.sin(q*h)/(2*1j*k*p*np.sin(p*h))

    UX = 1j*k*np.cos(p*thick_vect)*A2 + q*np.cos(q*thick_vect)
    UY = -p*np.sin(p*thick_vect)*A2 - 1j*k*np.sin(q*thick_vect)

    thick_array = np.tile(thick_vect, (len(omega), 1))
    return np.hstack((thick_array, UX, UY))


def antisymmetric_disp_multiple(omega, k, thickness, num_points, rho, lamb, mu):
    """
    Parameters
    ----------
    omega: Type, float or complex array
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}

    k: Type, complex array
    DESCRIPTION. wavenumber {rad/m}

    thickness: Type, float
    DESCRIPTION. thickness of the plate waveguide. {m}

    num_points: Type, float
    DESCRIPTION. number of points along thickness.

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}
    Returns
    -------
    antisymmetric Lamb wave displacement profiles for a given omega and k are arranged in a row: UX, UY
    """
    omega.shape = (len(omega), 1)  # formatting shape
    k.shape = (len(k), 1)  # formatting shape
    cL, cT = es.wave_speed_calculator(rho=rho, lamb=lamb, mu=mu)
    h = thickness / 2
    q = np.sqrt(omega ** 2 / (cT ** 2) - k ** 2)  # rad/m
    p = np.sqrt(omega ** 2 / (cL ** 2) - k ** 2)  # rad/m
    thick_vect = np.linspace(start=-h, stop=h, num=num_points, dtype=np.complex128)
    thick_vect.shape = (1, len(thick_vect))  # formatting shape

    A1 = -(k ** 2 - q ** 2) * np.cos(q * h) / (2 * 1j * k * p * np.cos(p * h))

    UX = 1j*k*np.sin(p*thick_vect)*A1 - q*np.sin(q*thick_vect)
    UY = p*np.cos(p*thick_vect)*A1 - 1j*k*np.cos(q*thick_vect)

    thick_array = np.tile(thick_vect, (len(omega), 1))
    return np.hstack((thick_array, UX, UY))


def symmetric_disp_single(omega, k, thickness, num_points, rho, lamb, mu):
    """
    Parameters
    ----------
    omega: Type, float or complex array
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}

    k: Type, complex array
    DESCRIPTION. wavenumber {rad/m}

    thickness: Type, float
    DESCRIPTION. thickness of the plate waveguide. {m}

    num_points: Type, float
    DESCRIPTION. number of points along thickness.

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}
    Returns
    -------
    antisymmetric Lamb wave displacement profiles for a given omega and k are arranged in a row: UX, UY
    """

    cL, cT = es.wave_speed_calculator(rho=rho, lamb=lamb, mu=mu)
    h = thickness / 2
    q = np.sqrt(omega ** 2 / (cT ** 2) - k ** 2)  # rad/m
    p = np.sqrt(omega ** 2 / (cL ** 2) - k ** 2)  # rad/m
    thick_vect = np.linspace(start=-h, stop=h, num=num_points, dtype=np.complex128)
    thick_vect.shape = (1, len(thick_vect))  # formatting shape

    A2 = (k ** 2 - q ** 2) * np.sin(q * h) / (2 * 1j * k * p * np.sin(p * h))

    UX = 1j * k * np.cos(p * thick_vect) * A2 + q * np.cos(q * thick_vect)  # symmetric
    UY = -p * np.sin(p * thick_vect) * A2 - 1j * k * np.sin(q * thick_vect)

    return thick_vect, UX, UY


def antisymmetric_disp_single(omega, k, thickness, num_points, rho, lamb, mu):
    """
    Parameters
    ----------
    omega: Type, float or complex array
    DESCRIPTION. Frequency multiplied by 2*pi {rad/s}

    k: Type, complex array
    DESCRIPTION. wavenumber {rad/m}

    thickness: Type, float
    DESCRIPTION. thickness of the plate waveguide. {m}

    num_points: Type, float
    DESCRIPTION. number of points along thickness.

    rho: Type, float
    DESCRIPTION. density of the plate waveguide. {kg/m^3}

    lamb: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}

    mu: Type, float
    DESCRIPTION. frequency value at which roots are solved at. {default: 1e6 Hz}
    Returns
    -------
    antisymmetric Lamb wave displacement profiles for a given omega and k are arranged in a row: UX, UY
    """

    cL, cT = es.wave_speed_calculator(rho=rho, lamb=lamb, mu=mu)
    h = thickness / 2
    q = np.sqrt(omega ** 2 / (cT ** 2) - k ** 2)  # rad/m
    p = np.sqrt(omega ** 2 / (cL ** 2) - k ** 2)  # rad/m
    thick_vect = np.linspace(start=-h, stop=h, num=num_points, dtype=np.complex128)
    thick_vect.shape = (1, len(thick_vect))  # formatting shape

    A1 = -(k ** 2 - q ** 2) * np.cos(q * h) / (2 * 1j * k * p * np.cos(p * h))

    UX = 1j*k*np.sin(p*thick_vect)*A1 - q*np.sin(q*thick_vect)
    UY = p*np.cos(p*thick_vect)*A1 - 1j*k*np.cos(q*thick_vect)

    return thick_vect, UX, UY