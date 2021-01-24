"""
Started on Thu May 14 2020

Description: A library of common elasticity tools

@author: Christopher Hakoda

"""
import numpy as np


def wave_speed_calculator(rho=None, lamb=None, mu=None, E=None, G=None, poisson=None):
    """
    Parameters
    ----------
    rho: Type, float
    DESCRIPTION. density of the isotropic material {units: kg/m^3}

    lamb: Type, float
    DESCRIPTION. 1st Lame parameter of the isotropic material {units: Pa}

    mu: Type, float
    DESCRIPTION. 2nd Lame parameter (a.k.a. the shear modulus) of the isotropic material {units: Pa}

    G: Type, float
    DESCRIPTION. shear modulus of the isotropic material {units: Pa}

    E: Type, float
    DESCRIPTION. Young's modulus of the isotropic material {units: Pa}

    poisson: Type, float
    DESCRIPTION. poisson's ratio of the isotropic material {unitless}

    Returns
    -------
    The longitudinal wave speed, cL,  and shear wave speed, cT
    """

    if G != None:
        mu = G
    else:
        pass  # do nothing

    if rho != None:
        pass
    else:
        print('Please put in a valid density!')
        return

    if lamb != None and mu != None:
        LAMB = lamb
        MU = mu

    elif E != None and mu != None:
        LAMB = mu * (E - 2 * mu) / (3 * mu - E)
        MU = mu

    elif E != None and poisson != None:
        LAMB = E * poisson / ((1 + poisson) * (1 - 2 * poisson))
        MU = E / (2 * (1 + poisson))

    elif lamb != None and E != None:
        LAMB = lamb
        R = np.sqrt(E ** 2 + 9 * lamb ** 2 + 2 * E * lamb)
        MU = (E - 3 * lamb + R) / 4

    elif lamb != None and poisson != None:
        LAMB = lamb
        MU = lamb * (1 - 2 * poisson) / (2 * poisson)

    elif mu != None and poisson != None:
        LAMB = 2 * mu * poisson / (1 - 2 * poisson)
        MU = mu

    else:
        print('Please put in a valid combination of two(2) elastic parameters!')
        return

    return np.sqrt((LAMB+2*MU)/rho), np.sqrt(MU/rho)  # cL & cT, respectively
