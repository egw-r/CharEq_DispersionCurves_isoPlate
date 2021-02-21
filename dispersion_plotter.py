"""
Started on Saturday January 23 2021

Description: This file is for plotting data collected in a single csv file

Usage:
1) Mainly for reviewing dispersion curve data that was calculated

@author: Christopher Hakoda
         christopherhakoda@gmail.com

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def phase_plot(sym_filepath, asym_filepath, fig_num):
    """
    Parameters
    ----------
    sym_filepath: Type, string
    DESCRIPTION. the filepath for the symmetric Lamb dispersion solutions to be plotted.

    asym_filepath: Type, string
    DESCRIPTION. the filepath for the symmetric Lamb dispersion solutions to be plotted.

    fig_num: Type, integer
    DESCRIPTION. the starting figure number that will be used for plotting.

    Returns
    -------

    """

    # this file is assumed to be comma-separated and formatted as [freq,real(wavenumber),imag(wavenumber),...whatever else]
    data1 = np.loadtxt(sym_filepath, delimiter=',')
    data2 = np.loadtxt(asym_filepath, delimiter=',')

    mpl.rcParams['font.size'] = 14

    # only plots using the first three columns assuming the three columns are freq, real(k), imag(k)
    f1 = plt.figure(fig_num)
    plt.plot(data1[:, 1] / (2 * np.pi), data1[:, 0] / (10 ** 6), '.k', label='symmetric')
    plt.plot(data2[:, 1] / (2 * np.pi), data2[:, 0] / (10 ** 6), '.r', label='antisymmetric')
    plt.axis([0, None, 0, 20])
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Wavenumber (1/m)')
    plt.legend()

    f2 = plt.figure(fig_num+1)
    plt.plot(data1[:, 0] / (10 ** 6), data1[:, 0] * 2 * np.pi / data1[:, 1] / 1000, '.k', label='symmetric')
    plt.plot(data2[:, 0] / (10 ** 6), data2[:, 0] * 2 * np.pi / data2[:, 1] / 1000, '.r', label='antisymmetric')
    plt.axis([0, 20, 0, 10])
    plt.ylabel('Phase Velocity (km/s)')
    plt.xlabel('Frequency (MHz)')
    plt.legend()

    plt.show()


def phase_and_group_plot(sym_filepath, asym_filepath, fig_num):
    """
    Parameters
    ----------
    sym_filepath: Type, string
    DESCRIPTION. the filepath for the symmetric Lamb dispersion solutions to be plotted.

    asym_filepath: Type, string
    DESCRIPTION. the filepath for the symmetric Lamb dispersion solutions to be plotted.

    fig_num: Type, integer
    DESCRIPTION. the starting figure number that will be used for plotting.

    Returns
    -------

    """

    # this file is assumed to be comma-separated and formatted as [freq,real(wavenumber),imag(wavenumber),...whatever else]
    data1 = np.loadtxt(sym_filepath, delimiter=',')
    data2 = np.loadtxt(asym_filepath, delimiter=',')

    mpl.rcParams['font.size'] = 14

    # only plots using the first three columns assuming the three columns are freq, real(k), imag(k)
    f1 = plt.figure(fig_num)
    plt.plot(data1[:, 1] / (2 * np.pi), data1[:, 0] / (10 ** 6), '.k', label='symmetric')
    plt.plot(data2[:, 1] / (2 * np.pi), data2[:, 0] / (10 ** 6), '.r', label='antisymmetric')
    plt.axis([0, None, 0, 20])
    plt.ylabel('Frequency (MHz)')
    plt.xlabel('Wavenumber (1/m)')
    plt.legend()

    f2 = plt.figure(fig_num+1)
    plt.plot(data1[:, 0] / (10 ** 6), data1[:, 0] * 2 * np.pi / data1[:, 1] / 1000, '.k', label='symmetric')
    plt.plot(data2[:, 0] / (10 ** 6), data2[:, 0] * 2 * np.pi / data2[:, 1] / 1000, '.r', label='antisymmetric')
    plt.axis([0, 20, 0, 10])
    plt.ylabel('Phase Velocity (km/s)')
    plt.xlabel('Frequency (MHz)')
    plt.legend()

    f3 = plt.figure(fig_num+2)
    plt.plot(data1[:, 0] / (10 ** 6), abs(data1[:, 3]) / 1000, '.k', label='symmetric')
    plt.plot(data2[:, 0] / (10 ** 6), abs(data2[:, 3]) / 1000, '.r', label='antisymmetric')
    plt.axis([0, 20, 0, 7])
    plt.ylabel('Group Velocity (km/s)')
    plt.xlabel('Frequency (MHz)')
    plt.legend()


    plt.show()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Functions Above   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Solving Below   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

phase_plot(sym_filepath="symmetric_solution.txt", asym_filepath="antisymmetric_solution.txt", fig_num=1)
phase_and_group_plot(sym_filepath="symmetric_solution_addGroupVelocity.txt", asym_filepath="antisymmetric_solution_addGroupVelocity.txt", fig_num=3)

# plt.savefig('exampleSVG', format="svg")  # for saving the figures as an svg file
