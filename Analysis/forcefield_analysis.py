import numpy as np
import matplotlib.pyplot as plt


NUM_AVOG = 6.022e23
DIEL_CNST = 1/(4*np.pi*1)
ELECTRON_CHARGE = 1.6022e-19
MAGIC_NUMBER = 10e10/(4.184*1000)


def energy_vdw(epsilon1, epsilon2, sigma1, sigma2, radius):
    energy = float(2*np.sqrt(epsilon1))*float(2*np.sqrt(epsilon2))*(((float(sigma1)**12)/(float(radius)**12)) -
                                              ((float(sigma2)**6)/(float(radius)**6)))
    return energy


def charge_function(charge1, charge2, radius):
    factor = np.sqrt((DIEL_CNST)*(NUM_AVOG)*(ELECTRON_CHARGE**2)*(MAGIC_NUMBER))
    charge1 = charge1 * factor
    charge2 = charge2 * factor
    energy = (charge1*charge2)/(1*radius)
    return energy


def plot_vdw_function(epsilon1, epsilon2, sigma1, sigma2):
    t1 = np.arange(0.01, 5.0, 0.001)
    energies = []
    for radii in t1:
        energy = energy_vdw(epsilon1=epsilon1, epsilon2=epsilon2, sigma1=sigma1, sigma2=sigma2, radius=radii)
        energies.append(energy)
    return plt.plot(t1, energies, lw=2)


def plot_nbon_function(epsilon1, epsilon2, sigma1, sigma2, charge1, charge2):
    t1 = np.arange(0.01, 5.0, 0.001)
    energies = []
    for radii in t1:
        energy = energy_vdw(epsilon1=epsilon1, epsilon2=epsilon2, sigma1=sigma1, sigma2=sigma2, radius=radii)
        energy += charge_function(charge1=charge1, charge2=charge2, radius=radii)
        energies.append(energy)
    return plt.plot(t1, energies, lw=2)


def plot_charge_function(charge1, charge2):
    t1 = np.arange(0.01, 5.0, 0.001)
    energies = []
    for radii in t1:
        energy = charge_function(charge1=charge1, charge2=charge2, radius=radii)
        energies.append(energy)
    print(energies)
    return plt.plot(t1, energies, lw=2)


def main(epsilon1, epsilon2, sigma1, sigma2, charge1, charge2):
    plt.figure(1)
    plt.subplot(221)
    plot_vdw_function(epsilon1, epsilon2, sigma1, sigma2)
    plt.ylim(-0.1, 0.001)
    plt.xlim(0.0, 5.0)
    plt.subplot(222)
    plot_charge_function(charge1, charge2)
    plt.subplot(223)
    plot_nbon_function(epsilon1, epsilon2, sigma1, sigma2, charge1, charge2)
    plt.ylim(-0.1, 0.001)
    plt.xlim(0.0, 5.0)
    plt.show()


main(epsilon1=0.066, epsilon2=0.066, sigma1=3.5, sigma2=3.5, charge1=-0.50, charge2=-0.50)