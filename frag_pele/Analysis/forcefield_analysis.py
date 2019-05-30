import numpy as np
import math
import matplotlib.pyplot as plt
from frag_pele.Growing.template_fragmenter import Atom, TemplateOPLS2005, Phi, Theta, Bond

NUM_AVOG = 6.022e23
DIEL_CNST = 1/(4*np.pi*1)
ELECTRON_CHARGE = 1.6022e-19  # C

MAGIC_NUMBER = 10e10/(4.184*1000)
Ein = 1  # ??
Esolv = 80.0  # ??
Evac = 1.257 * 10**(-6)  # H/m
I = 0.15
T = 1000  # K
R = 8.314  # J / K**-1 mol**-1


class ComputeNBONEnergies:
    def __init__(self, epsilon1=0, epsilon2=0, sigma1=0, sigma2=0, charge1=0, charge2=0, radnpSGB1=0, sgbnpType1=0,
                 radnpSGB2=0, sgbnpType2=0, radii=0):
        self.epsilon1 = epsilon1
        self.epsilon2 = epsilon2
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.charge1 = charge1
        self.charge2 = charge2
        self.radnpSGB1 = radnpSGB1
        self.sgbnpType1 = sgbnpType1
        self.radnpSGB2 = radnpSGB2
        self.sgbnpType2 = sgbnpType2
        self.radii = radii

    def energy_vdw(self):
        energy = float(2*np.sqrt(self.epsilon1))*float(2*np.sqrt(self.epsilon2))*(((float(self.sigma1)**12) /
                 (float(self.radii)**12)) - ((float(self.sigma2)**6)/(float(self.radii)**6)))
        return energy

    def charge_function(self):
        factor = np.sqrt((DIEL_CNST)*(NUM_AVOG)*(ELECTRON_CHARGE**2)*(MAGIC_NUMBER))
        charge1 = self.charge1 * factor
        charge2 = self.charge2 * factor
        energy = (charge1*charge2)/(1*self.radii)
        return energy

    def solv_pol_energy(self):
        K = 0.73*10**-10*(math.sqrt((2*NUM_AVOG*ELECTRON_CHARGE**2*I)/(Esolv*Evac*R*T)))
        alpha = math.sqrt(self.sgbnpType1*self.sgbnpType2)
        D = (self.radnpSGB1**2/(2*alpha)**2)
        fGB = math.sqrt(self.radnpSGB1**2 + alpha**2 * math.exp(-D))
        energy = ((1/Ein) - (math.exp(-K*fGB) / Esolv)) * ((self.charge1 * self.charge2) / fGB) \
                 - ((1/Ein) - (math.exp(-K*self.sgbnpType1)/Esolv))*(self.charge1**2/self.sgbnpType1)
        return energy


class PlotEnergies:
    def __init__(self, template, atom1, atom2):
        self.template = template
        self.atom1 = atom1
        self.atom2 = atom2
        self.nbon_energy = ComputeNBONEnergies(epsilon1=self.template.list_of_atoms[self.atom1].epsilon,
                                               epsilon2=self.template.list_of_atoms[self.atom2].epsilon,
                                               sigma1=self.template.list_of_atoms[self.atom1].sigma,
                                               sigma2=self.template.list_of_atoms[self.atom2].sigma,
                                               charge1=self.template.list_of_atoms[self.atom2].charge,
                                               charge2=self.template.list_of_atoms[self.atom2].charge,
                                               radnpSGB1=self.template.list_of_atoms[self.atom1].radnpSGB,
                                               radnpSGB2=self.template.list_of_atoms[self.atom2].radnpSGB,
                                               sgbnpType1=self.template.list_of_atoms[self.atom1].sgbnpType,
                                               sgbnpType2=self.template.list_of_atoms[self.atom2].sgbnpType)

    def plot_vdw_function(self):
        t1 = np.arange(0.01, 5.0, 0.001)
        energies = []
        for radii in t1:
            self.nbon_energy.radii = radii
            energy = self.nbon_energy.energy_vdw()
            energies.append(energy)
        return plt.plot(t1, energies, lw=2)

    def plot_charge_function(self):
        t1 = np.arange(0.01, 5.0, 0.001)
        energies = []
        for radii in t1:
            self.nbon_energy.radii = radii
            energy = self.nbon_energy.charge_function()
            energies.append(energy)
        return plt.plot(t1, energies, lw=2)

    def plot_solv_function(self):
        t1 = np.arange(0.01, 5.0, 0.001)
        energies = []
        for radii in t1:
            self.nbon_energy.radii = radii
            energy = self.nbon_energy.solv_pol_energy()
            energies.append(energy)
        return plt.plot(t1, energies, lw=2)

    def plot_nbon_function(self):
        t1 = np.arange(0.01, 5.0, 0.001)
        energies = []
        for radii in t1:
            self.nbon_energy.radii = radii
            energy = self.nbon_energy.energy_vdw()
            energy += self.nbon_energy.charge_function()
            energy += self.nbon_energy.solv_pol_energy()
            energies.append(energy)
        return plt.plot(t1, energies, lw=2)

    def show_plots(self):
        plt.figure(1)
        plt.subplot(221)
        self.plot_vdw_function()
        plt.subplot(222)
        self.plot_charge_function()
        plt.subplot(223)
        self.plot_solv_function()
        plt.subplot(224)
        self.plot_nbon_function()
        plt.show()
