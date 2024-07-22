import numpy as np
import matplotlib.pyplot as plt
import astropy.modeling.physical_models as models
from astropy import units as u

# Photochemical model of V1298 Tau b, made using 1D PAC by Robin Baeyens (June 2024)
#pressure_bar, temperature_K, mean_molecular_mass_amu, H2, He, H2O, CO, CH4, CO2, HCN, C2H2, NH3, SO2

data = np.loadtxt("/root/Scripts/Models/photochem_model.csv",delimiter=',').transpose()

pressures = data[0]
print(pressures)
temperatures = data[1]
mean_molar_masses = data[2]
H2fraction = data[3]
Hefraction = data[4]
H2Ofraction = data[5]
COfraction = data[6]
CH4fraction = data[7]
CO2fraction = data[8]
HCNfraction = data[9]
C2H2fraction = data[10]
NH3fraction = data[11]
SO2fraction = data[12]

figure,axes = plt.subplots(1)
axes.plot(COfraction,pressures,label='CO')
axes.plot(CO2fraction,pressures,label="CO2")
axes.plot(SO2fraction,pressures,label="SO2")
axes.plot(H2fraction,pressures,label="H2",linestyle="dashed")
axes.plot(Hefraction,pressures,label="He",linestyle="dashed")
axes.plot(H2Ofraction,pressures,label="H2O",linestyle="dashed")

axes.plot(CH4fraction,pressures,label="CH4",linestyle="dashed")

axes.plot(HCNfraction,pressures,label="HCN",linestyle="dashed")
axes.plot(C2H2fraction,pressures,label="C2H2",linestyle="dashed")
axes.plot(NH3fraction,pressures,label="NH3",linestyle="dashed")


axes.legend()
axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlim(10**-12,10**-1.5)
axes.set_ylim(10**2,10**-6)
axes.set_ylabel("Pressure [bar]")
axes.set_xlabel("Molar Fraction")
axes.set_title("Simulated abundance of CO, CO2 and SO2 in \n the atmosphere of V1298tau b")
plt.show()

figure,axes = plt.subplots(1)
axes.plot(temperatures,pressures)
axes.set_yscale('log')
axes.set_ylim(10**2,10**-7)
axes.set_ylabel("Pressure [bar]")
axes.set_xlabel("Temperature [K]")
axes.set_title("Simulated T-P profile of V1298tau b")

plt.show()