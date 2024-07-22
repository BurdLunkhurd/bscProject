import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS import physical_constants as cst
from petitRADTRANS.physics import temperature_profile_function_guillot_global
from petitRADTRANS.plotlib import plot_radtrans_opacities

from astropy import units as u
from astropy import constants as c
from astropy.modeling import models
matplotlib.use('tkAgg')

data = np.loadtxt('/home/michel/code/data/photochem_model.csv', delimiter=',')
data = data.transpose()

# Stellar paramaters
Tstar = 4970 * u.K
starradius = 1.33 * u.Rsun
distanceToStar = 354 * u.lyr
# Planet parameters
planetradius = 0.916 * u.Rjup


def set_Radtrans(pressures,linespecieslist,linefractionlist,rayleighspecieslist,rayleighfractionlist,wavelength_boundaries):
    radtrans = Radtrans(
        pressures=pressures,
        line_species=linespecieslist,
        rayleigh_species=rayleighspecieslist,
        gas_continuum_contributors=['H2-H2', 'H2-He'],
        wavelength_boundaries=wavelength_boundaries
        
    )
    mass_fractions = {}
    for i in range(0,len(linespecieslist)):
        mass_fractions[linespecieslist[i]] = (linefractionlist[1]*np.ones(pressures.size))
    for i in range(0,len(rayleighspecieslist)):
        mass_fractions[rayleighspecieslist[i]] = (rayleighfractionlist[1]*np.ones(pressures.size))

    mean_molar_masses = 2.33 * np.ones(pressures.size)
    


    return radtrans, mass_fractions

def getBB(T,wavelengths):
    # Give temperature in Kelvin, and wavelengths as an array in Microns
    #scale=10*u.erg / (u.cm ** 2 * u.s * u.AA * u.sr)
    Blist = []
    for L in wavelengths:
        L = L * u.um
        L=L.to(u.m)
        B=2 * c.h * c.c**2 / (L**5 * (np.exp((c.h * c.c / (c.k_B * L * T)).to_value()) - 1))
        Blist.append(B.to_value())
    return Blist

# Params for planet
pressures = data[0]
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

# Not necessary with Robins model
T_eq = 600
T_int = 335 

linespecieslist = ['H2O','CO','CH4','CO2','HCN','C2H2','NH3']
linefractionlist = [H2Ofraction,COfraction,CH4fraction,CO2fraction,HCNfraction,C2H2fraction,NH3fraction]
rayleighspecieslist = ['H2','He']
rayleighfractionlist = [H2fraction,Hefraction]

wavelength_boundaries = [0.1,30]

reference_gravity = 10 ** 3.5
reference_pressure = 0.01

pressures_bar = pressures

# Create the radtrans object
radtrans,mass_fractions = set_Radtrans(pressures,linespecieslist,linefractionlist,rayleighspecieslist,rayleighfractionlist,wavelength_boundaries)

# Simulate the spectrum
wavelengths, flux, other_outputs = radtrans.calculate_flux(
        temperatures=temperatures,
        mass_fractions=mass_fractions,
        mean_molar_masses=mean_molar_masses,
        reference_gravity=reference_gravity,
        return_contribution=True
        )

wavelengths = wavelengths*10**4

def plot(temperatures,pressures_bar,wavelengths,flux,other_outputs):
    # Create axes
    fig, axes = plt.subplots(2,2,figsize = (8,8))

    # Plot temperatures
    axes[0,0].set_title("T-P profile")
    axes[0,0].plot(temperatures, pressures_bar,label='T_equ = {}, T_int= {}'.format(T_eq,T_int))
    axes[0,0].set_xlabel("T [Kelvin]")
    axes[0,0].set_ylabel("P [Bar]")
    axes[0,0].set_yscale('log')
    axes[0,0].set_ylim([np.max(pressures), np.min(pressures)])


    # Plot fluxes
    axes[0,1].set_title("Emitted flux")
    axes[0,1].plot(wavelengths,flux)
    axes[0,1].set_xscale('log')
    axes[0,1].set_xlabel('Wavelength [microns]')
    axes[0,1].set_ylabel(r'Planet flux, $F_{\lambda}$ [erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$]')

    # Plot contributions
    wlen_mu= wavelengths
    X,Y = np.meshgrid(wlen_mu,pressures_bar)
    axes[1,0].set_title("Emission contribution")
    axes[1,0].contour(X,Y,other_outputs['emission_contribution'],50,cmap=plt.cm.bone_r)
    axes[1,0].set_xscale('log')
    axes[1,0].set_yscale('log')
    axes[1,0].set_ylim([1e2,1e-7])
    axes[1,0].set_xlim([np.min(wlen_mu),np.max(wlen_mu)])
    axes[1,0].set_ylabel("P [Bar]")
    axes[1,0].set_xlabel('Wavelength [microns]')

    # Nirspec
    # axes[1].set_xlim(0.6,5)
    # axes[2].set_xlim(0.6,5)
    # MIRI
    #axes[1].set_xlim(4.9,27.9)
    #axes[2].set_xlim(4.9,27.9)
    # Label the axes
    plt.tight_layout()
    plt.show()

#plot(temperatures,pressures,wavelengths,flux,other_outputs)

#interpolation onto specific wavelengths
waverange = np.logspace(np.log10(0.1),np.log10(30),500)


# Transform the flux to units fp/fs
#planet flux comes in as erg/cm**2/s/cm, lets leave it like that
CGSFp = np.interp(waverange,wavelengths,flux) * u.erg / u.cm**2 / u.s / u.cm 

stellarflux=getBB(Tstar,waverange)
print(stellarflux[0])
SIFs = np.array(stellarflux)
SIFs = SIFs* u.J/u.m**3/u.s

CGSFs = SIFs.to(u.erg/u.cm**3/u.s)
print(CGSFs[1])

#apply sphere sizes
Fs = CGSFs * 4 * np.pi * starradius.to(u.cm)**2
Fp = CGSFp * 4 * np.pi * planetradius.to(u.cm)**2


FpFs = Fp/Fs

figure,axes = plt.subplots(2,2,figsize = (8,8))
axes[0,0].plot(waverange,CGSFp)
axes[0,0].set_xscale('log')
axes[0,0].set_xlim([0.6,25])
axes[0,0].plot(wavelengths,flux)
axes[0,0].set_xlabel("Wavelength [micron]")
axes[0,0].set_ylabel("Planet emission [ {} ]".format(CGSFp.unit))

axes[0,1].plot(waverange,CGSFs)
axes[0,1].set_xscale('log')
axes[0,1].set_xlim([0.6,25])

axes[0,1].set_xlabel("Wavelength [micron]")
axes[0,1].set_ylabel("Star emission [ {} ]".format(CGSFs.unit))

axes[1,0].plot(waverange,FpFs)
axes[1,0].set_xscale('log')
axes[1,0].set_xlim([0.6,25])
axes[1,0].set_xlabel("Wavelength [micron]")
axes[1,0].set_ylabel("Fp/Fs")

axes[1,1].plot(temperatures,pressures)
axes[1,1].set_yscale('log')
axes[1,1].set_ylim([np.max(pressures), np.min(pressures)])
axes[1,1].set_xlabel('T [K]')
axes[1,1].set_ylabel('P [bar]')
plt.show()

output_data = np.array([waverange,FpFs]).transpose()


np.savetxt("/home/michel/code/data/v1298taub_FpFs_photochem1_noso2.txt",output_data,delimiter=' ')


plot_radtrans_opacities(radtrans,
                        ['SO2','CO','CO2'],
                        temperature=700.,
                        pressure_bar=0.1)

plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-10,1e10])
plt.xlim([0.3,15.])
plt.ylabel('Opacity (cm$^2$ g$^{-1}$)')
plt.xlabel('Wavelength (micron)')
plt.legend()
plt.show()

###Get solar from claire et al