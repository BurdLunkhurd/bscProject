import numpy as np
import matplotlib.pyplot as plt

specall = np.loadtxt("/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_all.txt",delimiter=" ").transpose()
specco = np.loadtxt("/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_noco.txt",delimiter=" ").transpose()
specco2 = np.loadtxt("/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_noco2.txt",delimiter=" ").transpose()
specso2 = np.loadtxt("/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_noso2.txt",delimiter=" ").transpose()

Tstar = 4970

def BlackbodyFunc(wave,temp):
    wave = wave*10**-6 #from micron to meters

    h = 6.6261e-34  
    c = 3e8         
    k_B = 1.3806e-23

    exponent = h*c/(wave*k_B*temp)

    B = (2*h*c**2)/wave**5/(np.exp(exponent)-1) #unit= J/(m2 s m)

    return B


figure, axes = plt.subplots(2,1)

axes[0].plot(specco[0],specco[1]/BlackbodyFunc(specco[0],Tstar),label="No CO",linewidth = .9,color = "blue")
axes[0].plot(specco2[0],specco2[1]/BlackbodyFunc(specco2[0],Tstar),label="No CO2",linewidth = .9,color = "orange")
axes[0].plot(specso2[0],specso2[1]/BlackbodyFunc(specso2[0],Tstar),label="No SO2",linewidth = .9,color = "green")
axes[0].plot(specall[0],specall[1]/BlackbodyFunc(specall[0],Tstar),label="All species",color = "black")

axes[0].set_title(r"Simulated $F_p/F_s$, pRT output")
axes[0].set_xlabel("Wavelength [micron]")
axes[0].set_ylabel(r"$F_p/F_s$")
axes[0].set_xlim(0.5,15.5)

axes[1].plot(specco[0],specco[1],label="No CO",linewidth = .9,color = "blue")
axes[1].plot(specco2[0],specco2[1],label="No CO2",linewidth = .9,color = "orange")
axes[1].plot(specso2[0],specso2[1],label="No SO2",linewidth = .9,color = "green")
axes[1].plot(specall[0],specall[1],label="All species",color = "black")

axes[1].set_title(r"Simulated $F_p$, pRT output")
axes[1].set_xlabel("Wavelength [micron]")
axes[1].set_ylabel(r"$F_p$ [erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$]")
axes[1].set_xlim(0.5,15.5)


plt.legend()
plt.show()