# bscProject
Scripts used for my bachelors project

Files:

photochem_model.csv
-- This file contains the atmospheric parameters needed to run the pRT emission spectrum simulations. THis was provided to me by the group 
Column headers: pressures, temperatures, mean_molar_mass, H2fraction, Hefraction, H2Ofraction, COfraction, CH4fraction, CO2fraction, HCNfraction, C2H2fraction, NH3fraction, SO2fraction


Scripts

singleEmission Import.py
-- Import the provided atmosphere model and generate a planetary emission spectrum using petitRADTRANS, for all chemical species present in the provided model. Saves the data in a .txt

singleEmission Import no co.py
-- Import the provided atmosphere model and generate a planetary emission spectrum using petitRADTRANS, but not taking into account CO. Saves the data in a .txt

singleEmission Import no co2.py
-- Import the provided atmosphere model and generate a planetary emission spectrum using petitRADTRANS, but not taking into account CO2. Saves the data in a .txt

singleEmission Import no so2.py
-- Import the provided atmosphere model and generate a planetary emission spectrum using petitRADTRANS, but not taking into account SO2. Saves the data in a .txt

pandexoComparisonGeneral2.py
-- Contains functions:
DoPandexo(instrument,source): 
provide instrument name and path to FpFs file, 
return pandexo result. 
It will simulate a JWST observation of that spectrum, and also save it to a file.

getPandexoResult(instrument,molecule): 
provide an instrument and the type of spectrum you want, 
return pandexo results, 
it will try to look for it in the directory, if it cant find it, it will call doPandexo

getModel(molecule): 
Provide the type of spectrum you want ("all","noco","noco2","noso2"), 
return pRT simualated emission spectrum.

plotComparison(mol1,mol2,instrument,save=False): 
provide two types of spectrum ("all","noco","noco2","noso2"), 
return shows or saves a figure depending on the state of save. 
it calls getPandexoResult for both spectra, and subtracts them. Then it plots them.

findDomain(model,observation):
provide: compared prt emission spectrum, compared pandexo result
return: a 2d array containing lists of indexes 
It will search for indexes of the pandexo results where the compared emission spectrum is above a certain threshold. it discards bumps that are shorter than 3 indexes.

detectWithSigma(mol1,mol2,instrument,plot=True):
provide: 2 types of spectra ("all","noco","noco2","noso2"), a JWST instrument.
return a list of snr's. one for each bump in the compared emission spectrum
looks up emission spectra and pandexo results for the two specified spectrum types. Subtracts them, calls findDomain, and analyses the snr on the given domains.

getChi2(observation,model):
provide: 1 types of spectra ("all","noco","noco2","noso2"), pandexo result,
return: chi2 value between a model and observation

modelChi2(molecule,instrument):
provide: 1 type of spectra ("all","noco","noco2","noso2"), a JWST instrument
return: returns chi2 between the prt emission spectrum and its pandexo result
calls getModel for the emission spectrum, and calls getPandexoResult for the pandexo result. then sends them through getChi2

linearFunc(x,a,b):
provide 1 variable and 2 parameters
return y
just a function to fit with

linearFit(molecule,instrument):
provide: type of spectrum ("all","noco","noco2","noso2"), a JWST instrument
return: array for x and y of fitted line, optimal parameters, chi2 for fit compared to input pandexo result, 

FItTest(molecule,instrument,plot=True,show=False):
provide: type of spectrum ("all","noco","noco2","noso2"), a JWST instrument
return: chi2s and plot for a comparison of linear fit and input model through the pandexo result for a given molecule.
calls linearFit, then calls modelChi2




