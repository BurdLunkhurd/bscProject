import os
os.environ['pandeia_refdata'] = '/mnt/d/Project/pandeia_data-3.2-jwst/pandeia_data-3.2-jwst' #THIS MUST GO BEFORE YOUR IMPORT STATEMENT
os.environ['PYSYN_CDBS'] = '/mnt/d/Project/hlsp_reference-atlases_hst_multi_everything_multi_v16_sed/grp/redcat/trds'
import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import pandexo.engine.justplotit as jpi
import matplotlib.pyplot as plt
import numpy as np
import os
import pandeia.engine
from scipy.optimize import curve_fit
from astropy import units as u

# Nmask = [55,54,53,52,51]
binning = None

def DoPandexo(instrument,source):

    print("---Start PandEXO observation of " + source + "with" + instrument)

    exo_dict = jdi.load_exo_dict()
    #print(exo_dict['star']['w_unit'])

    exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                                #without having to redo the calcualtion
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 60*60*2 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                                #to a wavelength dependent noise floor solution (units are ppm)


    # exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['type'] = 'phoenix'
    exo_dict['star']['mag'] = 8.687               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 4970             #in K
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs


    exo_dict['planet']['type'] ='user'                       #tells pandexo you are uploading your own spectrum
    exo_dict['planet']['exopath'] = source
    exo_dict['planet']['w_unit'] = 'um'                      #other options include "um","nm" ,"Angs", "sec" (for phase curves)
    exo_dict['planet']['f_unit'] = 'fp/f*'               #other options are 'fp/f*'
    exo_dict['planet']['transit_duration'] = 60 * 60   #transit duration
    exo_dict['planet']['td_unit'] = 's'                      #Any unit of time in accordance with astropy.units can be added


    result = jdi.run_pandexo(exo_dict,[instrument])
    if binning == None:
        if(instrument=="MIRI LRS Slitless"):
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True)
            x=ex[0]
            y=ey[0]
            e=ee[0]
        elif(instrument == "NIRSpec G395H"):
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True)

            Nmask = []
            for i in range(0,len(ex[0])):
                if ee[0][i]>5*10**-4:
                    Nmask.append(i)
            print(Nmask)

            x=np.delete(ex,Nmask)
            y=np.delete(ey,Nmask)
            e=np.delete(ee,Nmask)
        else:
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True)
            x=ex[0]
            y=ey[0]
            e=ee[0]
    else:
        if(instrument=="MIRI LRS Slitless"):
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True)
            x=ex[0]
            y=ey[0]
            e=ee[0]
            print(x)
        elif(instrument == "NIRSpec G395H"):
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True,R=binning)
            Nmask = []
            for i in range(0,len(ex[0])):
                if ee[0][i]>5*10**-4:
                    Nmask.append(i)
            print("NMASK + {}".format(Nmask))

            x=np.delete(ex,Nmask)
            y=np.delete(ey,Nmask)
            e=np.delete(ee,Nmask)
        else:
            ex,ey,ee=jpi.jwst_1d_spec(result,model=True,R=binning)
            x=ex[0]
            y=ey[0]
            e=ee[0]

    np.savetxt('/root/Scripts/results/pandexocomparisons/'+instrument+"/"+instrument+"_"+str(binning)+"BINS_"+source.replace("/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_",""), np.array([x,y,e]))
    
    print("---PandEXO observation of " + source + "with" + instrument + "-DONE-" )
    
    return np.array([x,y,e])

def getPandexoResult(instrument,molecule):

    print("---Looking for PandEXO observation of " + molecule + " with " + instrument + " with " + str(binning) + " bins...")

    try:
        result = np.loadtxt("/root/Scripts/results/pandexocomparisons/"+instrument+"/"+instrument+"_"+str(binning)+"BINS_"+molecule+".txt")
        print("---File Found!")
    except:
        print("---No file found, generating one...")
        result = DoPandexo(instrument,"/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_"+molecule+".txt")
    return result

def getModel(molecule):
    model = np.loadtxt('/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_'+molecule+'.txt').transpose()
    return model

def plotComparison(mol1,mol2,instrument,save=False):
    print("---Comparing PandEXO observations of {} and {} with {}".format(mol1,mol2,instrument))
    xlims = {
        "MIRI LRS Slitless":    (4.7,14.7),
        "NIRSpec G395H":        (2.7,5.3),
        "NIRSpec G395M":        (2.7,5.3),
        "NIRSpec Prism":        (0.5,5.3)
    }
    mols = {
        "all":      "all",
        "noco":     "no CO",
        "noco2":    "no CO$_2$",
        "noso2":    "no SO$_2$"
            }

    OD1 = np.loadtxt('/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_'+mol1+'.txt',delimiter=' ').transpose()
    OD2 = np.loadtxt('/root/Scripts/results/PRTOutputs/v1298taub_FpFs_photochem1_'+mol2+'.txt',delimiter=' ').transpose()

    pand1 = getPandexoResult(instrument,mol1)
    pand2 = getPandexoResult(instrument,mol2)

    figure,axes=plt.subplots(2,1,figsize=(12,6))

    axes[0].set_xlim(xlims[instrument])
    axes[0].set_ylim(-0.0001,0.002)
    axes[0].set_title("Comparison spectra with {} and {}".format(mols[mol1],mols[mol2]))
    axes[0].set_xlabel("wavelength [micron]")
    axes[0].set_ylabel("Fp/Fs")
    axes[0].plot(OD1[0],OD1[1],color='firebrick',label='Emission, {}'.format(mols[mol1]))
    axes[0].errorbar(pand1[0],pand1[1],yerr=pand1[2],label='Pandexo output, {}'.format(mols[mol1]),fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8)

    axes[0].plot(OD2[0],OD2[1],color='mediumblue',label='Emission - {}'.format(mols[mol2]))
    axes[0].errorbar(pand2[0],pand2[1],yerr=pand2[2],label='Pandexo output, {}'.format(mols[mol2]),fmt = 's',color='b',markersize = 3,ecolor='royalblue',elinewidth=0.8)

    axes[0].legend()

    axes[1].set_xlim(xlims[instrument])
    axes[1].set_ylim(-0.0001,0.002)
    axes[1].set_title("differences, {} - {}".format(mols[mol2],mols[mol1]))
    axes[1].set_xlabel("wavelength [micron]")
    axes[1].set_ylabel("Fp/Fs")
    axes[1].plot(OD1[0],OD2[1]-OD1[1],color='black',label='{} - {}'.format(mols[mol2],mols[mol1]))
    axes[1].errorbar(pand1[0],pand2[1]-pand1[1],yerr=np.sqrt(pand1[2]**2+pand2[2]**2),label='Difference, pandexo output',fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8)
    axes[1].legend()


    plt.tight_layout()
    if(save==True):
        print("-Figure Saved!")
        plt.savefig('/root/Scripts/results/pandexocomparisons/'+instrument+"/"+instrument+"_"+ str(binning)+ "BINS_"+mol1+mol2+".png")
    else:
        print("-Figure shown!")
        plt.show()
    print("-----------------------------------------------------------------------------------------------------------------------")

def findDomain(model,observation):
    indexarray = []
    indexlist = []
    started = False
    Running = True
    i=0
    while Running==True:
        i+=1
        if model[1,i] >= 1*10**-4 and model[0,i]>min(observation[0]) and model[0,i]<max(observation[0]) and model[1,i] <= 5*10**-4: #0.5*10**-5
            indexlist.append(i)
            if started != True:
                print("Start Bump!")
                print("Added {}, {}".format(i,round(model[1,i],6)))
            else:
                print("Added {}, {}".format(i,round(model[1,i],6)))
            started =True

        if model[1,i] < 1*10**-4 and started:
            print("Stopped at {}, {}, because it is too small".format(i,round(model[1,i],6)))
            if len(indexlist)>=3:
                indexarray.append(indexlist)
            indexlist = []
            started = False
            
        if model[0,i+1] >= max(observation[0]):
            print("Stopped at {}, {}, because we hit the end of the spectrum".format(i,round(model[1,i],6)))
            if len(indexlist)>=3:
                indexarray.append(indexlist)
            started = False
            Running = False
        
    print("---Done scanning for bumps, found {} hits".format(len(indexarray)))
    print("----------------------------------------------------------")
    resultarray = []
    for i in range(0,len(indexarray)):
        try:
            minimum = min(indexarray[i])
            maximum = max(indexarray[i])
            resultarray.append([model[0,minimum],model[0,maximum]])
        except:
            print("NO RESULT")
            resultarray= [[5,5.3]]
        print(resultarray)
            
    return resultarray

def detectWithSigma(mol1,mol2,instrument,plot=True):

    mols = {
        "all":      "all",
        "noco":     "no CO",
        "noco2":    "no CO$_2$",
        "noso2":    "no SO$_2$"
            }

    observationin1 = getPandexoResult(instrument,mol1)
    observationin2 = getPandexoResult(instrument,mol2)
    
    observationin = np.array([observationin1[0],observationin2[1]-observationin1[1],np.sqrt(observationin2[2]**2+observationin1[2]**2)])

    modelin1 = getModel(mol1)
    modelin2 = getModel(mol2)
    modelin = np.array([modelin1[0],modelin2[1]-modelin1[1]])

    domains = findDomain(modelin,observationin)
    
    snrlist = []

    for domain in domains:
        obsMask = []
        modelMask = []
        
        for wavelength in observationin[0]:
            if wavelength<domain[0] or wavelength>domain[1]:
                obsMask.append(np.where(observationin[0]==wavelength)[0])
        
        for wavelength in modelin[0]:
            if wavelength<domain[0] or wavelength>domain[1]:
                modelMask.append(np.where(modelin[0]==wavelength)[0])

        observationx = np.delete(observationin[0],obsMask)
        observationy = np.delete(observationin[1],obsMask)
        observatione = np.delete(observationin[2],obsMask)
        observation = np.array([observationx,observationy,observatione])

        modelx = np.delete(modelin[0],modelMask)
        modely = np.delete(modelin[1],modelMask)
        model = np.array([modelx,modely])


        meansignal = np.mean(observation[1])
        sigma = np.sqrt(np.sum(observation[2]**2))/len(observation[2])
        SNR = meansignal/sigma
        snrlist.append(SNR)


        print("SNR = {} for {} with {}".format(SNR,mol2,instrument))

        if plot:
            figure, axes = plt.subplots(1)
            axes.set_title("SNR of simulated detections for {}, \n observed by {}".format(mols[mol2],instrument))
            axes.set_xlabel("Wavelength [micron]")
            axes.set_ylabel("Fp/Fs")
            axes.annotate("SNR = {}".format(round(SNR,1)),(2,10),xycoords='axes points',fontsize=12)
            axes.errorbar(observation[0],observation[1],yerr=observation[2],label='Pandexo differences',fmt = 's',color='b',markersize = 3,ecolor='royalblue',elinewidth=0.8)
            axes.errorbar(np.mean(observation[0]),meansignal,yerr = sigma,label='Mean Signal',fmt = 's',color='r',markersize = 6,ecolor='lightcoral',elinewidth=1.2)
            axes.plot([min(observation[0]),max(observation[0])],[0,0],color='black',linestyle="solid")
            axes.plot(model[0],model[1],color='black',label="Model difference",linestyle="dashed")

            axes.legend()
            plt.show()
    return snrlist

def getChi2(observation,model):
    # Make sure they are of the same length
    return np.sum((observation[1] - model[1]) ** 2 / (observation[2] ** 2))

def modelChi2(molecule,instrument):
    premodel = getModel(molecule)
    observation = getPandexoResult(instrument,molecule)

    modelx = observation[0]
    modely = np.interp(observation[0],premodel[0],premodel[1])
    model = np.array([modelx,modely])


    chi2 = getChi2(observation,model)
    return model, chi2

def linearFunc(x,a,b):
    return x*a+b

def linearFit(molecule,instrument):
    observation = getPandexoResult(instrument,molecule)
    popt, pcov = curve_fit(linearFunc,observation[0],observation[1])
    model = np.array([observation[0],linearFunc(observation[0],popt[0],popt[1])])
    chi2 = getChi2(observation,model)
    return model, popt, chi2, observation

def FitTest(molecule,instrument,plot=True,show=False):
    
    mols = {
        "all":      "all",
        "noco":     "no CO",
        "noco2":    "no CO$_2$",
        "noso2":    "no SO$_2$"
            }
    
    linearfit,params,lchi2, observation = linearFit(molecule,instrument)
    model, mchi2 = modelChi2(molecule,instrument)

    print("------------------------------------------------")
    print("Pandexo results for {} observed with {}".format(molecule,instrument))
    print(" Linear fit:")
    print("  Parameters: a= {}".format(params[0]))
    print("  Parameters: b= {}".format(params[1]))
    print("  Chi2 = {}".format(round(lchi2,3)))
    print(" Input model:")
    print("  Chi2 = {}".format(round(mchi2,3)))

    if plot == True:
        figure, axes = plt.subplots(1)

        axes.set_title('$\\chi^2$ comparison for emissions with {}, \n observed by {}'.format(mols[molecule],instrument))
        axes.set_ylabel("Fp/Fs")
        axes.set_xlabel("Wavelength [micron]")
        
        axes.annotate("Linear fit params: a={}, b={} \n    $\\chi^2$= {}\nInput model:\n    $\\chi^2$={}".format(round(params[0],5),round(params[1],5),round(lchi2,1),round(mchi2,1)),(2,10),xycoords='axes points',fontsize=12)

        axes.plot(linearfit[0],linearfit[1],label="Linear fit",color = 'blue')
        axes.errorbar(observation[0],observation[1],yerr = observation[2],label="Pandexo output",fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8)
        axes.plot(model[0],model[1],label="Input model",color = 'black')
        axes.legend()
        plt.savefig("/root/Scripts/results/chi2 test/{}/{}_{}_CHI2.png".format(instrument,instrument,molecule))
        if show:
            plt.show()

def BlackbodyFunc(wave,temp):
    wave = wave*10**-6 #from micron to meters

    h = 6.6261e-34  
    c = 3e8         
    k_B = 1.3806e-23

    exponent = h*c/(wave*k_B*temp)

    B = (2*h*c**2)/wave**5/(np.exp(exponent)-1) #unit= J/(m2 s m)

    return B

def fitBlackbody(instrument,plot=True):
    Tstar = 4970
    
    observation = getPandexoResult(instrument,"all")
    model = getModel("all")

    Is = np.array([observation[0],BlackbodyFunc(observation[0],Tstar)])

    rs = 1.33 * 6.957*10**8     #meters
    rp = 0.92 * 69911000        #meters

    Ip = np.array([observation[0],observation[1] * Is[1] * (rs**2/rp**2),observation[2] * Is[1] * (rs**2/rp**2)])

    Is = np.array([model[0],BlackbodyFunc(model[0],Tstar)])
    Imodel = np.array([model[0],model[1] * Is[1] * (rs**2/rp**2)])
    popt, pcov = curve_fit(BlackbodyFunc,Ip[0],Ip[1],bounds=(100,6000))
    fit=np.array([observation[0],BlackbodyFunc(observation[0],popt[0])])
    
    dof = len(Ip[1])-len(popt)
    chi2=getChi2(Ip,fit)
    Rchi2=chi2/dof
    residuals = Ip[1]-BlackbodyFunc(Ip[0],popt[0])
    residual_variance = np.sum(residuals**2)/dof
    perr = np.sqrt(np.diag(pcov))
    print(Rchi2)

    if plot:
        xlims= {
           "MIRI LRS Slitless":    (4.7,14.7),
            "NIRSpec G395H":        (2.7,5.3),
            "NIRSpec G395M":        (2.7,5.3),
            "NIRSpec Prism":        (0.5,5.3) 
        }
        figure,axes = plt.subplots(3,1,figsize = (10,10),gridspec_kw={'height_ratios': [2, 2, 1]})

        axes[0].set_title("Fp/Fs spectrum and observation by {}".format(instrument))
        unit = 1*u.J/u.s/u.m/u.m/u.m
        axes[0].set_ylabel("FpFs")
        axes[0].set_xlabel("wavelength (um)")
        axes[0].errorbar(observation[0],observation[1],yerr=observation[2],fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8,label='Simulated {} Detection'.format(instrument))
        axes[0].plot(model[0],model[1],linewidth=1,color="b",label="input Fp/Fs given by pRT")
        axes[0].set_xlim(xlims[instrument])
        axes[0].set_ylim(-0.0001,0.002)
        axes[0].legend()

        axes[1].set_title("Blackbody fit on simulated {} observations \nof V1298Tau b".format(instrument))
        unit = 1*u.J/u.s/u.m/u.m/u.m
        axes[1].set_ylabel("Irradiance [{}]".format(unit.unit))
        axes[1].set_xlabel("wavelength (um)")
        axes[1].errorbar(Ip[0],Ip[1],yerr=Ip[2],fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8,label='Simulated {} Detection'.format(instrument))
        axes[1].plot(fit[0],fit[1],color=(0,0,0),label='Blackbody fit')
        axes[1].plot(Imodel[0],Imodel[1],linewidth=1,color="b",label="planetary emission calculated from original FpFs")
        axes[1].set_xlim(xlims[instrument])
        axes[1].set_ylim(-4e9,4e9)
        axes[1].annotate("Fitted temperature: {}+-{}K \nReduced chi2: {}".format(round(popt[0]),round(perr[0],1),round(Rchi2,4)),(10,10),xycoords='axes points',size = 15)
        axes[1].legend()
        axes[2].set_title("Residuals")
        axes[2].errorbar(Ip[0],Ip[1]-BlackbodyFunc(Ip[0],popt[0]),yerr = Ip[2],fmt = 's',color='r',markersize = 3,ecolor='lightcoral',elinewidth=0.8)
        axes[2].plot((0,15),(0,0),color=(0,0,0))
        axes[2].set_ylabel("Irradiance [{}]".format(unit.unit))
        axes[2].set_xlabel("wavelength (um)")
        axes[2].set_xlim(xlims[instrument])
        axes[2].set_ylim(-2e9,2e9)
        plt.tight_layout()
        plt.show()
        print("T={}+-{} with Rchi2 = {} for {} with binning {}".format(popt[0],perr[0],Rchi2,instrument,binning))
    return [popt[0],perr[0],Rchi2,instrument,binning]

#------------------------------------------
def CHI2testAll():
    for instrument in ["MIRI LRS Slitless","NIRSpec G395H","NIRSpec G395M","NIRSpec Prism"]:
        for molecule in ["all","noco", "noco2", "noso2"]:
            FitTest(molecule, instrument,plot=True,show=True)
    print("")

def plotAllComparisons():
    SNRlist1 = []
    SNRlist2 = []
    for instrument in ["MIRI LRS Slitless","NIRSpec G395H","NIRSpec G395M","NIRSpec Prism"]:
        for molecule in ["noco", "noco2", "noso2"]:
            plotComparison("all", molecule, instrument,save=True)
            SNRlist1.append(detectWithSigma("all",molecule,instrument,plot=False))
        SNRlist2.append(SNRlist1)

    i=0
    j=0
    for instrument in ["MIRI LRS Slitless","NIRSpec G395H","NIRSpec G395M","NIRSpec Prism"]:
        print(instrument)
        for molecule in ["noco", "noco2", "noso2"]:
           
            print("SNR = {} for {} with {}".format(SNRlist2[i][j],molecule,instrument))
            j+=1
        i+=1

# CHI2testAll()
# plotAllComparisons()
Tprofile = []
for instrument in ["MIRI LRS Slitless","NIRSpec G395H","NIRSpec G395M","NIRSpec Prism"]:
    Tprofile.append(fitBlackbody(instrument,plot=True))
