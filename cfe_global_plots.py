# CFE

from CFEClass import CFEClass
from CFEClass import PCFE
import numpy as np
from astropy.utils.console import ProgressBar

cfe = CFEClass()
pcfe = PCFE()

nmc = 1000
cfes = np.zeros(nmc)
fbound = np.zeros(nmc)
fcce = np.zeros(nmc)

qvir_arr = 10**np.random.normal(np.log10(1.1),0.4/1.1/np.log(10),nmc)
tview_arr = np.random.normal(0.74,0.16,nmc)
surfGMC_arr = 10.**np.random.normal(np.log10(4.1e3),1.7e3/4.1e3/np.log(10.),nmc)
beta0_arr = 10.**np.random.normal(np.log10(0.34),0.25/0.34/np.log(10.),nmc)
surfg_arr = 10.**np.random.normal(np.log10(1.e3),.5e3/1.e3/np.log(10.),nmc)*cfe.MSun/cfe.pc**2.
sigma_arr = 10.**np.random.normal(np.log10(10.e3),3.e3/10.e3/np.log(10.),nmc)
Omega_arr = np.random.normal(1.8,.25,nmc)/cfe.Myr
qT_arr = sigma_arr*Omega_arr*np.sqrt(2.*1.7)/np.pi/cfe.G/surfg_arr

for i in ProgressBar(range(nmc)):
    pcfe.sflaw = 0
    pcfe.qvir = 1.1
    pcfe.tsn = 3 
    pcfe.tview = tview_arr[i]
    pcfe.surfGMC = surfGMC_arr[i]
    pcfe.ecore = 0.5
    pcfe.beta0 = beta0_arr[i]
    pcfe.radfb = 0    
    surfg = surfg_arr[i]
    sigma = sigma_arr[i]
    Omega = Omega_arr[i]
    qT = qT_arr[i]
    gamma = cfe.f_cfe(pcfe, surfg, qT, Omega) # CFE
    cfes[i] = gamma[0]*100.
    fbound[i] = gamma[1]*100.
    fcce[i] = gamma[2]*100.

cfes_median = np.nanmedian(cfes)
cfes_errmin = cfes_median - np.nanpercentile(cfes, 16.)
cfes_errmax = np.nanpercentile(cfes, 84.) - cfes_median

print('CFE: %.2f+%.2f-%.2f' % (np.around(cfes_median,2),np.around(cfes_errmax,2),np.around(cfes_errmin,2)))

fbound_median = np.nanmedian(fbound)
fbound_errmin = fbound_median - np.nanpercentile(fbound, 16.)
fbound_errmax = np.nanpercentile(fbound, 84.) - fbound_median

print('fbound: %.2f+%.2f-%.2f' % (np.around(fbound_median,2),np.around(fbound_errmax,2),np.around(fbound_errmin,2)))

fcce_median = np.nanmedian(fcce)
fcce_errmin = fcce_median - np.nanpercentile(fcce, 16.)
fcce_errmax = np.nanpercentile(fcce, 84.) - fcce_median

print('fcce: %.2f+%.2f-%.2f' % (np.around(fcce_median,2),np.around(fcce_errmax,2),np.around(fcce_errmin,2)))









