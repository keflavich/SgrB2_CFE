# CFE

import numpy as np
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import constants

from cfemodel import cfeglobal

from parameters import (nmc, qvir_arr, tview_arr, surfGMC_arr, beta0_arr,
                        surfg_arr, sigma_arr, Omega_arr, qT_arr,
                        mass_density_arr, cs_arr
                       )


cfes = np.zeros(nmc)
fbound = np.zeros(nmc)
fcce = np.zeros(nmc)


# check to see what Mach numbers we're getting
phi_p = 10 - 8*(1+0.025*(surfg_arr/(100*u.Msun.si.scale/u.pc.si.scale**2))**-2)**-1
#Q_arr = 2**0.5 * Omega_arr * sigma_arr / (np.pi * cfe.G * surfg_arr)
Mach_number_global = 2.82 * phi_p**(1/8.) * qT_arr * (Omega_arr/u.Myr.si.scale**-1)**-1 * (surfg_arr/(100*u.Msun.si.scale/u.pc.si.scale**2))

#from CFEClass import CFEClass
#from CFEClass import PCFE
#cfe = CFEClass()
#pcfe = PCFE()
# for i in ProgressBar(range(nmc)):
#     pcfe.sflaw = 0
#     pcfe.qvir = 1.1
#     pcfe.tsn = 3 
#     pcfe.tview = tview_arr[i]
#     pcfe.surfGMC = surfGMC_arr[i]
#     pcfe.ecore = 0.5
#     pcfe.beta0 = beta0_arr[i]
#     pcfe.radfb = 0    
#     surfg = surfg_arr[i]
#     sigma = sigma_arr[i]
#     Omega = Omega_arr[i]
#     qT = qT_arr[i]
#     gamma = cfe.f_cfe(pcfe, surfg, qT, Omega) # CFE
#     cfes[i] = gamma[0]*100.
#     fbound[i] = gamma[1]*100.
#     fcce[i] = gamma[2]*100.
# 
# cfes_median = np.nanmedian(cfes)
# cfes_errmin = cfes_median - np.nanpercentile(cfes, 16.)
# cfes_errmax = np.nanpercentile(cfes, 84.) - cfes_median
# 
# print('CFE: %.2f+%.2f-%.2f' % (np.around(cfes_median,2),np.around(cfes_errmax,2),np.around(cfes_errmin,2)))
# 
# fbound_median = np.nanmedian(fbound)
# fbound_errmin = fbound_median - np.nanpercentile(fbound, 16.)
# fbound_errmax = np.nanpercentile(fbound, 84.) - fbound_median
# 
# print('fbound: %.2f+%.2f-%.2f' % (np.around(fbound_median,2),np.around(fbound_errmax,2),np.around(fbound_errmin,2)))
# 
# fcce_median = np.nanmedian(fcce)
# fcce_errmin = fcce_median - np.nanpercentile(fcce, 16.)
# fcce_errmax = np.nanpercentile(fcce, 84.) - fcce_median
# 
# print('fcce: %.2f+%.2f-%.2f' % (np.around(fcce_median,2),np.around(fcce_errmax,2),np.around(fcce_errmin,2)))





for ii in ProgressBar(range(nmc)):
    # we write the parameters to file for each call, since the fortran code
    # reads these parameters from file
    with open('parameters.in', 'w') as fh:
        fh.write("""
0                star formation law
{qvir}              GMC virial parameter
{tsn}                time of first SN (Myr)
{tview}               time of determining the CFE (Myr)
{gmcsurfdens}              GMC surface density (Msun/pc^2)
{sfemax}              maximum (protostellar core) SFE
{beta0}            turbulent/magnetic pressure ratio
0                SN/radiative feedback mode
""".strip().lstrip().format(qvir=qvir_arr[ii],
                            tview=tview_arr[ii],
                            tsn=3,
                            gmcsurfdens=surfGMC_arr[ii],
                            sfemax=0.5,
                            beta0=beta0_arr[ii],
        ))

    rslt = cfeglobal.cfemod.f_cfe(surfg_arr[ii],
                                        qT_arr[ii],
                                        Omega_arr[ii])
    cfe_,fbound_,fcce_,fcce2 = rslt


    surfg = surfg_arr[ii]
    sigma = sigma_arr[ii]
    Omega = Omega_arr[ii]

    cfes[ii] = cfe_*100
    fbound[ii] = fbound_*100
    fcce[ii] = fcce_*100

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
