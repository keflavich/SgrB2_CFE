# CFE

import numpy as np
from astropy.utils.console import ProgressBar
from astropy import units as u
import shutil

from cfemodel import cfelocal

nmc = 10000
cfes = np.zeros(nmc)
fbound = np.zeros(nmc)
fcce = np.zeros(nmc)

# Set up the parameters
qvir_arr = 10**np.random.normal(np.log10(1.1),0.4/1.1/np.log(10),nmc)
tview_arr = np.random.normal(0.74,0.16,nmc)
surfGMC_arr = 10.**np.random.normal(np.log10(4.1e3),1.7e3/4.1e3/np.log(10.),nmc)
beta0_arr = 10.**np.random.normal(np.log10(0.34),0.25/0.34/np.log(10.),nmc)
surfg_arr = (10.**np.random.normal(np.log10(1.e3),.5e3/1.e3/np.log(10.),nmc)*u.Msun/u.pc**2).to(u.kg/u.m**2).value
sigma_arr = 10.**np.random.normal(np.log10(10.e3),3.e3/10.e3/np.log(10.),nmc)
Omega_arr = (np.random.normal(1.8,.25,nmc) * u.Myr).to(u.s).value

mean_mass_density = (1e4 * 2.8 * u.Da / u.cm**3).to(u.kg/u.m**3)
e_mass_density = 0.5 * mean_mass_density
mass_density_arr = 10.**np.random.normal(np.log10(mean_mass_density.value),
                                         e_mass_density.value/mean_mass_density.value/np.log(10),
                                         nmc)
mean_cs = (0.53*u.km/u.s).to(u.m/u.s)
e_cs = (0.07*u.km/u.s).to(u.m/u.s)

cs_arr = np.random.normal((mean_cs.value), (e_cs.value), nmc)

# BEGIN test code
# Run the "test code"
shutil.copy('cfemodel/parameters.in','.')

# sanity check
print(cfelocal.cfelocalmod.f_cfelocal((0.03*u.M_sun/u.pc**3).to(u.kg/u.m**3).value,
                                      7000,
                                      200))
print("The first number above should be {0}".format(8.827569))
# END test code


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

    rslt = cfelocal.cfelocalmod.f_cfelocal(mass_density_arr[ii],
                                           sigma_arr[ii],
                                           cs_arr[ii])
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









