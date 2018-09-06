# CFE

import os
import numpy as np
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.table import Table
import shutil

from cfemodel import cfelocal

from parameters import (nmc, qvir_arr, tview_arr, surfGMC_arr, beta0_arr,
                        surfg_arr, sigma_arr, mass_density_arr, cs_arr
                       )

cfes = np.zeros(nmc)
fbound = np.zeros(nmc)
fcce = np.zeros(nmc)



# check to see what Mach numbers we're getting
Mach_number_local = sigma_arr / cs_arr

# BEGIN test code
# Run the "test code"
shutil.copy('cfemodel/parameters.in','.')

# sanity check
print(cfelocal.cfelocalmod.f_cfelocal((0.03*u.M_sun/u.pc**3).to(u.kg/u.m**3).value,
                                      7000,
                                      200))
print("The first number above should be {0}".format(8.827569))
# END test code

if os.path.exists("cfe_local_table.txt"):
    tbl = Table.read("cfe_local_table.txt", format='ascii.csv')
    fbound = tbl['fbound']

else:

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

    tbl = Table(data=[sigma_arr, surfg_arr, fbound], names=['sigma','surfg','fbound'])

    tbl.write("cfe_local_table.txt", format='ascii.csv')
