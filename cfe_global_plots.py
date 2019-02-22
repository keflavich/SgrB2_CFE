# CFE
import os

import numpy as np
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import constants
from astropy import table
from astropy.table import Table
import textwrap
import shutil

from cfemodel import cfeglobal

from parameters import (nmc, qvir_arr, tview_arr, surfGMC_arr, beta0_arr,
                        surfg_arr, sigma_arr, Omega_arr, qT_arr,
                        mass_density_arr, cs_arr
                       )


# basic test
shutil.copy('cfemodel/parameters.in', 'parameters.in')
result1 = cfeglobal.cfemod.f_xcce(surfg=(9.3*u.M_sun/u.pc**2).to(u.kg/u.m**2).value,
                                  qt=2.0, omega=(0.026/u.Myr).to(u.s**-1).value,
                                  surfgmc=(100*u.M_sun/u.pc**2).to(u.kg/u.m**2).value,
                                  qvir=1.3, tview=(10*u.Myr).to(u.s).value
                                 )
result = cfeglobal.cfemod.f_cfe(
    surfg=(9.3*u.M_sun/u.pc**2).to(u.kg/u.m**2).value, qt=2.0,
    omega=(0.026/u.Myr).to(u.s**-1).value,
)
# if either of these fail, it indicates that the code was not properly compiled
assert np.abs(result1 - 138.797882) / 138.797 < 0.01
assert np.abs(result[0]*100 - 6.70294094) / 6.702 < 0.01

cfes = np.zeros(nmc)
fbound = np.zeros(nmc)
fcce = np.zeros(nmc)



if os.path.exists("cfe_global_table.txt"):
    tbl = Table.read("cfe_global_table.txt", format='ascii.csv')
    fbound = tbl['fbound']

else:

    for ii in ProgressBar(range(nmc)):
        # we write the parameters to file for each call, since the fortran code
        # reads these parameters from file
        with open('parameters.in', 'w') as fh:
            fh.write(
                textwrap.dedent("""
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
            )+"\n"))

        rslt = cfeglobal.cfemod.f_cfe(surfg_arr[ii],
                                      qT_arr[ii],
                                      Omega_arr[ii])
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

    tbl = table.Table(data=[sigma_arr, surfg_arr, fbound],
                      names=['sigma','surfg','fbound'])

    tbl.write("cfe_global_table.txt", format='ascii.csv')


if __name__ == "__main__":

    from mpl_plot_templates import adaptive_param_plot
    import pylab as pl
    pl.clf()

    bins = np.array([np.logspace(1.5,4,15), np.logspace(0.5,2,15)])
    rslt = adaptive_param_plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value,
                               fbound, marker='none',
                               bins=bins,
                               percentilelevels=[0.05, 0.32],
                              )

    # attempt to plot covariance
    # this is to enable a (theoretical) sanity check because of a heisenbug
    # that has cropped up dozens of times in which the global data become
    # uncorrelated in the surfg-fbound axis, which is not correct.
    data = np.log10(np.array([(surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value, fbound]))
    cov = np.cov(data)
    var = np.diag(cov)
    rot = cov/var
    mn = np.mean(data, axis=1)

    vecs = [[-1,0],
            [ 1,0],
            [ 0,-1],
            [ 0, 1],]
    vecs = np.array(vecs)# * var**0.5
    rotvec = np.dot(cov/var**0.5, np.transpose(vecs))
    lines = rotvec.T + mn

    angle = (np.arctan2(cov[0,0], cov[0,1])*u.rad).to(u.deg)
    print("rotation angle of fit to data",angle)

    pl.plot(lines.T[0,:2], lines.T[1,:2])
    pl.plot(lines.T[0,2:], lines.T[1,2:])
