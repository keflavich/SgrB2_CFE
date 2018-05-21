import pylab as pl
import numpy as np
from astropy import units as u

# Adamo+ 2015 for most of the data
# Kruijssen & Bastian 2016 for gas surface densities
galaxy_data = {'SMC': {'SigSFR':0.001, 'Gamma':4.2, 'eGamma':(0.3,0.2), 'SigGas': 0.96, 'distance': 60*u.kpc},
               'LMC': {'SigSFR':1.52e-3, 'Gamma':5.8, 'eGamma':(0.5,0.5), 'SigGas': 1.04, 'distance': 50*u.kpc},
               'NGC 3256': {'SigSFR': 0.62, 'Gamma': 22.9, 'eGamma':(9.8, 7.3), 'SigGas':0, 'distance': 19.8*u.Mpc},
               'Solar Neighborhood': {'SigSFR': 0.012, 'Gamma': 7.0, 'eGamma': (3,7), 'SigGas':0, 'distance':500*u.pc},
               'NGC 45': {'SigSFR': 1.02e-3, 'Gamma': 5.2, 'eGamma':(0.3,0.3), 'SigGas':0},
               'NGC 1313': {'SigSFR':0.011, 'Gamma':3.2, 'eGamma':(0.2,0.2), 'SigGas': 0},
               'NGC 4395': {'SigSFR':4.66e-3, 'Gamma':1, 'eGamma':(0.6,0.6), 'SigGas': 0},
               'NGC 7793': {'SigSFR':6.51e-3, 'Gamma':2.5, 'eGamma':(0.3,0.3), 'SigGas': 0},
               'NGC 4449': {'SigSFR':0.04, 'Gamma':9.0, 'eGamma':(0,0), 'SigGas': 0},
               'NGC 1569': {'SigSFR':0.03, 'Gamma':13.9, 'eGamma':(0.8,0.8), 'SigGas': 0},
               #'Dwarf Sample': {'SigSFR':, 'Gamma':, 'eGamma':, 'SigGas': 0},
               'IC 10': {'SigSFR':0.03, 'Gamma':4.2, 'eGamma':(0,0), 'SigGas': 0},
               'ESO 338': {'SigSFR': 1.55, 'Gamma':50, 'eGamma':(10,10), 'SigGas': 0},
               'Haro 11': {'SigSFR': 2.16, 'Gamma':50, 'eGamma':(15,13), 'SigGas': 0},
               'ESO 185-IG13': {'SigSFR': 0.52, 'Gamma':26, 'eGamma':(5,5), 'SigGas': 0},
               'MRK 930': {'SigSFR':0.59, 'Gamma':25, 'eGamma':(10,10), 'SigGas': 0},
               'SBS 0335-052E': {'SigSFR':0.95, 'Gamma':49, 'eGamma':(15,15), 'SigGas': 0},
               'NGC 2997': {'SigSFR':9.4e-3, 'Gamma':10, 'eGamma':(2.6,2.6), 'SigGas': 0},
               'M83 center': {'SigSFR':0.54, 'Gamma':26.7, 'eGamma':(4,5.3), 'SigGas': 0, 'distance':4.85*u.Mpc},
               'M83 middle': {'SigSFR':0.013, 'Gamma':18.2, 'eGamma':(3,3), 'SigGas': 1.70, 'distance':4.85*u.Mpc},
               'M83 outer': {'SigSFR':0.013, 'Gamma':5.6, 'eGamma':(0.6,0.6), 'SigGas': 0, 'distance':4.85*u.Mpc},
               'NGC 6946': {'SigSFR':4.6e-3, 'Gamma':12.5, 'eGamma':(2.5,1.8), 'SigGas': 0},
              }


sigsfr = np.array([x['SigSFR'] for x in galaxy_data.values()])*u.M_sun/u.yr
gamma = np.array([x['Gamma'] for x in galaxy_data.values()])
egamma_low,egamma_high = np.array([x['eGamma'] for x in galaxy_data.values()]).T


pl.clf()
ax = pl.gca()
ax.errorbar(sigsfr.value, gamma, yerr=np.array([egamma_low, egamma_high]),
            marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax.set_xscale('log')
ax.set_yscale('log')
