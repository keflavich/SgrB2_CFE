import pylab as pl
import numpy as np
from astropy import units as u
from astropy.table import Table, Column

pl.matplotlib.rcParams['xtick.direction'] = 'in'
pl.matplotlib.rcParams['ytick.direction'] = 'in'
pl.rc('font', family='serif')
pl.matplotlib.rcParams['font.serif'] = 'Times New Roman'
pl.matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'
pl.matplotlib.rcParams['mathtext.rm'] = 'dejavuserif'
pl.rc('text', usetex=True)

# Adamo+ 2015 for most of the data
# Kruijssen & Bastian 2016 for gas surface densities
galaxy_data = {'SMC': {'SigSFR':0.001, 'Gamma':4.2, 'eGamma':(0.3,0.2), 'SigGas': 10**0.96, 'distance': 60*u.kpc},
               'LMC': {'SigSFR':1.52e-3, 'Gamma':5.8, 'eGamma':(0.5,0.5), 'SigGas': 10**1.04, 'distance': 50*u.kpc},
               'NGC 3256': {'SigSFR': 0.62, 'Gamma': 22.9, 'eGamma':(9.8, 7.3), 'SigGas': 10**2.1, 'distance': 19.8*u.Mpc},
               'Solar Neighborhood': {'SigSFR': 0.012, 'Gamma': 7.0, 'eGamma': (3,7), 'SigGas': 10**0.97, 'distance':500*u.pc},
               'NGC 45': {'SigSFR': 1.02e-3, 'Gamma': 5.2, 'eGamma':(0.3,0.3), 'SigGas': 10**0.67, 'distance': 6.8*u.Mpc},
               'NGC 1313': {'SigSFR':0.011, 'Gamma':3.2, 'eGamma':(0.2,0.2), 'SigGas': 10**0.88, 'distance': 4.15*u.Mpc},
               'NGC 4395': {'SigSFR':4.66e-3, 'Gamma':1, 'eGamma':(0.6,0.6), 'SigGas': 10**0.59, 'distance': 4.23*u.Mpc},
               'NGC 7793': {'SigSFR':6.51e-3, 'Gamma':2.5, 'eGamma':(0.3,0.3), 'SigGas': 10**0.84, 'distance': 3.775*u.Mpc},
               'NGC 4449': {'SigSFR':0.04, 'Gamma':9.0, 'eGamma':(0,0), 'SigGas': np.nan, 'distance': 3.98*u.Mpc},
               'NGC 1569': {'SigSFR':0.03, 'Gamma':13.9, 'eGamma':(0.8,0.8), 'SigGas': 10**1.33, 'distance': 2.63*u.Mpc},
               #'Dwarf Sample': {'SigSFR':, 'Gamma':, 'eGamma':, 'SigGas': 0},
               'IC 10': {'SigSFR':0.03, 'Gamma':4.2, 'eGamma':(0,0), 'SigGas': np.nan, 'distance': 0.792*u.Mpc},
               'ESO 338': {'SigSFR': 1.55, 'Gamma':50, 'eGamma':(10,10), 'SigGas': np.nan, 'distance': np.nan*u.Mpc},
               'Haro 11': {'SigSFR': 2.16, 'Gamma':50, 'eGamma':(15,13), 'SigGas': 10**1.98, 'distance': 92*u.Mpc}, # wikipedia
               'ESO 185-IG13': {'SigSFR': 0.52, 'Gamma':26, 'eGamma':(5,5), 'SigGas': np.nan, 'distance': 80*u.Mpc}, # Adamo+ 2011
               'MRK 930': {'SigSFR':0.59, 'Gamma':25, 'eGamma':(10,10), 'SigGas': np.nan, 'distance': 72*u.Mpc}, # Adamo+ 2011
               'SBS 0335-052E': {'SigSFR':0.95, 'Gamma':49, 'eGamma':(15,15), 'SigGas': np.nan, 'distance': 54.3*u.Mpc}, # Pustilnik+ 2000
               'NGC 2997': {'SigSFR':9.4e-3, 'Gamma':10, 'eGamma':(2.6,2.6), 'SigGas': np.nan, 'distance': 12.2*u.Mpc}, # Hess+ 2009
               'M83 center': {'SigSFR':0.54, 'Gamma':26.7, 'eGamma':(4,5.3), 'SigGas': np.nan, 'distance': 4.85*u.Mpc},
               'M83 middle': {'SigSFR':0.013, 'Gamma':18.2, 'eGamma':(3,3), 'SigGas': 10**1.70, 'distance': 4.85*u.Mpc},
               'M83 outer': {'SigSFR':0.013, 'Gamma':5.6, 'eGamma':(0.6,0.6), 'SigGas': np.nan, 'distance': 4.85*u.Mpc},
               'NGC 6946': {'SigSFR':4.6e-3, 'Gamma':12.5, 'eGamma':(2.5,1.8), 'SigGas': 10**1.30, 'distance': 5.9*u.Mpc}, # Karachentsev+ 2000
              }

m83surf = Table.read('m83_profiles.csv')
m83cfe = Table.read('m83_profiles_adamo.txt', format='ascii')

surfs = []
for row in m83cfe:
    sel = (m83surf['Radius'] < row['outer']) & (m83surf['Radius'] > row['inner'])
    avsurf = m83surf[sel]['Surfdens_H2'].mean()
    surfs.append(avsurf)

m83cfe.add_column(Column(name='surfdens_h2', data=surfs))




siggas = np.array([x['SigGas'] for x in galaxy_data.values()])*u.M_sun/u.pc**2
sigsfr = np.array([x['SigSFR'] for x in galaxy_data.values()])*u.M_sun/u.yr
gamma = np.array([x['Gamma'] for x in galaxy_data.values()])
egamma_low,egamma_high = np.array([x['eGamma'] for x in galaxy_data.values()]).T
distance = u.Quantity([x['distance'] for x in galaxy_data.values()])

sgrb2_data = Table.read('../tables/cluster_mass_estimates_cfe.csv')
cfe_sb2 = sgrb2_data['$M_{inferred}$'][-1] / sgrb2_data['$M_{inferred}$'][-2] * 100
cfe_sb2_low = 37
cfe_sb2_high = 43
cfe_sb2_yerr = np.array([[cfe_sb2-cfe_sb2_low, cfe_sb2_high-cfe_sb2]]).T
sgrb2_surfdens = 1e3 * u.Msun/u.pc**2
# this is not right...
sgrb2_sfr_surfdens = (0.062 * u.Msun/u.yr/(15*u.pc)**2).to(u.Msun/u.yr/u.kpc**2)
esgrb2_sfr_surfdens = (0.1 * 0.062 * u.Msun/u.yr/(15*u.pc)**2).to(u.Msun/u.yr/u.kpc**2)
sgrb2_distance = (8.5*u.kpc).to(u.Mpc)
stream_area = (120**2-60**2)*np.pi*u.pc**2
sgrb2_sfr_surfdens = (0.1 * u.Msun/u.yr/(stream_area)).to(u.M_sun/u.yr/u.kpc**2)
min_stream_area, max_stream_area = (100**2-80**2)*np.pi*u.pc**2, (140**2-60**2)*np.pi*u.pc**2
esgrb2_sfr_surfdens = np.array(((sgrb2_sfr_surfdens - (0.1 * u.Msun/u.yr/(max_stream_area)).to(u.M_sun/u.yr/u.kpc**2)).value,
                               ((0.1 * u.Msun/u.yr/(min_stream_area)).to(u.M_sun/u.yr/u.kpc**2) - sgrb2_sfr_surfdens).value))


fig = pl.figure(1)
fig.clf()
ax = fig.gca()
ax.errorbar(sigsfr.value, gamma, yerr=np.array([egamma_low, egamma_high]),
            marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax.errorbar(sgrb2_sfr_surfdens.value, cfe_sb2, xerr=np.array([esgrb2_sfr_surfdens]).T,
            yerr=cfe_sb2_yerr,
            marker='o', linestyle='none', markeredgecolor='r',
            markerfacecolor='k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("$\Sigma_{\mathrm{SFR}}$ [M$_\odot$ yr$^{-1}$ pc$^{-2}$]", fontsize=18)
ax.set_ylabel("$\Gamma$", fontsize=18)
ax.set_ylim(1,100)
ax.tick_params(labelsize=16)
fig.savefig('GammaVsSigmaSFR.pdf', bbox_inches='tight')

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()
ax2.errorbar(distance.to(u.Mpc).value, gamma, yerr=np.array([egamma_low, egamma_high]),
             marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax2.errorbar(sgrb2_distance.value, cfe_sb2,
             yerr=cfe_sb2_yerr,
             marker='o', linestyle='none', markeredgecolor='r',
             markerfacecolor='k')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("Distance [Mpc]", fontsize=18)
ax2.set_ylabel("$\Gamma$", fontsize=18)
ax2.tick_params(labelsize=16)
ax2.set_ylim(1,100)
fig.savefig('GammaVsSigmaSFR.pdf', bbox_inches='tight')
fig2.savefig('GammaVsDistance.pdf', bbox_inches='tight')

fig3 = pl.figure(3)
fig3.clf()
ax3 = fig3.gca()
ax3.errorbar(siggas.value, gamma, yerr=np.array([egamma_low, egamma_high]),
             marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax3.errorbar(sgrb2_surfdens.value, cfe_sb2,
             xerr=np.array([[sgrb2_surfdens.value-sgrb2_surfdens.value/1.65, sgrb2_surfdens.value*1.65-sgrb2_surfdens.value]]).T,
             yerr=cfe_sb2_yerr,
            marker='o', linestyle='none', markeredgecolor='r',
            markerfacecolor='orange')
ax3.errorbar(m83cfe['surfdens_h2'][m83cfe['zone'] == 'eqarea'],
             m83cfe['cfe'][m83cfe['zone'] == 'eqarea'],
             yerr=m83cfe['ecfe'][m83cfe['zone'] == 'eqarea'],
            )
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel("$\Sigma_{\mathrm{gas}}$ [M$_\odot$ pc$^{-2}$]", fontsize=18)
ax3.set_ylabel("$\Gamma$", fontsize=18)
ax3.tick_params(labelsize=16)
ax3.set_ylim(1,100)
ax3.set_xlim(0.5, 7e3)
fig3.savefig('GammaVsSigmaGas.pdf', bbox_inches='tight')



import imp
#import parameters
#import cfe_global_plots
#import cfe_local_plots
#imp.reload(parameters)
#imp.reload(cfe_global_plots)
#imp.reload(cfe_local_plots)
from cfe_global_plots import sigma_arr, surfg_arr, fbound as global_fbound
from cfe_local_plots import fbound as local_fbound

from matplotlib.colors import LinearSegmentedColormap

import sys
imp.reload(sys.modules['mpl_plot_templates.adaptive_param_plot'])
imp.reload(sys.modules['mpl_plot_templates'])
from mpl_plot_templates import adaptive_param_plot

bins = 10**np.mgrid[1.5:3.5:0.25, 0.5:2:0.1]
bins = np.array([np.logspace(1.5,4,15), np.logspace(0.5,2,15)])

cdict1 = {'red':   ((0.0, 1.0, 1.0),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
cdictblue = {'red':   ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 1.0, 1.0),
                   (1.0, 1.0, 1.0))
        }
cm_red = LinearSegmentedColormap('red', cdict1)
cm_blue = LinearSegmentedColormap('blue', cdictblue)

#ax3.plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value, cfes, 'k.', alpha=0.25, zorder=-5)
rslt = adaptive_param_plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value,
                           global_fbound, marker='none',
                           #levels=[1-0.95,1-0.68])
                            #colors=['b']*15,
                            cmap=cm_red,
                           bins=bins,
                           percentilelevels=[0.05, 0.32],
                          )
rslt2 = adaptive_param_plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value,
                            local_fbound, marker='none',
                            #levels=[1-0.95,1-0.68])
                            #colors=['r']*15,
                            cmap=cm_blue,
                            linestyles='dotted',
                           bins=bins,
                           percentilelevels=[0.05, 0.32],
                           )
ax3.set_ylim(1,100)
ax3.set_xlim(0.5, 7e3)

fig3.savefig('GammaVsSigmaGas_withTheory.pdf', bbox_inches='tight')


fbound_median = np.nanmedian(global_fbound)
fbound_errmin = fbound_median - np.nanpercentile(global_fbound, 16.)
fbound_errmax = np.nanpercentile(global_fbound, 84.) - fbound_median

print('fbound_global: %.2f+%.2f-%.2f' % (np.around(fbound_median,2),np.around(fbound_errmax,2),np.around(fbound_errmin,2)))

fbound_median = np.nanmedian(local_fbound)
fbound_errmin = fbound_median - np.nanpercentile(local_fbound, 16.)
fbound_errmax = np.nanpercentile(local_fbound, 84.) - fbound_median

print('fbound_local: %.2f+%.2f-%.2f' % (np.around(fbound_median,2),np.around(fbound_errmax,2),np.around(fbound_errmin,2)))





fig4 = pl.figure(4, figsize=(10,4))
fig4.clf()
ax = fig4.add_subplot(1,3,2)
ax.errorbar(sigsfr.value, gamma, yerr=np.array([egamma_low, egamma_high]),
            marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax.errorbar(sgrb2_sfr_surfdens.value, cfe_sb2, xerr=np.array([esgrb2_sfr_surfdens]).T,
            yerr=cfe_sb2_yerr,
            marker='o', linestyle='none', markeredgecolor='r',
            markerfacecolor='k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("$\Sigma_{\mathrm{SFR}}$ [M$_\odot$ yr$^{-1}$ pc$^{-2}$]", fontsize=12)
ax.set_ylabel("$\Gamma$", fontsize=12)
ax.set_ylim(1,100)
ax.set_xlim(6e-4,15)
ax.tick_params(labelsize=10)
[xx.set_visible(False) for xx in ax.get_yticklabels()]

ax2 = fig4.add_subplot(1,3,3)
ax2.errorbar(distance.to(u.Mpc).value, gamma, yerr=np.array([egamma_low, egamma_high]),
             marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax2.errorbar(sgrb2_distance.value, cfe_sb2,
             yerr=cfe_sb2_yerr,
             marker='o', linestyle='none', markeredgecolor='r',
             markerfacecolor='k')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("Distance [Mpc]", fontsize=12)
#ax2.set_ylabel("$\Gamma$", fontsize=18)
ax2.tick_params(labelsize=10)
ax2.set_ylim(1,100)
[xx.set_visible(False) for xx in ax2.get_yticklabels()]

ax3 = fig4.add_subplot(1,3,1)
ax3.errorbar(siggas.value, gamma, yerr=np.array([egamma_low, egamma_high]),
             marker='s', linestyle='none', markeredgecolor='k', alpha=0.6)
ax3.errorbar(sgrb2_surfdens.value, cfe_sb2,
             xerr=np.array([[sgrb2_surfdens.value-sgrb2_surfdens.value/1.65, sgrb2_surfdens.value*1.65-sgrb2_surfdens.value]]).T,
             yerr=cfe_sb2_yerr,
            marker='o', linestyle='none', markeredgecolor='r',
            markerfacecolor='orange')
ax3.errorbar(m83cfe['surfdens_h2'][m83cfe['zone'] == 'eqarea'],
             m83cfe['cfe'][m83cfe['zone'] == 'eqarea'],
             yerr=m83cfe['ecfe'][m83cfe['zone'] == 'eqarea'],
            )
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel("$\Sigma_{\mathrm{gas}}$ [M$_\odot$ pc$^{-2}$]", fontsize=12)
ax3.tick_params(labelsize=10)
ax3.set_ylim(1,100)
ax3.set_xlim(0.5, 7e3)
#[xx.set_visible(False) for xx in ax3.get_yticklabels()]


rslt = adaptive_param_plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value,
                           global_fbound, marker='none',
                           #levels=[1-0.95,1-0.68])
                            #colors=['b']*15,
                            cmap=cm_red,
                           bins=bins,
                           percentilelevels=[0.05, 0.32],
                          )
rslt2 = adaptive_param_plot((surfg_arr*u.kg/u.m**2).to(u.M_sun/u.pc**2).value,
                            local_fbound, marker='none',
                            #levels=[1-0.95,1-0.68])
                            #colors=['r']*15,
                            cmap=cm_blue,
                            linestyles='dotted',
                           bins=bins,
                           percentilelevels=[0.05, 0.32],
                           )
ax3.set_ylim(1,100)
ax3.set_xlim(0.5, 7e3)
ax3.set_ylabel("$\Gamma$", fontsize=12)

fig4.subplots_adjust(wspace=0)

fig4.savefig('GammaVsEverything_3panel.pdf', bbox_inches='tight')
