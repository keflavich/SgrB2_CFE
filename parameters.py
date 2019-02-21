import numpy as np
from astropy import units as u
from astropy import constants


nmc = 100000

qvir_arr = 10**np.random.normal(np.log10(1.1),0.4/1.1/np.log(10),nmc)
tview_arr = np.random.normal(0.74,0.16,nmc)
surfGMC_arr = 10.**np.random.normal(np.log10(4.1e3),1.7e3/4.1e3/np.log(10.),nmc)
beta0_arr = 10.**np.random.normal(np.log10(0.34),0.25/0.34/np.log(10.),nmc)
surfg_arr = (10.**np.random.normal(np.log10(1.e3),0.5e3/1.e3/np.log(10.),nmc)*u.Msun/u.pc**2).to(u.kg/u.m**2).value
sigma_arr = 10.**np.random.normal(np.log10(10.e3),1.5e3/10.e3/np.log(10.),nmc)
Omega_arr = (np.random.normal(1.8,.25,nmc) * 1/u.Myr).to(1/u.s).value
qT_arr = sigma_arr*Omega_arr*np.sqrt(2.*1.7)/np.pi/constants.si.G.value/surfg_arr
qT_arr_ = (sigma_arr * u.m/u.s) * (Omega_arr/u.s) * np.sqrt(2*1.7) / (np.pi * constants.G * surfg_arr*u.kg/u.m**2)


mean_mass_density = (1e4 * 2.8 * u.Da / u.cm**3).to(u.kg/u.m**3)
e_mass_density = 0.5 * mean_mass_density
mass_density_arr = 10.**np.random.normal(np.log10(mean_mass_density.value),
                                         e_mass_density.value/mean_mass_density.value/np.log(10),
                                         nmc)
mean_cs = (0.53*u.km/u.s).to(u.m/u.s)
e_cs = (0.07*u.km/u.s).to(u.m/u.s)

cs_arr = np.random.normal((mean_cs.value), (e_cs.value), nmc)


# check to see what Mach numbers we're getting
phi_pbar = 10 - 8*(1+0.025*(surfg_arr/(100*u.Msun.si.scale/u.pc.si.scale**2))**-2)**-1
#Q_arr = 2**0.5 * Omega_arr * sigma_arr / (np.pi * cfe.G * surfg_arr)
Mach_number_global = 2.82 * phi_pbar**(1/8.) * qT_arr * (Omega_arr/u.Myr.si.scale**-1)**-1 * (surfg_arr/(100*u.Msun.si.scale/u.pc.si.scale**2))

phi_p = 3
surfg_equilibrium = (2 * mass_density_arr * sigma_arr**2 / (np.pi * constants.G.si.value * phi_p))**0.5
