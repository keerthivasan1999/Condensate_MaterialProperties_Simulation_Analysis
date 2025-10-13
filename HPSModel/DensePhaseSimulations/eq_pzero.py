#importing the libraries
import sys, os, numpy as np
import hoomd, hoomd.md as md
from hoomd import azplugins
#import gsd, gsd.hoomd, gsd.pygsd
nsteps = 500000000 #simulation runtime length

temp=300 #simulation temperature

hoomd.context.initialize() #Initialize the execution context
system = hoomd.init.read_gsd("400.0x400.0x400.0_box.gsd") #Read initial system state from an GSD file.

nl = hoomd.md.nlist.cell() #Cell list based neighbor list
#### BOND DATA ####
harmonic = hoomd.md.bond.harmonic() #Harmonic bond potential and its parameters
harmonic.bond_coeff.set('1', k=20.000000, r0=3.800000)

nl.reset_exclusions(exclusions=['1-2', 'body']) #setting the exclusions from short range pair interactions
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl) #Ashbaugh-Hatch potential and its parameters
nb.pair_coeff.set('1', '1', epsilon = 0.200000, sigma = 5.920000, lam = 0.459459, r_cut = 23.680000)
nb.pair_coeff.set('1', '2', epsilon = 0.200000, sigma = 6.140000, lam = 0.486486, r_cut = 24.560000)
nb.pair_coeff.set('2', '2', epsilon = 0.200000, sigma = 6.360000, lam = 0.513514, r_cut = 25.440000)

### ELECTROSTATICS ###

yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
yukawa.pair_coeff.set(['1', '2'], ['1', '2'], epsilon=0.000000, kappa=0.000000, r_cut=0.000000)
yukawa.pair_coeff.set('1','1', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('1','2', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','1', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','2', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)

### MAKE PARTICLE GROUPS ###
solute = hoomd.group.all()


### NPT Integration ###
integrator1 = hoomd.md.integrate.npt(group=solute, kT=temp*0.0019872067,tau=100*0.2045814925,P=1.0*(1/68568.96063),tauP=1000*0.2045814925)

## Zeroing the linear momentum to avoid any possible COM drift
fix_mum=hoomd.md.update.zero_momentum(period=1)

## Running at dt = 1fs and P = 1 atm for short time
hoomd.md.integrate.mode_standard(dt = 0.2045814925/10)
preeqthermo = hoomd.analyze.log(filename='preeq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = 100000, overwrite = True, header_prefix = '#')
hoomd.run(tsteps=2500000)
preeqthermo.disable()

## Switching to dt = 10fs and linearly decreasing pressure from 1 atm to 0 atm for short time
hoomd.md.integrate.mode_standard(dt = 0.2045814925)
integrator1.set_params(P = hoomd.variant.linear_interp([(0,1.0*(1/68568.96063)),(2500000,0.0*(1/68568.96063))], zero='now'))
eqthermo = hoomd.analyze.log(filename='eq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = 100000, overwrite = True, header_prefix = '#')
hoomd.run(tsteps=2500000)
eqthermo.disable()

## Running at P = 0 atm for short time
hoomd.md.integrate.mode_standard(dt = 0.2045814925)
integrator1.set_params(P = 0.0*(1/68568.96063))
longthermo = hoomd.analyze.log('posteq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = 100000, overwrite = True, header_prefix = '#')
hoomd.run(tsteps=15000000)
longthermo.disable()
integrator1.disable()

## Running 1 microsecond long LD simulation
hoomd.md.integrate.mode_standard(dt=0.2045814) #time step is 10fs
integrator2 = hoomd.md.integrate.langevin(group=solute, kT=temp*0.0019872067,seed=np.random.randint(1,1000000)) #Langevin integrator
#setting the friction factor for the Langevin integrator [mass/1000fs]
integrator2.set_gamma('1',gamma=0.006312)
integrator2.set_gamma('2',gamma=0.006268)

box_data = np.loadtxt('posteq_run.log', usecols = [10], skiprows = 40, comments='#')
box_avg = np.mean(box_data)
lxf, lyf, lzf = box_avg, box_avg, box_avg
print(lxf, lyf, lzf)
print(system.box.Lx,system.box.Ly,system.box.Lz)

#Updating box size
hoomd.update.box_resize(Lx=hoomd.variant.linear_interp([(0,system.box.Lx),(1000000-100000,lxf)],  zero='now'),
Ly=hoomd.variant.linear_interp([(0,system.box.Ly),(1000000-100000,lyf)], zero='now'),
Lz=hoomd.variant.linear_interp([(0,system.box.Lz),(1000000-100000,lzf)],  zero='now'), scale_particles=True)
boxthermo = hoomd.analyze.log(filename='box_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = 10000, overwrite = True, header_prefix = '#')
gsdfile_forprod = hoomd.dump.gsd('prod.gsd', period = 100000, dynamic = ['property','momentum'], group = solute, overwrite=True, truncate=True)
hoomd.run(1000000)
boxthermo.disable()
