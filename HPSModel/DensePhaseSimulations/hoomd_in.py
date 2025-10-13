import hoomd
hoomd.context.initialize()
from hoomd import md
from hoomd import azplugins
import numpy as np

## Simulation informaiton
T = 300.0 #Kelvin
kB = 0.0019872067 #Kcal/mol/Kelvin
kT = kB*T #Kcal/mol
dt = 0.2045814 #10fs
run_eq = 4000000 #40 nanoseconds
veryshort_period = 10
short_period = 500 
long_period = 50000 
veryshortrun_nve = 15000 
shortrun_nve = 750000 
longrun_nve = 100000000 #1 microsecond

## Initialization
system = hoomd.init.read_gsd('prod.gsd')
sys_snapshot = system.take_snapshot(all = True)
lbox = sys_snapshot.box.Lx #Assuming the simulation box is cubic

## Interactions
### Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])

### Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('1', k=20.000000, r0=3.800000)

### Pairwise interactions
nb = azplugins.pair.ashbaugh(r_cut = 0.0, nlist = nl)
nb.pair_coeff.set('1', '1', epsilon = 0.200000, sigma = 5.920000, lam = 0.459459, r_cut = 23.680000)
nb.pair_coeff.set('1', '2', epsilon = 0.200000, sigma = 6.140000, lam = 0.486486, r_cut = 24.560000)
nb.pair_coeff.set('2', '2', epsilon = 0.200000, sigma = 6.360000, lam = 0.513514, r_cut = 25.440000)

### Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut = 0.0, nlist = nl)
yukawa.pair_coeff.set(['1', '2'], ['1', '2'], epsilon=0.000000, kappa=0.000000, r_cut=0.000000)
yukawa.pair_coeff.set('1','1', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('1','2', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','1', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','2', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)

## Integration
#Particles MD Integrator
solute = hoomd.group.all()

#Zeroing the linear momentum to avoid any possible COM drift
fix_mum=hoomd.md.update.zero_momentum(period=1)

## Running 2 microsecond long LD simulation
hoomd.md.integrate.mode_standard(dt=dt) #time step is 10fs
integrator2 = hoomd.md.integrate.langevin(group = solute, kT = kT, seed = np.random.randint(1,1000000)) #Temp is kT/0.0019872067
integrator2.set_gamma('1',gamma=0.006312)
integrator2.set_gamma('2',gamma=0.006268)

## writing outputs and run the simulation
eqthermo = hoomd.analyze.log(filename='afteq_run.log', quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
hoomd.run(run_eq)
eqthermo.disable()

#Very short, short, and long simulation runs
longthermo = hoomd.analyze.log("long_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = long_period*5, overwrite = True, header_prefix = '#')
longpress = hoomd.analyze.log(filename='long_pressure.log', quantities=['pressure_xx', 'pressure_yy', 'pressure_zz', 'pressure_xy', 'pressure_xz', 'pressure_yz'], period=500, overwrite=True, header_prefix='#')
longgsd_file = hoomd.dump.gsd("long_dumps.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
writegsd = hoomd.dump.gsd("long_restart.gsd", period = long_period, dynamic = ['property','momentum'], group = solute, overwrite=True, truncate=True)
hoomd.run(tsteps=longrun_nve)
longthermo.disable()
longpress.disable()
longgsd_file.disable()

shortthermo = hoomd.analyze.log("short_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = short_period*5, overwrite = True, header_prefix = '#')
shortgsd_file = hoomd.dump.gsd("short_dumps.gsd", period = short_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
hoomd.run(tsteps=shortrun_nve)
shortthermo.disable()
shortgsd_file.disable()

veryshortthermo = hoomd.analyze.log("veryshort_run.log", quantities = ['temperature','pressure','potential_energy','kinetic_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','lx','ly','lz'], period = veryshort_period*10, overwrite = True, header_prefix = '#')
veryshortgsd_file = hoomd.dump.gsd("veryshort_dumps.gsd", period = veryshort_period, dynamic = ['property','momentum'], group = solute, overwrite=True)
hoomd.run(tsteps=veryshortrun_nve)
veryshortthermo.disable()
veryshortgsd_file.disable()
