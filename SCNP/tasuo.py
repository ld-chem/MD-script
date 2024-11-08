from hoomd_script import *
from sys import argv
fl = argv[1]
# read in the file
init.read_xml(filename=fl)
# force field setup
harmonic = bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=1.0)
harmonic.bond_coeff.set('cross', k=330.0, r0=1.0)

#fene = bond.fene()
#fene.bond_coeff.set('graft', k=30.0, r0=1.5, sigma=1.0, epsilon=1.0)
#fene.bond_coeff.set('polymer', k=30.0, r0=1.5, sigma=1.0, epsilon=1.0)
lj = pair.lj(r_cut=2**(1/6))
lj.set_params(mode="shift")
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0,r_cut = 3.0)
lj.pair_coeff.set('C', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('C', 'C', epsilon=1.0, sigma=1.0,r_cut = 3.0)
lj.pair_coeff.set('B', 'C', epsilon=1.0, sigma=1.0,r_cut = 3.0)
lj.pair_coeff.set('B', 'D', epsilon=1.0, sigma=1.0,r_cut = 3.0)
lj.pair_coeff.set('D', 'D', epsilon=1.0, sigma=1.0,r_cut = 3.0)
nlist.reset_exclusions(exclusions = ['bond','body'])


xml=dump.xml(filename='particles', period = 1e7)
xml.set_params( position=True, type=True, bond=True, body=True, image=True)
#analyze.log(filename='momentum.log',quantities=['momentum'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='temperature.log',quantities=['temperature_b-particles'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='kinetic_energy.log',quantities=['kinetic_energy_b-particles'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='potential_energy.log',quantities=['potential_energy_b-particles'],period=100, header_prefix='#',overwrite=True)

# integrate NVT for a bunch of time steps


# integrate at constant temperature,momentum
#update.zero_momentum(period=100)
integrate.mode_standard(dt=0.005)
integrator = integrate.langevin(group=group.nonrigid(), T=1.0, seed=5)
run (1e7)
run(1)




