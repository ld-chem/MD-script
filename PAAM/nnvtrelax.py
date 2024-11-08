from hoomd_script import *
from sys import argv
fl = argv[1]
# read in the file
init.read_xml(filename=fl)
# force field setup
harmonic = bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=1.0)
harmonic.bond_coeff.set('cross', k=330.0, r0=1.0)
harmonic.bond_coeff.set('A-A', k=330.0, r0=1.0)
#harmonic.bond_coeff.set('A-B', k=1111.0, r0=1.0)
#fene = bond.fene()
#fene.bond_coeff.set('cross', k=30.0, r0=1.5, sigma=1.0, epsilon=1.0)
#fene.bond_coeff.set('polymer', k=30.0, r0=1.5, sigma=1.0, epsilon=1.0)
lj = pair.lj(r_cut=2.5)
lj.set_params(mode="shift")
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'C', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('C', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('C', 'C', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('D', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'C', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('C', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('D', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('S', 'S', epsilon=1.0, sigma=1.0)
nlist.reset_exclusions(exclusions = ['bond','body'])

xml=dump.xml(filename='nvt', period = 5e6)
xml.set_params(position=True, type=True, body = True,bond=True, image=True)
analyze.log(filename='log.log',quantities=['temperature','pressure','lx'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='temperature.log',quantities=['temperature_b-particles'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='kinetic_energy.log',quantities=['kinetic_energy_b-particles'],period=100, header_prefix='#',overwrite=True)
#analyze.log(filename='potential_energy.log',quantities=['potential_energy_b-particles'],period=100, header_prefix='#',overwrite=True)

# integrate NVT for a bunch of time steps
typeA = group.type(name='a-particles', type='A')
typeB = group.type(name='b-particles', type='B')
typeC = group.type(name='c-particles', type='C')
typeD = group.type(name='d-particles', type='D')
#groupBCD = group.union(name="bcd-particles", a=typeB, b=typeB, c=typeC)
#groupBC = group.union(name="bc-particles", a=typeB, b=typeC)

integrate.mode_standard(dt=0.005)
#nve = integrate.nve(group=group.all(),limit=0.01)
#run(2e7)
#nve.disable()

# integrate at constant temperature,momentum
#update.zero_momentum(period=100)

relax = integrate.nvt(group=group.nonrigid(),T = 1, tau= 0.65)
run(2e9+1)

