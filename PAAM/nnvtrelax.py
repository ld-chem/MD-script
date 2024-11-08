from hoomd_script import *
from sys import argv

fl = argv[1]
init.read_xml(filename=fl)

harmonic = bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=1.0)
harmonic.bond_coeff.set('cross', k=330.0, r0=1.0)
harmonic.bond_coeff.set('A-A', k=330.0, r0=1.0)

lj = pair.lj(r_cut=2*2**(1.0/6.0))
lj.set_params(mode="shift")
lj.pair_coeff.set('B', 'B', epsilon=2.0, sigma=1.0)
lj.pair_coeff.set('B', 'C', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('B', 'D', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('C', 'D', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('C', 'C', epsilon=1.2, sigma=1.0)
lj.pair_coeff.set('D', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'B', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('A', 'C', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('A', 'D', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('A', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('B', 'S', epsilon=0.8, sigma=1.0)
lj.pair_coeff.set('C', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('D', 'S', epsilon=1.0, sigma=1.0)
lj.pair_coeff.set('S', 'S', epsilon=1.0, sigma=1.0)
nlist.reset_exclusions(exclusions = ['bond','body'])

xml=dump.xml(filename='nvt', period = 5e6)
xml.set_params(position=True, type=True, body = True,bond=True, image=True)
analyze.log(filename='log.log',quantities=['temperature','pressure','lx'],period=100, header_prefix='#',overwrite=True)

# integrate NVT for a bunch of time steps
typeA = group.type(name='a-particles', type='A')
typeB = group.type(name='b-particles', type='B')
typeC = group.type(name='c-particles', type='C')
typeD = group.type(name='d-particles', type='D')


integrate.mode_standard(dt=0.005)
nve = integrate.nve(group=group.all(),limit=0.01)
run(2e5)
nve.disable()

# integrate at constant temperature,momentum
update.zero_momentum(period=100)

relax = integrate.nvt(group=group.nonrigid(),T = 1, tau= 0.65)
run(2e8+1)
