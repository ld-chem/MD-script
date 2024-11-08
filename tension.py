from hoomd_script import *
import os
import numpy as np
from sys import argv

context.initialize()
filename = argv[1]
s = init.read_xml(filename, time_step=0)

Lx = L = s.box.Lx
Ly = s.box.Ly
Lz = s.box.Lz
sL = L * (0.85/1.05)**(1/3)

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

import math

harmonic = bond.harmonic(name="mybond")
harmonic.bond_coeff.set('A-A', k=330.0, r0=1.0)
harmonic.bond_coeff.set('polymer', k=330.0, r0=1.0)
harmonic.bond_coeff.set('cross', k=330.0, r0=1.0)

integrate.mode_standard(dt=0.001)
nvt = integrate.nvt(group=group.nonrigid(), T=0.2, tau=0.6)

sorter = update.sort()
sorter.set_period(100)
zeroer= update.zero_momentum(period=10)

eq = int(1e6)
pre0 = analyze.log(filename='pre0.log', quantities=['lx', 'ly', 'lz', 'pressure','pressure_xx','pressure_yy','pressure_zz','pressure_xy', 'pressure_xz', 'pressure_yz'], header_prefix='#', period=10000, overwrite=True)

ay = [ (i, (1e-6 * float((i-eq)) + 1) ** (-0.5) * Ly) for i in range(eq, eq+ int(1e6+1), 1)]
az = [ (i, (1e-6 * float((i-eq)) + 1) ** (-0.5) * Lz) for i in range(eq, eq+ int(1e6+1), 1)]

eq = 0
ay = [ (i, (1e-6 * float((i-eq)) + 1) ** (-0.5) * Ly) for i in range(eq, eq+ int(1e6+1), 1)]
az = [ (i, (1e-6 * float((i-eq)) + 1) ** (-0.5) * Lz) for i in range(eq, eq+ int(1e6+1), 1)]
box_resize = update.box_resize(Lx = variant.linear_interp([(0, Lx), (eq, Lx), (eq+int(1e6), 2*Lx)]), Ly = variant.linear_interp(ay), Lz= variant.linear_interp(az), period = 1)
box_resize.set_params(scale_particles=True)

xml = dump.xml(filename="tension", period=2e3)
xml.set_params(position=True, mass=True, diameter=True,type=True, image=True, bond=True, angle=True, velocity=True)

pre = analyze.log(filename='pre.log', quantities=['volume', 'lx', 'ly', 'lz', 'pressure','pressure_xx','pressure_yy','pressure_zz','pressure_xy', 'pressure_xz', 'pressure_yz'], header_prefix='#', period=100, overwrite=True)

a = s.take_snapshot()
a_idx = a.particles.typeid == a.particles.types.index('A')
b_idx = a.particles.typeid == a.particles.types.index('B')
c_idx = a.particles.typeid == a.particles.types.index('C')
d_idx = a.particles.typeid == a.particles.types.index('D')
chain_idx = a_idx + b_idx + c_idx + d_idx
box = np.array([a.box.Lx, a.box.Ly, a.box.Lz], dtype=np.float64)

idx_a = np.arange(a.particles.N)[chain_idx]

def break_bond(ts):
    a = s.take_snapshot()
    bond_list = s.bonds
    box = np.array([a.box.Lx,a.box.Ly,a.box.Lz])

    for bondi in bond_list:
        if bondi.a in idx_a or bondi.b in idx_a:
            posa = a.particles.position[bondi.a]+a.particles.image[bondi.a]*box
            posb = a.particles.position[bondi.b]+a.particles.image[bondi.b]*box
            bond_d = np.sqrt(np.sum((posa-posb)**2,axis=-1))
            if bond_d > 1.8:
                s.bonds.remove(bondi.tag)
                print(bond_d)
    print(len(s.bonds))

analyze = analyze.callback(callback=break_bond,period=5e4)
run(int(1e6+1))
