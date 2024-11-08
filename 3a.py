#!/usr/bin/python
####script head####
from sys import argv
from poetry import cu_gala as gala
from poetry import _options
filename = argv[1]

build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.004
app = gala.Application(all_info,dt)
neighbor_list = gala.NeighborList(all_info, 2*2**(1.0/6.0), 0.1)
neighbor_list.addExclusionsFromBonds()
lj = gala.LJForce(all_info, neighbor_list, 2*2**(1.0/6.0)) #( , , rcut)
lj.setParams('B', 'B', 2.0, 1.0, 1.0)
lj.setParams('B', 'C', 0.8, 1.0, 1.0)
lj.setParams('B', 'D', 0.5, 1.0, 1.0)
lj.setParams('C', 'D', 0.5, 1.0, 1.0)
lj.setParams('C', 'C', 1.2, 1.0, 1.0)
lj.setParams('D', 'D', 1.0, 1.0, 1.0)
lj.setParams('A', 'A', 1.0, 1.0, 1.0)
lj.setParams('A', 'B', 0.5, 1.0, 1.0)
lj.setParams('A', 'C', 0.5, 1.0, 1.0)
lj.setParams('A', 'D', 1.0, 1.0, 1.0)
lj.setParams('A', 'S', 1.0, 1.0, 1.0)
lj.setParams('B', 'S', 0.5, 1.0, 1.0)
lj.setParams('C', 'S', 0.8, 1.0, 1.0)
lj.setParams('D', 'S', 1.0, 1.0, 1.0)
lj.setParams('S', 'S', 1.0, 1.0, 1.0)
lj.setEnergy_shift()
app.add(lj)

bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('A-A', 330.0, 1.0)
bondforce.setParams('polymer', 330.0, 1.0) #(,K0, R0)
bondforce.setParams('cross', 330.0, 1.0)
app.add(bondforce)

group = gala.ParticleSet(all_info,'all')
comp_info = gala.ComputeInfo(all_info,group)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(10000)

zm = gala.ZeroMomentum(all_info)
zm.setPeriod(10000)

Temperature = 298
T = Temperature/120.27236  #reduced unit T
#P = 0.0610
P = 1

v = gala.VariantLinear()
v.setPoint(0, 80) # timesteps, temperature length
v.setPoint(20000000, 400)
v1 = gala.VariantRsqrt()
v1.setPoint(0, 80) # timesteps, temperature length
v1.setPoint(20000000, 400)
v1.setFactor(80)

axs = gala.AxialStretching(all_info, group)
axs.setBoxLength(v, 'Z')
axs.setBoxLength(v1, 'X')
axs.setBoxLength(v1, 'Y')
#axs.setRigidBody(True)
axs.setPeriod(1000)
app.add(axs)

#NH = gala.NoseHooverNVT(all_info, group, comp_info, T, 0.228)#( ,temperature, tau)
#app.add(NH)
NPT = gala.NPTMTK(all_info, group, comp_info, comp_info, T, P, 0.5, 5.0)# temperature,pressure,tauT,tauP
NPT.setAnisotropic(P, P, P)
#NPT.setSemiisotropic(P, P)
app.add(NPT)

dinfo = gala.DumpInfo(all_info, comp_info, 'data.log')
dinfo.setPeriod(5000)

##set stretching paramater

lx=all_info.getBasicInfo().getGlobalBox().getL().x
ly=all_info.getBasicInfo().getGlobalBox().getL().y
lz=all_info.getBasicInfo().getGlobalBox().getL().z

t = 0
step = int(1e5)

restart = gala.XMLDump(all_info,'restart')
restart.setPeriod(int(1e3))
restart.setOutputPosition(True)
restart.setOutputMass(True)
restart.setOutputType(True)
restart.setOutputImage(True)
restart.setOutputBond(True)
restart.setOutputAngle(True)
restart.setOutputDihedral(True)
#app.add(restart)

for i in range(1):
## load
    dinfo = gala.DumpInfo(all_info, comp_info, '%dth_data.log' % int(i+1))
    dinfo.setPeriod(100)
    dinfo.dumpPressTensor()
    dinfo.dumpBoxSize()
    load = gala.Application(all_info, dt)
    v = gala.VariantLinear()
    v.setPoint(t, lz)#time step, box length.
    v.setPoint(t+step, lz*13)
    axs = gala.AxialStretching(all_info,group)
    axs.setBoxLength(v,'Z')
    axs.setPeriod(100)
    load.add(lj)
    load.add(bondforce)
#    load.add(angle)
#    load.add(dihedral)
    load.add(sort_method)
    load.add(zm)
    load.add(dinfo)
    load.add(axs)
    load.add(restart)
    load.add(NPT)
    load.run(step)
    t = t + step
## unload
    unload = gala.Application(all_info, dt)
    v = gala.VariantLinear()
    v.setPoint(t, lz*13)#time step, box length.
    v.setPoint(t+step, lz)
    axs = gala.AxialStretching(all_info,group)
    axs.setBoxLength(v,'Z')
    axs.setPeriod(100)
    unload.add(lj)
    unload.add(bondforce)
#    unload.add(angle)
#    unload.add(dihedral)
    unload.add(sort_method)
    unload.add(zm)
    unload.add(dinfo)
    unload.add(axs)
    unload.add(NPT)
    unload.add(restart)
    unload.run(step)
    t = t + step
