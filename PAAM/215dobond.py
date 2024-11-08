#!/usr/bin/python
####script head#####
from sys import argv
from poetry import cu_gala as gala
from poetry import _options
filename = argv[1]

build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.001
app = gala.Application(all_info,dt)
neighbor_list = gala.NeighborList(all_info, 2*2**(1.0/6.0), 0.1) #(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds() # 1-2 exclusion
#neighbor_list.addExclusionsFromAngles() # 1-3 exclusion
#neighbor_list.addExclusionsFromDihedrals() #1-4 exclusion
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
T = 1.0
bd=gala.BerendsenNVT(all_info,group,comp_info,T,0.5)
app.add(bd)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(100)
app.add(sort_method)

dInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
dInfo.setPeriod(5000) # (period)
app.add(dInfo)

xml = gala.XMLDump(all_info, 'PAAM')
xml.setPeriod(int(5e6)) # (period)
xml.setOutputPosition(True)
xml.setOutputImage(True)
xml.setOutputMass(True)
xml.setOutputType(True)
xml.setOutputBody(True)
xml.setOutputBond(True)
#xml.setOutputAngle(True)
#xml.setOutputDihedral(True)
#xml.setOutputVelocity(True)
#xml.setOutputDiameter(True)
xml.setOutputInit(True)
xml.setOutputCris(True)
app.add(xml)

app.run(int(5e7))

reaction =gala.Polymerization(all_info,neighbor_list,2.0,16361)
reaction.setPr('A', 'A', 0.99)
#reaction.setNewBondTypeByPairs()
reaction.setPeriod(2000)
reaction.setInitInitReaction(True)
app.add(reaction)

app.run(int(2e9+1))
neighbor_list.printStats()
