import numpy as np
from sys import argv
from hoomd_script import *

def pbc(r,d):
	return r-d * np.round(r/d)

a = init.read_xml(argv[1])
box = np.array([a.box.Lx,a.box.Ly,a.box.Lz])
#gb = group.type(type = 'B')
#pos = []
#for i in gb:
#	pos.append(i.position)
#pos = np.array(pos)



lj = pair.lj(r_cut=2 ** (1/6))
wca = 2**(1.0/6.0)
lj.set_params(mode='shift')
lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0,r_cut = 3.0);
lj.pair_coeff.set('C', 'C', epsilon=1.0, sigma=1.0,r_cut = 3.0);
lj.pair_coeff.set('D', 'D', epsilon=1.0, sigma=1.0,r_cut = 3.0);
lj.pair_coeff.set('B', 'C', epsilon=1.0, sigma=1.0,r_cut = 3.0);
lj.pair_coeff.set('B', 'D', epsilon=1.0, sigma=1.0,r_cut = 3.0);
lj.pair_coeff.set('C', 'D', epsilon=1.0, sigma=1.0);
nlist.reset_exclusions(exclusions = ['bond','body'])

harmonic = bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0, r0=1.0)
harmonic.bond_coeff.set('cross', k=330.0, r0=1.0)
harmonic.bond_coeff.set('graft', k=330.0, r0=1.0)
#fen = bond.fene()
#fen.bond_coeff.set('polymer',k=30,r0=1.5,epsilon = 1.0,sigma = 1.0)
#fen.bond_coeff.set('cross',k=30,r0=1.5,epsilon = 1.0,sigma = 1.0)
# angle force
def cosine(theta, kappa, esp):
	V = kappa * esp * (1-np.cos(theta)) 
	T = kappa * esp * np.sin(theta)
	return (V, T)

#btable = angle.table(width=1000)
#btable.angle_coeff.set('an', func=cosine, coeff=dict(kappa=8, esp=1.0))

xml=dump.xml(filename='cross', period = 1e5)
xml.set_params( position=True, angle = True,type=True, body = True,bond=True, image=True)
integrate.mode_standard(dt=0.005)
update.zero_momentum(period=50000)
integrator = integrate.langevin(group= group.nonrigid(), T=1.0, seed=5)
#log = analyze.log(filename="analyze.log", period=1000, overwrite=False, quantities=['temperature', 'pressure', 'potential_energy'])

#run(1e7)

# ====================================================================
# begin crosslink


#r_balance =0.96`/n` 
crnum = 100 # the max value of cross bond
gd = group.type(type='D')
gb = group.type(type='B')
gc = group.type(type='C')

#nato = len(ga)
box = np.array([a.box.Lx,a.box.Ly,a.box.Lz])
def crosslink(timestep):
	used= []
	cr = 0
	ab = []
	for i in a.bonds:
		ab.append([i.type,i.a,i.b])
		if i.type == 'cross':
			if i.a ==0 and i.b ==1:
				pass
			else:
				cr+=1
				used.append(i.a)
				used.append(i.b)	
	pos = []
	di = []
	dii = []
	for i in range(len(gc)):
		p1 = np.array(gc[i].position)
		t1 = gc[i].tag
		#di = []
		if t1 not in used:
			for j in range(len(gb)):
				p2 = np.array(gb[j].position)
				t2 = gb[j].tag
				if t2 not in used:
					dis = np.sqrt(np.sum(pbc(p1-p2,box)**2))
					if dis < 1.2 and abs(t2-t1)>=4:
						di.append([dis,t1,t2])
	di.sort()
	for i in range(len(di)):
		if di[i][1] not in used and di[i][2] not in used:
			a.bonds.add('cross',di[i][1],di[i][2])
			used.append(di[i][1])
			used.append(di[i][2])
	
	for i in range(len(gd)):
				p3 = np.array(gd[i].position)
				t3 = gd[i].tag
				#di = []
				if t3 not in used:
						for j in range(len(gb)):
								p4 = np.array(gb[j].position)
								t4 = gb[j].tag
								if t4 not in used:
										dis = np.sqrt(np.sum(pbc(p3-p4,box)**2))
										if dis < 1.2 and abs(t4-t3)>=4:
												dii.append([dis,t3,t4])
	dii.sort()
	for i in range(len(dii)):
		if dii[i][1] not in used and dii[i][2] not in used:
			a.bonds.add('cross',dii[i][1],dii[i][2])
			used.append(dii[i][1])
			used.append(dii[i][2])
analyze.callback(callback=crosslink, period=1e4)
run(1e5)
run(1)
