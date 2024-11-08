from hoomd_script import *
from sys import argv
import math
import numpy as np 

def pbc(r,d):
    return r - d * np.round(r / d)

def unwrap_pbc(r,d):
    ref_atom = r[0]
    unwrapped_traj = np.copy(r)
    n_atoms, _ = r.shape
    for j in range(n_atoms):
        delta = unwrapped_traj[j] - ref_atom
        delta -= np.round(delta / d) * d
        unwrapped_traj[j] = ref_atom + delta
    return unwrapped_traj

#def qiucm(pos,box):
#    cm = pbc(pos[...,0,:]+pbc(np.diff(pos, prepend=pos[...,:1,:],axis=-2),box).cumsum(axis=-2).mean(axis=-2),box)
#    return cm

def qiucm(pos,box):
    dpos = np.append(np.zeros(((1,3))),pbc(np.diff(pos,axis=0),box),axis = 0)
    rpos = np.cumsum(dpos,axis = 0)
    cm = rpos.mean(axis = 0)
    cm = cm +pos[0]
    cm = pbc(np.mean(pos,axis = 0),box)
    return cm
    
a = init.read_xml(argv[1])
pos = []
ima = []
ty = []
bond = []
body = []
mass = []
for i in a.bonds:
	bond.append([i.type,i.a,i.b])
Lx,Ly,Lz = a.box.Lx,a.box.Ly,a.box.Lz

for i in a.particles:
    pos.append(i.position)
    ima.append(i.image)
    ty.append(i.type)
    mass.append(i.mass)
    body.append(i.body)

pos = np.array(pos)
box = np.array([a.box.Lx, a.box.Ly, a.box.Lz])

qiu = 300
gb = group.type(type = 'B')
gc = group.type(type = 'C')
gd = group.type(type = 'D')
gbc = group.union(name='gbc',a=gb,b=gc)
gg = group.union(name='combined',a=gbc,b=gd)

qiupos = []
for i in range(qiu):
    qiupos.append(gg[i].position)
qiupos = np.array(qiupos)
qiupos = unwrap_pbc(qiupos,box) 
print(qiupos)
cm = qiucm(qiupos,box)
cm = np.array(cm)
print(cm)

new_pos = []
for i in range(len(pos)):
    new_pos.append(pbc(pos[i] - cm,box))

natoms = len(ty)
o = open('cm_dd_remake_SCNP_PAAM.xml', 'w')
o.write("<?xml version = '1.0' encoding = 'UTF-8'?>\n<galamost_xml version = '1.4'>\n<configuration time_step = '0' dimensions = '3' natoms = '%s'>\n<box lx = '%s' ly = '%s' lz = '%s'/>\n<position num = '%s'>\n" % (natoms, Lx, Ly, Lz, natoms))
o.write('\n'.join(["%s %s %s" % (str(p[0]), str(p[1]), str(p[2])) for p in new_pos]))
o.write("\n</position>\n")
o.write("<body num='%s'>\n" % (natoms))
o.write('\n'.join([ str(x) for x in body ]))
o.write("\n</body>\n")
o.write("<type num='%s'>\n" % (natoms))
o.write('\n'.join(ty))
o.write("\n</type>\n")
o.write("<image num='%s'>\n" % (natoms))
o.write('\n'.join([ "%s %s %s" % (str(p[0]), str(p[1]), str(p[2])) for p in ima ]))
o.write("\n</image>\n")
o.write("<bond num='%s'>\n" % (len(bond)))
o.write('\n'.join([ "%s %s %s" % (str(p[0]), str(p[1]), str(p[2])) for p in bond ]))
o.write("\n</bond>\n")
#o.write("<h_init num='%s'>\n" % (len(hi)))
#o.write('\n'.join(hi))
#o.write("\n</h_init>\n")
#o.write("<h_cris num='%s'>\n" % (len(hc)))
#o.write('\n'.join(hc))
#o.write("\n</h_cris>\n")
o.write("\n</configuration>\n</galamost_xml>")
o.close()

