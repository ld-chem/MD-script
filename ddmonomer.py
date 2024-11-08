from sys import argv
from hoomd_script import *
import numpy as np
##delete unbonded monomer##

a = init.read_xml(argv[1])

pos = []
ima = []
ty = []
bond = []
body = []

for i in a.bonds:
	bond.append([i.type,i.a,i.b])
Lx,Ly,Lz = a.box.Lx,a.box.Ly,a.box.Lz

bond = np.array(bond)
tb = bond.transpose((1,0))
tb_0 = tb[0]
tb_1 = tb[1]
tb_2 = tb[2]
tbbc = np.hstack((tb[1],tb[2]))
tbbcu = np.unique(tbbc)
nb = tbbcu.astype(int)

c = 0
new_i = 0
num = 118477 # sum of particles && modifiable
new = np.zeros(num, dtype=int)

for i in a.particles:
    if c in nb:
        pos.append(i.position)
        ima.append(i.image)
        ty.append(i.type)
        body.append(-1)
        new[c] = new_i
        new_i += 1
        c += 1
    else:
        c += 1

new = new.astype(int)
newtb_1 = []
newtb_2 = []

for i in tb_1:
    newtb_1.append(new[int(i)])

for i in tb_2:
    newtb_2.append(new[int(i)])

tb_0 = np.array(tb_0)
newtb_1 = np.array(newtb_1)
newtb_2 = np.array(newtb_2)
new_bond = np.vstack((tb_0,newtb_1,newtb_2))
bond = new_bond.transpose((1,0))
			
natoms = len(ty)
o = open('ddmn_noSCNPPAAM_treu.xml', 'w')
o.write("<?xml version = '1.0' encoding = 'UTF-8'?>\n<galamost_xml version = '1.4'>\n<configuration time_step = '0' dimensions = '3' natoms = '%s'>\n<box lx = '%s' ly = '%s' lz = '%s'/>\n<position num = '%s'>\n" % (natoms, Lx, Ly, Lz, natoms))
o.write('\n'.join(["%s %s %s" % (str(p[0]), str(p[1]), str(p[2])) for p in pos]))
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
