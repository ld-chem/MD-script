from sys import argv
from hoomd_script import *

a = init.read_xml(argv[1])
hi = []
hc = []
pos = []
ima = []
ty = []
bond = []
body = []
for i in a.bonds:
        bond.append([i.type,i.a,i.b])
Lx,Ly,Lz = a.box.Lx,a.box.Ly,a.box.Lz
dian = 59
chu = 0
for i in a.particles:
        pos.append(i.position)
        ima.append(i.image)
        #ty.append(i.type)
        if i.type == 'S' or i.type == 'B' or i.type =='C' or i.type =='D':
                hi.append('0')
                hc.append('0')
                body.append(-1)
                ty.append(i.type)
        else:
                if chu < dian:
                        hi.append('1')#initial concentration
                        hc.append('0')
                        body.append(-1)
                        ty.append('A')
                        chu+=1
                else:
                        hi.append('0')
                        hc.append('0')
                        body.append(-1)
                        ty.append('A')

natoms = len(ty)
o = open('pre_SCNP_AAM_start.xml', 'w')
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
o.write("<h_init num='%s'>\n" % (len(hi)))
o.write('\n'.join(hi))
o.write("\n</h_init>\n")
o.write("<h_cris num='%s'>\n" % (len(hc)))
o.write('\n'.join(hc))
o.write("\n</h_cris>\n")
o.write("\n</configuration>\n</galamost_xml>")
o.close()
