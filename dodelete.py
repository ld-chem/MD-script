from sys import argv
from hoomd_script import *
##删除hi hr##
a = init.read_xml(argv[1])


pos = []
ima = []
ty = []
bond = []
body = []
for i in a.bonds:
	bond.append([i.type,i.a,i.b])
Lx,Ly,Lz = a.box.Lx,a.box.Ly,a.box.Lz

for i in a.particles:
    pos.append(i.position)
    ima.append(i.image)
    ty.append(i.type)
 
    body.append(i.body)
    
   
			
natoms = len(ty)
o = open('dd_SCNPPAAM_re_init.xml', 'w')
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
