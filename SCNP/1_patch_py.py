# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/bin/env python
from math import copysign
from math import sqrt
from math import sin
from math import cos
from math import pi
from random import random
from random import randint
from random import shuffle
import re
from functools import reduce
from numpy import loadtxt


#
## Functions
#

def pbc(r, d):
    return r - d * round(r / d)


class position(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        return None
    def __str__(self):
        return "(%s, %s, %s)" % (self.x, self.y, self.z)

def icell(ix, iy, iz, ibx, iby, ibz):
    return (ix + ibx) % ibx + ((iy + iby) % iby) * ibx + ((iz + ibz) % ibz) * iby * ibx


def build_box(r_cut, box):
    ibx = int(box.x / r_cut)
    iby = int(box.y / r_cut)
    ibz = int(box.z / r_cut)
    bx = box.x / float(ibx)
    by = box.y / float(iby)
    bz = box.z / float(ibz)
    fb = position(bx, by, bz)
    ib = position(ibx, iby, ibz)
    box_hash = {}
    box_map_hash = {}
    for ix in range(ibx):
        for iy in range(iby):
            for iz in range(ibz):
                imap = icell(ix, iy, iz, ibx, iby, ibz)
                box_hash[imap] = []
                box_map_hash[imap] = []
                for iix in range(-1,2):
                    for iiy in range(-1,2):
                        for iiz in range(-1,2):
                            box_map_hash[imap].append(icell(ix + iix, iy+iiy, iz+iiz, ibx, iby, ibz))
    print(max(box_map_hash.keys()))
    return box_hash, box_map_hash, ib
    

def wcell(x, y, z, ib, box):
    return int((x/box.x + 0.5) * ib.x) + int((y/box.y + 0.5) * ib.y) * ib.x + int((z/box.z + 0.5) * ib.z) * ib.y * ib.x


def distance(a, b, Lx, Ly, Lz):
    return sqrt(pbc(a.x - b.x, Lx) ** 2 + pbc(a.y - b.y, Ly) ** 2 + pbc(a.z - b.z, Lz) ** 2)

def checkbefore_dict(p, imap, box_hash, box_map_hash, check, Lx, Ly, Lz):
    for box in box_map_hash[imap]:
        for pos in box_hash[box]:
            if distance(p, pos, Lx, Ly, Lz) < check:
                return 1
    return 0

def checkbefore(p, f, check, Lx, Ly, Lz):
    for pos in f:
        if distance(p, pos, Lx, Ly, Lz) < check:
            return 1
    return 0


def rand():
    return random() - 0.5
    
def readxyz(file):
    res = []
    line = file.readline()
    while line:
        line = line.strip()
        l = re.split('\s', line)
        res.append(position(float(l[0]), float(l[1]), float(l[2])))
        line = file.readline()
    return res

def randw(radi):
    x = rand()
    y = rand()
    z = rand()
    r = sqrt(x*x + y*y + z*z)
    x /= r
    y /= r
    z /= r
    return position(x * radi, y * radi, z * radi)

from random import randint
#
## Varibles
#
#BoxDim = [ float(x) for x in re.split('\s', input("Box dimention (Lx Ly Lz): ")) ]
BoxDim  = [100,100,100]
BoxVol = reduce((lambda x, y : x * y), BoxDim)
[Lx, Ly, Lz] = BoxDim
box = position(Lx, Ly, Lz)
#NumberDensity = 0.85
#MolLength_AVG = 450
#percent = 0.5
#A = int(MolLength_AVG * percent)
#B = MolLength_AVG - A
A = 0
B = 360
per_ABP = 0.05
per_AAM = 0.2
us = []
nt = 0
ind = []
numd = int(B *per_ABP)
import numpy as np
ids = np.arange(360)
ind_ABP = np.random.choice(ids,18)
print(ind_ABP,len(ind_ABP))
ind_AAM = []
for i in range(72):
    ind_AAM_rand = np.random.randint(0,359)
    while 1:
        if ind_AAM_rand in ind_ABP or ind_AAM_rand in ind_AAM:
            ind_AAM_rand = np.random.randint(0, 359)
        else:
            break
    ind_AAM.append(ind_AAM_rand)
print(ind_AAM)	#这里判断一下ind2不能在ind里面
print(len(ind_AAM))
MolLength_AVG= int(A * 2 + B)
from random import randint
import numpy as np
#ind = np.loadtxt('kkkk.log')+A
#print(np.array(ind)-A)
#MolNum = int(BoxVol * NumberDensity )
#Sels = len(NanoParticleConf)
#GraftMolLength = round(MolLength / pdn)
#GraftSiteNum = round(sigmasqrtn / sqrt(7 * GraftMolLength) * NanoParticleSurf)
#MolNum1 = int(BoxVol * 0.85/MolLength_AVG)
#MolNum2 = MolNum - MolNum1 * MolLength_AVG
MolNum1 = 1
box_hash, box_map_hash, ib = build_box(5.0, box)


# Rock with dispersion


#LdSn = 1/3
#LdSl = 3
#ns = round(GraftSiteNum / (1 + LdSn))
#S = round(GraftSiteNum * GraftMolLength / (1 + LdSn * LdSl) / ns)
#L = round(S * LdSl)
#nl = round(ns * LdSn)


#
## Put NP into box
#


#
## Decorate
#


## Matrix
#
NPCoor = []
Bodies = []
Types = []
GraftCoor = []
Images = []
Bonds = []
pidx = 0
MolCoor = []
mass=""
NanoParticleMorph = []
check = 0.5

m=0
while 1:
    if m>= MolNum1:
        break
    print("%sth Mol / %s" % (m, MolNum1))
    TP = position(0,0,0)
    imap = wcell(TP.x, TP.y, TP.z, ib, box)
    MolCoor.append(TP)
    box_hash[wcell(TP.x, TP.y, TP.z, ib, box)].append(TP)
    if 0 in ind_ABP:
        Types.append('B')
    elif 0 in ind_AAM:
        Types.append('D')
    else:
        Types.append('C')
    Bodies.append(-1)
    mass+="1\n"
    Images.append(position(0,0,0))
    pidx += 1
    n = 1
    mollen = MolLength_AVG
    while 1:
        if n >= mollen :
            break
        step = randw(1)
        if n == 1:
            MRC = position(TP.x + step.x, TP.y + step.y, TP.z + step.z)
        else:
            MRC = position(MRC.x + step.x, MRC.y + step.y, MRC.z + step.z)
        STG = position(0,0,0)
        STG.x = pbc(MRC.x, Lx)
        STG.y = pbc(MRC.y, Ly)
        STG.z = pbc(MRC.z, Lz)
        imap = wcell(STG.x, STG.y, STG.z, ib, box)
        if checkbefore_dict(STG, imap, box_hash, box_map_hash, 0.5, Lx, Ly, Lz) == 1:
            MRC = position(MRC.x - step.x, MRC.y - step.y, MRC.z - step.z)
            continue
        MolCoor.append(STG)
        box_hash[wcell(STG.x, STG.y, STG.z, ib, box)].append(STG)
        Bodies.append(-1)
        if n in ind_ABP :
            Types.append('B')
        elif n in ind_AAM:
            Types.append('D')
        else:
            Types.append('C')
        ix = round(MRC.x / Lx)
        iy = round(MRC.y / Ly)
        iz = round(MRC.z / Lz)
        Images.append(position(ix, iy, iz))
        mass+="1\n"
        Bonds.append(position('polymer', pidx -1, pidx))
        pidx += 1
        n += 1
    m+=1
Bonds.append(position('cross',0,1)) 
# Flush 2 file:


natoms = len(Types)
o = open('out.xml', 'w')
o.write("<?xml version = '1.0' encoding = 'UTF-8'?>\n<hoomd_xml version = '1.4'>\n<configuration time_step = '0' dimensions = '3' natoms = '%s'>\n<box lx = '%s' ly = '%s' lz = '%s'/>\n<position num = '%s'>\n" % (natoms, Lx, Ly, Lz, natoms))
o.write('\n'.join(["%s %s %s" % (str(p.x), str(p.y), str(p.z)) for p in NPCoor + GraftCoor + MolCoor]))
o.write("\n</position>\n")
o.write("<body num='%s'>\n" % (natoms))
o.write('\n'.join([ str(x) for x in Bodies ]))
o.write("\n</body>\n")
o.write("<type num='%s'>\n" % (natoms))
o.write('\n'.join(Types))
o.write("\n</type>\n")
o.write("<mass num='%s'>\n" % (natoms))
o.write(mass)
o.write("\n</mass>\n")
o.write("<image num='%s'>\n" % (natoms))
o.write('\n'.join([ "%s %s %s" % (str(p.x), str(p.y), str(p.z)) for p in Images ]))
o.write("\n</image>\n")
o.write("<bond num='%s'>\n" % (len(Bonds)))
o.write('\n'.join([ "%s %s %s" % (str(p.x), str(p.y), str(p.z)) for p in Bonds ]))
o.write("\n</bond>")
o.write("\n</configuration>\n</hoomd_xml>")
o.close()
