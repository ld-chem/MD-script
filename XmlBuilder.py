# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 21:40:33 2021

@author: Greenbridge
"""

import numpy as np
from io import StringIO
from xml.etree import cElementTree  # pypy will be a bit slower than python
from pandas import read_csv
import sys
sys.dont_write_bytecode = True

class hoomd_xml(object):
    '''
    usage: 
        1. xml = hoomd_xml("particle000.xml",["position","bond"]) # put needed attrib in the args 
        2. Pos = xml.nodes["position"]  # get a multidimension numpy array shape(atoms, columns), for positon it is (atoms,3)
        3. box = xml.box
    '''

    @staticmethod
    def _get_attrib(dd):
        dt = eval('[' + ','.join(["('%s', int)" % key for key in dd.keys()]) + ']')
        values = [tuple(dd.values())]
        return np.array(values, dtype=dt)

    def __init__(self, filename, needed=[]):
        tree = cElementTree.ElementTree(file=filename)
        root = tree.getroot()
        configuration = root[0]
        self.configure = self._get_attrib(configuration.attrib)
        self.nodes = {}
        for e in configuration:
            if e.tag == 'box':
                #print(e.attrib)
                self.box = np.array([ float(e.attrib['lx']), float(e.attrib['ly']), float(e.attrib['lz']) ])
                continue
            if (len(needed) != 0) and (not e.tag in needed):
                continue
            self.nodes[e.tag] = read_csv(StringIO(e.text), delim_whitespace=True, squeeze=1, header=None).values



def bond_hash_dualdirect(bond, natoms): 
    '''
    :param bond: bond data in hoomdxml format (name, id1, id2)
    :param natoms: total number of particles
    :return: hash table of natoms keys with value in {bondname1: [idxes], bondname2:[idxes]...} for each particle (in dual direct)
    '''
	#b_t_name = list(set(list(bond[:,0])))
    bond_hash_wn = {i:{} for i in range(natoms)}
    bond_hash_nn = {i:[] for i in range(natoms)}
    for b in bond:
        nm = b[0]
        idx = int(b[1])
        jdx = int(b[2])
        bond_hash_wn[idx][jdx]=nm
        bond_hash_nn[idx].append(jdx)
        bond_hash_nn[jdx].append(idx)
        bond_hash_wn[jdx][idx]=nm

    #print('Done.')
    return bond_hash_nn, bond_hash_wn
	
def grab_iter(i, bond_hash, mol_used, body_hash=None):#from shirui
    #int of mol_used = {i:False for i in range(natoms)}
    #dict((i,False) for i in range(natoms))
    #get mol list(contain atom idx in list)
    S = [i]
    R = []
    while S:
        v = S.pop()
        if mol_used[v] == False:
            R.append(v)
            mol_used[v] = True
            for w in bond_hash[v]:
                S.append(w)
            if not body_hash:
                continue
            for w in body_hash[v]:
                S.append(w)
                for x in bond_hash[w]:
                    S.append(x)
    return R
	
def get_mol_idxes(natoms,bond_hash,body_hash=None):
    #the is modifed and add by lsj based on shirui's hoomd_mols
    molecular_hash = {i:False for i in range(natoms)}
    molecular_list=[]
    for i in range(natoms):
        molecular_idxs = grab_iter(i, bond_hash, mol_used=molecular_hash, body_hash=body_hash)
        if not len(molecular_idxs)==1:
            molecular_list.append(molecular_idxs)
        for atom in molecular_idxs:
            molecular_hash[atom] = True
    while [] in molecular_list:
        molecular_list.remove([])
    return(np.array([ np.array(x) for x in molecular_list ]))	