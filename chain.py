from xml_parser import XmlParser
import numpy as np
from sys import argv
import pandas as pd

f = argv[1]
xml = XmlParser(f,needed = ['type', 'bond'])
types = set()
chain_hash = []
for t in xml.data['type']:
    types.add(t)
    if t in ['A']:
        chain_hash.append(True)
natoms = len(types)
bond_info = (xml.data['bond']=='polymer')
bond_info = bond_info[bond_info[:,0] == 'polymer']

idx = np.arange(len(xml.data['type']))[xml.data['type'] == 'A']

len_mol_fre = np.array([len(mol) for mol in idx])

length = []
for t in types:
    _mass = 1
    length.append(_mass)
length = np.asarray(length)
for i, mol in enumerate(idx):
    mol_mass = np.array([length[index] for index in mol])

denominator = len(len_mol_fre)
Data = pd.Series(len_mol_fre)
Fre = Data.value_counts()
Fre_sort = Fre.sort_index(axis=0, ascending=True)
DF = Fre_sort.reset_index()
DF[0] = DF[0] / denominator
DF = DF[DF[0] != 0]
DF.to_csv('chain.txt', index=0, header=None, sep=' ')

sorted_indices = np.argsort(len_mol_fre)
with open('output.txt','w') as file:
    file.write("Chain Index: {}, Particle Indices: {}".format(sorted_indices, mol_index[sorted_indices]))

