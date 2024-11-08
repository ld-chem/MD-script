import pickle
from sys import argv
import networkx as nx
import numpy as np
from xml_parser import XmlParser
from cg_topology import read_cg_topology
from scipy.spatial import cKDTree
import tqdm
from itertools import combinations

def get_combinations(n):
    indices = list(range(n))
    return list(combinations(indices, 3))

f = argv[1]
xml = XmlParser(f)
types = set()
scnp_hash = []
for t in xml.data['type']:
    types.add(t)
    if t in ['A','S']:
        scnp_hash.append(False)
    else:
        scnp_hash.append(True)
mols = {t:{'smiles':'C','file':None} for t in types}
cg, cg_mols = read_cg_topology(xml,mols)
scnp = cg_mols[0]
cs = nx.cycle_basis(scnp)
position = xml.data['position'][xml.data['type'] == 'A']
idx = np.arange(len(xml.data['type']))[xml.data['type'] == 'A']
CM = xml.data['position'][scnp_hash].mean(axis=0)
tree = cKDTree(position)
distance = 5
indices = tree.query_ball_point(CM, distance)
p_CM = position[indices]
p_CM_idx = idx[indices]
cross = set()
for c in cs:
    N = len(c)
    combination = get_combinations(N)
    for com in tqdm.tqdm(combination,total=len(combination)):
        i,j,k = com
        xi,xj,xk = scnp.nodes[c[i]]['x'],scnp.nodes[c[j]]['x'],scnp.nodes[c[k]]['x']
        ij = np.array(xi-xj)
        jk = np.array(xj-xk)
        n = np.cross(ij,jk)
        for cid,pi in zip(p_CM_idx ,p_CM):
            planeOcid = np.dot(n,(pi - xi))
            neis = cg_mols[1].neighbors(cid)
            for nei in neis:
                if (cid,nei) in cross:
                    continue
                planeOnid = np.dot(n,(cg_mols[1].nodes[nei]['x'] - xi))
                if planeOnid*planeOcid < 0:
                    cross.add((cid,nei))
pickle.dump(cross,open('cross_chain.pkl','wb'))
cross = pickle.load(open('cross_chain.pkl','rb'))
print(cross)
