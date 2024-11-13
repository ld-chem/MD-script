import numpy as np
from sys import argv
import glob
from tqdm import tqdm
import re

def pbc(pos, box):
    return pos - box * np.round(pos / box)

input_folder = argv[1]
output_file = 'output.txt'
file_paths = glob.glob(f"{input_folder}")

with open(output_file, 'w') as outfile:
    for file_path in tqdm(file_paths, desc='Processing files', unit='file'):
        with open(file_path, 'r') as infile:
            lines = infile.readlines()

        start_reading_pos = False
        start_reading_bond = False
        pos = []
        bond = []
        bond_vec = []
        box = []

        for line in lines:
            if 'box' in line and not box:
                box_values = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                box = list(map(float, box_values[:3]))

        for line in lines:
            if '<position' in line:
                start_reading_pos = True
                continue
            if '</position' in line and start_reading_pos:
                break
            if start_reading_pos:
                parts = line.strip().split()
                if len(parts) >= 3:
                    x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                    pos.append([x, y, z])

        pos = np.array(pos)
        box = np.array(box)
        for line in lines:
            if '<bond' in line:
                start_reading_bond = True
                continue
            if '</bond' in line and start_reading_bond:
                break
            if start_reading_bond:
                parts = line.strip().split()
                if len(parts) >= 3:
                    particle1 = int(parts[1])
                    particle2 = int(parts[2])
                    bond.append((particle1, particle2))
                    coord1 = np.array(pbc(pos[particle1],box))
                    coord2 = np.array(pbc(pos[particle2],box))
                    vector = coord2 - coord1
                    bond_vec.append(vector)

#        for item in bond:
#            outfile.write(f"{item[0]} {item[1]}\n")
#        for coord in pos:
#            outfile.write(f"{coord[0]} {coord[1]} {coord[2]}\n")       
#        for vec in bond_vec:
#            outfile.write(f"{vec[0]} {vec[1]} {vec[2]}\n")

        bond_vec = np.array(bond_vec)
        x_axis = np.array([1, 0, 0])
        cos_theta = np.dot(bond_vec, x_axis) / np.linalg.norm(bond_vec, axis=1)
        P2 = 0.5 * (3 * cos_theta**2 - 1)

        P2_mean = np.mean(P2)
#        angles = np.degrees(np.arccos(cos_theta))
        outfile.write(f"{file_path} {P2_mean}\n")
