from tqdm import tqdm
import re
import sys
import numpy as np
import glob

def pbc(r, d):
    return r - d * np.round(r / d)

def unwrap_pbc(r, d):
    ref_atom = r[0]
    unwrapped_traj = np.copy(r)
    n_atoms, _ = r.shape

    for j in range(n_atoms):
        delta = unwrapped_traj[j] - ref_atom
        delta -= np.round(delta / d) * d
        unwrapped_traj[j] = ref_atom + delta

    return unwrapped_traj

def qiucm(pos, box):
    dpos = np.append(np.zeros(((1, 3))), pbc(np.diff(pos, axis=0), box), axis=0)
    rpos = np.cumsum(dpos, axis=0)
    cm = rpos.mean(axis=0) + pos[0]
    cm = pbc(cm, box)
    return cm

def qiurg(pos, cm):
    rg_squared = np.mean(np.sum((pos - cm) ** 2, axis=1))
    rg = np.sqrt(rg_squared)
    return rg

def qiurg_axes(pos, cm):
    rg_squared_axes = np.mean((pos - cm) ** 2, axis=0) 
    rg_axes = np.sqrt(rg_squared_axes) 
    return rg_axes

def readpos(file_path, num=300):
    pos = []
    box = []
    start_reading = False
    count = 0

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for line in lines:
            if 'box' in line and not box:
                box_values = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                box = list(map(float, box_values[:3]))

            if 'position' in line:
                start_reading = True
                continue

            if start_reading and count < num:
                coords = list(map(float, line.strip().split()))
                pos.append(coords)
                count += 1

            if count >= num:
                break

    return np.array(pos), np.array(box)

input_files = glob.glob(sys.argv[1])
output_file = 'output_Rgvec.txt'

with open(output_file, 'w') as outfile:
    outfile.write("file_name\tRg\tRg_x\tRg_y\tRg_z\n")

    for file_path in tqdm(input_files, desc="Processing files", unit="file"):
        pos, box = readpos(file_path)

        qiu = 300
        qiupos = pos[:min(qiu, len(pos))]

        qiupos = unwrap_pbc(qiupos, box)

        cm = qiucm(qiupos, box)
        rg = qiurg(qiupos, cm)  
        rg_axes = qiurg_axes(qiupos, cm)  

        outfile.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\n" % (file_path, rg, rg_axes[0], rg_axes[1], rg_axes[2]))
