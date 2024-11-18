import numpy as np
import glob
import re
import numba as nb
from tqdm import tqdm
from sys import argv

@nb.njit
def pbc(pos, box):
    """Apply periodic boundary conditions."""
    return pos - box * np.round(pos / box)

@nb.njit
def calculate_bond_vectors(pos, bonds, box):
    """Calculate bond vectors with PBC applied."""
    bond_vec = np.zeros((len(bonds), 3))
    for i in range(len(bonds)):
        particle1, particle2 = bonds[i]
        coord1 = pbc(pos[particle1], box)
        coord2 = pbc(pos[particle2], box)
        bond_vec[i] = coord2 - coord1
    return bond_vec

@nb.njit
def compute_p2(bond_vec):
    """Compute P2 order parameter."""
    x_axis = np.array([1, 0, 0])
    cos_theta = np.dot(bond_vec, x_axis) / np.linalg.norm(bond_vec, axis=1)
    P2 = 0.5 * (3 * cos_theta**2 - 1)
    return np.mean(P2)

def parse_input_file(file_path):
    """Parse box, positions, and bonds from a single input file."""
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    box = []
    pos = []
    bonds = []

    # Parse box dimensions
    for line in lines:
        if 'box' in line and not box:
            box_values = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            box = np.array(list(map(float, box_values[:3])))
            break

    # Parse positions
    start_reading_pos = False
    for line in lines:
        if '<position' in line:
            start_reading_pos = True
            continue
        if '</position' in line:
            break
        if start_reading_pos:
            parts = line.strip().split()
            if len(parts) >= 3:
                pos.append([float(parts[0]), float(parts[1]), float(parts[2])])
    pos = np.array(pos)

    # Parse bonds
    start_reading_bond = False
    for line in lines:
        if '<bond' in line:
            start_reading_bond = True
            continue
        if '</bond' in line:
            break
        if start_reading_bond:
            parts = line.strip().split()
            if len(parts) >= 3:
                bonds.append((int(parts[1]), int(parts[2])))
    bonds = np.array(bonds)

    return box, pos, bonds

def main(input_folder, output_file):
    """Main function to process all files and compute P2 order parameters."""
    file_paths = glob.glob(f"{input_folder}/*")
    
    with open(output_file, 'w') as outfile:
        for file_path in tqdm(file_paths, desc='Processing files', unit='file'):
            # Parse data from the file
            box, pos, bonds = parse_input_file(file_path)
            
            # Compute bond vectors and P2
            bond_vec = calculate_bond_vectors(pos, bonds, box)
            P2_mean = compute_p2(bond_vec)

            # Write results to output file
            outfile.write(f"{file_path} {P2_mean:.6f}\n")

if __name__ == "__main__":
    input_folder = argv[1]
    output_file = "output.txt"
    main(input_folder, output_file)
