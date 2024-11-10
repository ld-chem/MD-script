import glob
import re
import sys
from tqdm import tqdm

input_file = '/home/lmy/ld/202410_hydrogel/PAAM_preparation/reaction/PAAM*'
output_file = 'output.txt'
file_paths = glob.glob(input_file)

column_headers = ['file_name', 'bond', 'monomer_index_count', 'reacted_monomer_count']

with open(output_file,'w') as outfile:
    outfile.write('  '.join(column_headers) + '\n')

    for file_path in tqdm(file_paths, desc='Processing files', unit='file'):
        with open(file_path,'r') as infile:
            lines = infile.readlines()

            natoms = None
            bond = None
            bond_index = set()
            bond_count = 0
            recording_bond_index = False

            for i, line in enumerate(lines):
                if 'position' in line and natoms is None:
                    match = re.search(r'(\d+)', line)
                    if match:
                        natoms = int(match.group(1))

                if 'bond' in line:
                    bond_count += 1

                    if bond_count == 1:
                        match = re.search(r'(\d+)', line)
                        if match:
                            bond = int(match.group(1))

                        if i + 1 < len(lines):
                            next_line = lines[i + 1].strip().split()
                            if len(next_line) == 3:
                                bond_index.add(int(next_line[1]))
                                bond_index.add(int(next_line[2]))

                        recording_bond_index = True

                    elif bond_count == 2:
                        break

                if recording_bond_index and bond_count == 1:
                    next_line = line.strip().split()
                    if len(next_line) == 3:
                        bond_index.add(int(next_line[1]))
                        bond_index.add(int(next_line[2]))

            if natoms is not None:
                full_index = list(range(natoms))
                monomer_index = [i for i in full_index if i not in bond_index]
                monomer_index_count = len(monomer_index)
            else:
                monomer_index_count = "N/A"

            file_name = file_path.split('/')[-1]
            bond_value = bond if bond is not None else "N/A"

            if isinstance(monomer_index_count, int):
                reacted_monomer_count = (59000 + 59000 + 59 * 2) - monomer_index_count
                reacted_monomer_percentage = (reacted_monomer_count / 59118) * 100
            else:
                reacted_monomer_percentage = "N/A"

            outfile.write(f"{file_name}  {bond_value}  {monomer_index_count}  {reacted_monomer_percentage:.2f}%\n")
        
