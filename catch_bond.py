import glob
import re 
from tqdm import tqdm
import sys

#input_folder = sys.argv[1]
#output_file = sys.argv[2]
input_folder = './PAAM.*'
output_file = 'output.txt' 
file_paths = glob.glob(input_folder)

column_headers = ['time','bond']

with open(output_file,'w') as outfile:
    outfile.write('  '.join(column_headers) + '\n')
    
    for file_path in tqdm(file_paths, desc='Processing files', unit='file'):
        with open(file_path,'r') as infile:
            lines = infile.readlines()

            if len(lines) >=3:
                time = lines[3-1]
                step = re.findall(r'\d+\.?\d*', time)
                if step:
                    first_num = step[0]
                    bond_numbers = []
                    for line in lines:
                        if 'bond' in line.lower():
                            bond_numbers = re.findall(r'\d+\.?\d*', line)
                            break
                    if bond_numbers:
                        outfile.write(f"{first_num} {' '.join(bond_numbers)}\n")
