import glob
import re
from tqdm import tqdm
import sys

input = sys.argv[1]
output_file = 'init_catch.txt'
output_modified_file = 'SCNPPAAM_re_init.xml'

with open(output_file,'w') as outfile:
   for file_path in glob.glob(input):
        with open(file_path, 'r') as infile:
            lines = infile.readlines()
 
        start = False
        line_numbers = [] 
        index = []

        for i,line in enumerate(lines):
            if 'init' in line:
                if start:
                    break
                start = True
                index_zero = i + 1
                continue
            
            if start and '1' in line:
                line_numbers.append(i + 1)
                index.append(i - index_zero)
   
        outfile.write(f"{line_numbers}\n{len(line_numbers)}\n{index}\n")

for file_path in glob.glob(input):
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    start_modifying = False
    zero_line_index = None

    for i, line in enumerate(lines):
        if 'type' in line:
            zero_line_index = i + 1
            start_modifying = True
            break

    if zero_line_index is not None and index:
        for idx in index:
            target_line = zero_line_index + idx
            if target_line < len(lines):
                lines[target_line] = re.sub(r'\bA\b', 'A1', lines[target_line])

        for i in range(zero_line_index, len(lines)):
            if 'type' in lines[i]:
                break

    with open(output_modified_file, 'w') as modified_file:
        modified_file.writelines(lines)
