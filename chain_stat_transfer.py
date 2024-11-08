import numpy as np
import ast

# Convert the format to a form that Ovito can recognize #
# The file "data.txt" was gained from chain_stat.py #

with open("data.txt", 'r') as file:
    content = file.read()

data = ast.literal_eval(f"[{content}]")

data1 = []
data2 = []

for i in range(len(data)):
    data1.append(data[i][0])
    data2.append(data[i][1])

datau1 = np.unique(data1)
datau2 = np.unique(data2)
combined_datau = np.concatenate((datau1, datau2))
datau = np.unique(combined_datau)

n_scnp = 300
# Additional particles that need to be input #

with open('output.txt','w') as file:
    for i in range(n_scnp):
        file.write(f"ParticleIndex == {i} || ")

    for index, item in enumerate(datau):
        if index < len(datau) - 1:
            file.write(f"ParticleIndex == {item} || ")
        else:
            file.write(f"ParticleIndex == {item}")