import numpy as np
import pandas as pd
from sys import argv

# get the output from uniaxial or triaxial tensile simulation #

data = pd.read_csv(argv[1], delim_whitespace=True)

lx  = data['lx'].values
pxx = data['pressure_xx'].values
pyy = data['pressure_yy'].values
pzz = data['pressure_zz'].values

hydrostatic_press = (pxx + pyy + pzz) / 3
init_lx = lx[0]
tensile_strain = (lx - init_lx) / init_lx

stress = 3*(-pxx + hydrostatic_press) / 2


np.savetxt('stress-strain.txt', np.c_[tensile_strain, stress], fmt='%.6f')