import argparse
from sys import argv
import numpy as np
import sys
import os
import random
#import tqdm
from XmlParser import XmlParser
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

def pbc(dr,L):
    dr=dr-np.rint(dr/L)*L
    return dr

def dr2(pos1,pos2):
    dr2 = np.sum((pbc(pos1 - pos2,boxl))**2,axis=-1)
    return dr2

def is_in(r2,counti,n_bins):
    for j in range(n_bins):
        condition1 = r2 > (delta_r * j) ** 2
        condition2 = r2 < (delta_r * (j + 1)) ** 2
        condition = condition1 & condition2
        counti[j] += np.sum(condition)
    return counti

def qiucm(pos,box):
    cm = pbc(pos[...,0,:]+pbc(np.diff(pos, prepend=pos[...,:1,:],axis=-2),box).cumsum(axis=-2).mean(axis=-2),box)
    return cm

def a_frame_calculate(pos1,pos2,pos3,pos,n_bins):
    '''pos1:scnp  pos2:shui  pos3:chain'''

    scnp_cm = qiucm(pos1,boxl)
    scnp_r2 = dr2(pos1,scnp_cm)
    h2o_r2 = dr2(pos2, scnp_cm)
    chain_r2 = dr2(pos3, scnp_cm)
    all_r2 = dr2(pos, scnp_cm)
    count_np = np.zeros(n_bins)
    count_h2o = np.zeros(n_bins)
    count_chain = np.zeros(n_bins)
    count_all = np.zeros(n_bins)
    count_np = is_in(scnp_r2,count_np, n_bins)
    count_h2o = is_in(h2o_r2, count_h2o, n_bins)
    count_chain = is_in(chain_r2, count_chain, n_bins)
    count_all = is_in(all_r2, count_all, n_bins)
    return count_h2o,count_chain,count_all,count_np
parse = argparse.ArgumentParser(description='acf')
parse.add_argument('-xml',dest='xml',type=str,nargs='+',help='XML Files')
args = parse.parse_args()
xmls = args.xml

xml = XmlParser(xmls[0],needed=None)
global boxl
with open(xmls[0],'r') as f:
    lines = f.readlines()
    boxl = float((lines[3].split()[1]).split('"')[-2])
print(boxl)
typied1 = (xml.data['type']=='S')         #selct type atom  水
typied2 = (xml.data['type']=='A')#chain
#typied3 = (xml.data['type']=='A')
#shui_pos = xml.data['position'][typied1]         #chain position
#chain_pos = xml.data['position'][typied2]
#N1 = len(cpos)
#N2 = len(npos)
npl = 300


delta_r = 1
rangemax = 20.0
n_bins = int(rangemax/delta_r)
bins_r = delta_r*np.arange(n_bins+1)
Vbins = 4/3*np.pi*((bins_r[1:])**3 - (bins_r[:-1])**3)
r_middle = delta_r*(np.arange(n_bins)+0.5)
d1 = []#水
d2 = []#链
d3 = []#全部
d4 = []#球
n=0
#for i in range(len(xmls)):
for i in range(1,5):
    print(f'{i}th')
    n+=1
    xmli = XmlParser(xmls[-i],needed=None)
    posi = xmli.data['position']
    posi = posi + xmli.data['image']*boxl
    shui_posi = xmli.data['position'][typied1] + xmli.data['image'][typied1]*boxl  # chain position
    chain_posi = xmli.data['position'][typied2] + xmli.data['image'][typied2]*boxl
    pos_scnp = posi[:npl,:] + xmli.data['image'][:300]*boxl
    print(posi.shape)
    print(shui_posi.shape)
    print(pos_scnp.shape)
    print(chain_posi.shape)
    ct_h2o,ct_cl,ct_all,ct_np = a_frame_calculate(pos_scnp,shui_posi,chain_posi,posi,n_bins)
    print(ct_h2o)
    print(ct_cl)
    print(ct_np)
    d1.append(ct_h2o)
    d2.append(ct_cl)
    d3.append(ct_all)
    d4.append(ct_np)
d1 = np.mean(np.asarray(d1),axis=0)/Vbins
d2 = np.mean(np.asarray(d2),axis=0)/Vbins
d3 = np.mean(np.asarray(d3),axis=0)/Vbins
d4 = np.mean(np.asarray(d4),axis=0)/Vbins

#fig = plt.figure()
#plt.plot(r_middle, d1, lw=2, alpha=.6, label='H2O_density')
#plt.plot(r_middle, d2, lw=2, alpha=.6, label='chain_density')
#plt.plot(r_middle, d3, lw=2, alpha=.6, label='total_density')
#plt.plot(r_middle, d4, lw=2, alpha=.6, label='NP_density')
#plt.legend(loc='upper right')
#plt.xlabel('\N{GREEK CAPITAL LETTER DELTA} r',fontsize=19)
#plt.ylabel('number density',fontsize=19)
#plt.legend(fontsize=19)
#fig.set_tight_layout('tight')
#plt.legend()
#plt.savefig('density.png', dpi = 300)
#np.savetxt('den.txt', np.vstack((r_middle, d1, d2, d3, d4)).T)

fig = plt.figure()

# 设置边框样式
for s in plt.gca().spines:
    plt.gca().spines[s].set_linewidth(2)

plt.tick_params('both', length=4, width=2, labelsize=15)

# 修改线条样式
plt.plot(r_middle, d1, lw=2, alpha=.6, label=r'H$_2$O density', color='blue', linestyle='-')
plt.plot(r_middle, d2, lw=2, alpha=.6, label='chain density', color='#FF5733', linestyle='-')
plt.plot(r_middle, d3, lw=2, alpha=.6, label='total density', color='red', linestyle='-')
plt.plot(r_middle, d4, lw=2, alpha=.6, label='SCNP density', color='green', linestyle='-')

# 坐标轴标签与字体大小
plt.xlabel(r'$\Delta r$', fontsize=19)
plt.ylabel('Number density', fontsize=19)

# 调整图例，放在右上角，并调整字体大小
plt.legend(fontsize=14, loc='upper right')

# 调整布局和分辨率
fig.set_tight_layout('tight')
plt.savefig('density.png', dpi=330)

# 保存数据
np.savetxt('den.txt', np.vstack((r_middle, d1, d2, d3, d4)).T)

