import pandas as pd
import matplotlib.pyplot as plt

# Integrate stress-strain data as .xlsx #

file_path = 'tensile.xlsx'
df = pd.read_excel(file_path, sheet_name=0)

x = df['X']
y1 = df['Y1']
y2 = df['Y2']

fig = plt.figure()
w, h = fig.get_figwidth(), fig.get_figheight()
fig.set_figwidth(w)

ax = fig.add_subplot(111)
for s in ax.spines:
    ax.spines[s].set_linewidth(2)
ax.tick_params('both', length=4, width=2, labelsize=15)

plt.plot(x, y1, lw=2, alpha=.6, label='SCNP_PAAM', color = 'g', linestyle='-')
plt.plot(x, y2, lw=2, alpha=.6, label='PAAM', color = 'm', linestyle='-')

#plt.title('triaxial tension',fontsize=19)
plt.xlabel('Strain $\lambda$',fontsize=19)
plt.ylabel('Stress $\sigma$',fontsize=19)
plt.legend(fontsize=19)

fig.set_tight_layout('tight')
plt.savefig('tension.png', dpi=330)
plt.legend()