
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cmasher as cmr

plt.style.use('http://stephenwilkins.co.uk/matplotlibrc.txt')




# define filters and colours
filters = [f'JWST/NIRCam.{f}' for f in ['F115W', 'F150W', 'F200W', 'F277W', 'F356W', 'F410M', 'F444W']]  # define a list of filter codes
colours = [[f1, f2] for f1, f2 in zip(filters[:-1], filters[1:])]

# define colour colours
colour_colours = cm.rainbow(np.linspace(0, 1, len(colours)))

# read synthesizer predictions
data = {}
for duration in [10,100]:
    data[duration] = pickle.load(open(f'data/{duration}Myr.pck','rb'))
z = data[10]['z']


# intialise figure
fig, axes = plt.subplots(len(colours), 1, figsize = (3.5,6), sharex = True)
plt.subplots_adjust(left=0.15, top=0.975, bottom=0.075, right=0.95, wspace=0.0, hspace=0.0)


for (f1, f2), ax, color in zip(colours, axes, colour_colours):

    for duration, lw, alpha in zip([10,100],[1,2],[1, 0.3]):

        for t, label, ls in zip(['stellar', 'intrinsic'],[r'pure\ stellar', r'stellar\ +\ nebular'],[':','-']):

            c = 2.5*np.log10(data[duration][t][f2]/data[duration][t][f1])

            # plot median line
            ax.plot(z, c, c=color, lw=lw, alpha=alpha, ls=ls, label = rf'$\rm {label}\ {duration}\ Myr$')

    ax.set_xlim([4.,10.])
    ax.set_ylim(-0.99, 1.49)
    ax.set_yticks(np.arange(-0.5, 1.5, 0.5))
    
    ax.set_ylabel(rf"$\rm {f1.split('.')[-1]} - {f2.split('.')[-1]}$", fontsize=7)

axes[-1].set_xlabel(rf"$\rm z$", fontsize=8)
axes[0].legend(loc = 'upper left', fontsize=8, labelspacing=0.0)

fig.savefig(f'figs/colour_evolution.pdf')










