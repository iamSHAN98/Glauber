#!/usr/bin/python3

import h5py
import matplotlib.pyplot as plt
import util
import datastream as ds

# Command line arguments
args = util.get_arguments()

# Binning parameter
nbin = 15

# Load data
lr = []
nev = 0

for entry in args.input :
    print("Reading data from ", entry)

    with h5py.File(entry, 'r') as file :
        dset_spec = ds.get_data(file, "Glauber/Spec")
        nev += ds.get_event(file, "Glauber/Spec")
    
    lr += list(abs(dset_spec[:,0] - dset_spec[:,1]))
    if nev % 10000 == 0 :
        print(nev, " events read")

print("Total events : " , nev)

# Collision system
with h5py.File(entry, 'r') as file :
    nu_a = ds.get_attribute(file, "Glauber", "NucleusA")
    nu_b = ds.get_attribute(file, "Glauber", "NucleusB")
    snn = ds.get_attribute(file, "Glauber", "sNN")

    nu_a = nu_a['Name'].decode('utf-8')
    nu_b = nu_b['Name'].decode('utf-8')

# Spectator asymmetry distribution
bin_center, bin_count = util.histogram_1d(lr, nbin)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(bin_center, bin_count, ls = '-', lw = 2, c = 'cornflowerblue', 
        marker = None, drawstyle = 'steps-mid')        
ax.fill_between(bin_center, bin_count, color = 'lightgrey', step = 'mid')

# Plot formatting
ax.set_xlabel(r"$|L - R|$", size = 25, labelpad = 5)
ax.set_ylabel(r"P($|L - R|$)", size = 25, labelpad = 5)

ax.set_xlim(0.,)
ax.set_yscale("log")

ax.tick_params(axis = 'both', labelsize = 20)
ax.ticklabel_format(axis = 'x', style = 'sci', scilimits = (-3, 3),
                    useMathText = True)
ax.locator_params(axis = 'x', nbins = 5, integer = False)
ax.minorticks_off()

sys = r"%s-%s, $\sqrt{s}$ = %.0f GeV" % (nu_a, nu_b, snn)
ax.legend(title = sys, title_fontsize = 15, prop = {'size' : 15},
          ncol = 2, framealpha = 1, loc = 'upper right')

# Output : plot
fig.tight_layout()
try :
    ext = args.output[0].split(".")[1]
    fig.savefig(args.output[0], format = ext, bbox_inches = 'tight',
                pad_inches = 0.25)
except :
    plt.show()