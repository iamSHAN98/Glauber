#!/usr/bin/python3

import h5py
import matplotlib.pyplot as plt
import util
import datastream as ds

# Command line arguments
args = util.get_arguments()

# Binning parameter
nbin = 20
hlen = 5

# Load data
nch, en = [], [ [] for i in range(hlen) ]
nev = 0

for entry in args.input :
    print("Reading data from ", entry)

    with h5py.File(entry, 'r') as file :
        dset_nch = ds.get_data(file, "Glauber/Nch")
        dset_en = ds.get_data(file, "Glauber/Ecc")
        nev += ds.get_event(file, "Glauber/b")

    nch += list(dset_nch)
    for i in range(hlen) :
        en[i] += list(dset_en[:,i])
    
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

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)

color = ['red', 'gold', 'limegreen', 'navy', 'fuchsia']

for i in range(hlen) :
    x, y, err = util.histogram_profile(nch, en[i], nbin)
    ax.errorbar(x, y, yerr = err, marker = 'o', ms = 8, mfc = color[i],
            mew = 1, mec = 'black', elinewidth = 1, capsize = 5,
            ecolor = 'black', lw = 0, label = r"$\varepsilon_{%s}$" % (i + 2))

# Plot formatting
ax.set_xlabel(r"$dN_{ch}/d\eta$", size = 25, labelpad = 5)
ax.set_ylabel(r"$\varepsilon_n$", size = 25, labelpad = 5)

ax.tick_params(axis = 'both', labelsize = 20)
ax.ticklabel_format(axis = 'both', style = 'sci', scilimits = (-3, 3),
                    useMathText = True)
ax.locator_params(axis = 'both', nbins = 5, integer = False)
ax.minorticks_off()

sys = r"%s-%s, $\sqrt{s}$ = %.0f GeV" % (nu_a, nu_b, snn)
ax.legend(title = sys, title_fontsize = 15, prop = {'size' : 15},
          ncol = 2, framealpha = 1, loc = 'upper right')

sys = r"%s-%s, $\sqrt{s}$ = %.0f GeV" % (nu_a, nu_b, snn)
ax.legend(title = sys, title_fontsize = 15, prop = {'size' : 15}, 
          ncol = 2, framealpha = 1, loc = 'best')

# Output
fig.tight_layout()
try :
    ext = args.output[0].split(".")[1]
    fig.savefig(args.output[0], format = ext, bbox_inches = 'tight',
                pad_inches = 0.25)
except :
    plt.show()