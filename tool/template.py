#!/usr/bin/python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import util
import datastream as ds

# Command line arguments
args = util.get_arguments()

# Binning parameter
nbin = 25

"""
    Declare storage array & load data
"""
nev = 0

for entry in args.input :
    print("Reading data from ", entry)

"""
    1. Read data from file

        with h5py.File(entry, 'r') as file :
            dset = ds.get_data(file, "Glauber/...")
            nev += ds.get_event(file, "Glauber/...")
    
    2. Apply numerics on raw data (if required) & append 
       to storage array

       arr += list(...)

"""

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

"""
    Apply binning : use util.histogram_1d or util.histogram_profile
"""

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)

"""
    Add type of plot : usually scatter, line or errorbar
"""

# Plot formatting
ax.set_xlabel(r"x", size = 25, labelpad = 5)
ax.set_ylabel(r"y", size = 25, labelpad = 5)

ax.tick_params(axis = 'both', labelsize = 20)
ax.ticklabel_format(axis = 'both', style = 'sci', scilimits = (-3, 3),
                    useMathText = True)
ax.locator_params(axis = 'both', nbins = 5, integer = False)
ax.minorticks_off()

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