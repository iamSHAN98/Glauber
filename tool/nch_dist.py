#!/usr/bin/python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import util
import datastream as ds

# Command line arguments
args = util.get_arguments()

# Binning parameter
nbin = 150

# Storage arrays
nch = []

# Load data
nev = 0
for entry in args.input :
    print("Reading data from ", entry)

    with h5py.File(entry, 'r') as file :
        dset_nch = ds.get_data(file, "Glauber/Nch")
        nev += ds.get_event(file, "Glauber/Nch")

    nch += list(dset_nch)
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

# Calculate : multiplicity distribution
bin_center, bin_count = util.histogram_1d(nch, nbin)

print("\nStarting centrality class determination\n")

# Calculate : centrality classes
cent_per = [0, 5, 10, 15, 20, 30, 40, 50, 60]

cent_data = []
for i in range(len(cent_per) - 1) :
    cent_data.append(util.Centrality(cent_per[i], cent_per[i + 1]))

event_count, index = 0, 0
for i in range(nbin - 1, -1, -1) :

    event_count += bin_count[i]
    cent_data[index].add('nch', bin_center[i])
    per = cent_per[index + 1] - cent_per[index]

    if event_count >= 0.01*per*nev :
        print(cent_data[index].get_name(), "bin done")
        event_count = 0

        index += 1
        if index == len(cent_data) :
            break

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)

# Plot : multiplicity distribution
ax.plot(bin_center, bin_count, ls = '-', lw = 3, c = 'cornflowerblue', 
        marker = None, drawstyle = 'steps-mid', label = "Model")        
ax.fill_between(bin_center, bin_count, color = 'lightgrey', step = 'mid')

# Plot : centrality classes
txt_y = 2*np.std(bin_count)/np.mean(bin_count)
for cent in cent_data :

    try :
        nch_min = cent.get_minimum('nch')
        nch_max = cent.get_maximum('nch')

        txt_x = cent.get_mean('nch') - 0.5*cent.get_error('nch')

        if cent_data.index(cent) % 2 == 0 :
            bound = (bin_center >= nch_min) & (bin_center <= nch_max)
            ax.fill_between(bin_center, bin_count, where = bound, 
                            step = 'mid', color = 'grey')

        ax.text(txt_x, txt_y, cent.get_name(), fontsize = 10, weight = 'bold', 
                color = 'black', rotation = 90)

    except :
        pass

# Plot formatting
ax.set_xlabel(r"$dN_{ch}/d\eta$", size = 25, labelpad = 5)
ax.set_ylabel(r"P($dN_{ch}/d\eta$)", size = 25, labelpad = 5)

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

# Calculate : quantities per centrality class
info = ["\nSystem : %s-%s @ %.0f GeV" % (nu_a, nu_b, snn)]
info.append("# Order : %-10s %-10s %-10s %-10s \n" % ("Centrality",
            "Minimum", "Maximum", "Average"))

# Quantity : multiplicity
info.append("# Multiplicity")
for cent in cent_data :    
    info.append("%-10s %-10.3f %-10.3f %-10.3f" % (cent.get_name(),
        cent.get_minimum('nch'), cent.get_maximum('nch'), cent.get_mean('nch')))

# Quantity : impact parameter
info.append("\n# Impact parameter [fm]")

# 1. Load data
b = []
for entry in args.input :
    with h5py.File(entry, 'r') as file :
        dset_b = ds.get_data(file, "Glauber/b")
    b += list(dset_b)

# 2. Sort into centrality classes
for cent in cent_data :
    cent.add_quantity('b')

    for i in range(nev) :
        event_nch, event_b = nch[i], b[i]
        if (event_nch >= cent.get_minimum('nch')) and (event_nch <= cent.get_maximum('nch')) :
            cent.add('b', event_b)

# 3. Append to output
for cent in cent_data :
    info.append("%-10s %-10.3f %-10.3f %-10.3f" % (cent.get_name(),
        cent.get_minimum('b'), cent.get_maximum('b'), cent.get_mean('b')))

# Output : file
try :
    with open(args.output[1], 'w') as file :
        for msg in info :
            file.write(msg + "\n")
except :
    for msg in info :
        print(msg)