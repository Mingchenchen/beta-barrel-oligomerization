# Packages either from standard library, or available from ppi
import csv
import scipy.stats
import math
import numpy as np

# modules that I wrote
import mfuncs

execfile('family_moments.py')

def mag_var_sheet(params, func):
    output = list()

    fams = families(params=params)
    for pdbid, family in fams.items():
        # Select exposed beta sheet residues
        sel = Selection(family, lambda x: x.rel_acc > .2 and x.ss == 'E')

        # Calculate moments
        moments = sel.all_moments(func)

        # Get magnitudes of each moment
        magnitudes = np.array([np.linalg.norm(moment) for moment in moments])
        # Get angles - all will be in range (-pi, pi)
        angles = np.array([math.atan2(moment[1], moment[0])\
                           for moment in moments])

        mean_mag = np.mean(magnitudes)
        circular_standard_deviation = scipy.stats.circstd(angles, low=-math.pi,
                                                          high=math.pi)

        output.append((pdbid, mean_mag, circular_standard_deviation))

    return output

# Name, parameter sheet, and function for ez-beta moments (made using published
# ez-beta parameters, not the newer asymmetric parameters, though that may have
# been better)
ez_b = ('published params', 'published params.csv', mfuncs.ez_b)
# Likewise, for moments where each residue is assigned a random energy between
# 0 and 1
rand_nrg = ('random energy', 'published params.csv', mfuncs.random_magnitude)
# And for ez-beta moments with parameters that I randomly generated (see March
# 6 entry of my log to understand how)
rand_ezb = ('random ezb parameters', 'normal random params.csv', mfuncs.ez_b)

for name, param_path, func in [ez_b, rand_nrg, rand_ezb]:
    with open('magnitude vs variance plots/{}.csv'.format(name), 'wb') as f:
        csv.writer(f).writerows(mag_var_sheet(param_path, func))
