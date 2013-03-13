from __future__ import division
import random
import csv
import math

def mean(sample):
    sample = [float(i) for i in sample]
    return sum(sample)/len(sample)

def sd(sample):
    sample = [float(i) for i in sample]
    return math.sqrt(sum((i-mean(sample))**2 for i in sample)\
                     /(len(sample)-1))

params = list(csv.reader(open('published params.csv')))

# Assign a random curve form to each residue
def rand_curve():
    if random.random() > 4/17:
        return 'sigmoidal'
    else:
        return 'gaussian'
params[1][2:] = [rand_curve() for i in params[1][2:]]

# Generate random E0/Emin's with a normal approximation, and replace
# the real E0/Emin's with them
energies = [float(i) for i in params[2][2:]]
rand_energies = [random.normalvariate(mean(energies), sd(energies))
                 for i in energies]
params[2][2:] = rand_energies

# Generate random zmid/zmins with a normal approximation, and replace
# the real zmid/zmins with them
zs = [float(i) for i in params[3][2:]]
rand_zs = [random.normalvariate(mean(zs), sd(zs)) for i in zs]
params[3][2:] = rand_zs

# Generate random n/sigmas with a normal approximation, and replace the
# real n/sigmas with them
steepnesses = [float(i) for i in params[4][2:]]
rand_steepnesses = [random.normalvariate(mean(steepnesses), sd(steepnesses))
                    for i in steepnesses]
params[4][2:] = rand_steepnesses

# Write the new, meaningless parameters to a file
with open('normal random params.csv', 'wb') as f:
    csv.writer(f).writerows(params)
