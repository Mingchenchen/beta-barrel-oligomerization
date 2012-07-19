import csv
import numpy as np
from sundries import CIDict
path = r'C:\cygwin\home\alex\beta barrels\omp ppi paper\displaying moments'

def coord_dict(reader):
    return CIDict((pdbid, np.array([float(x),float(y),float(z)])) \
                    for pdbid, x, y, z in reader)

with open(path + '/exc centers.csv', 'rb') as f:
    centers = coord_dict(csv.reader(f))

with open(path + '/exclusive_moment.csv', 'rb') as f:
    moments = coord_dict(csv.reader(f))
    
for pdbid, group in groupdict.items():
    cmd.pseudoatom(pdbid.upper() + '.center',
                   pos = tuple(centers[pdbid] + np.array([0, 0, 25])))
    cmd.pseudoatom('temp2',
                   pos = tuple(centers[pdbid] + moments[pdbid] + \
                               np.array([0,0,25])))
    cmd.dist(pdbid.upper() + '.exc_moment', pdbid + '.center', 'temp2')
    cmd.delete('temp*')
cmd.show('spheres', '*.center')
