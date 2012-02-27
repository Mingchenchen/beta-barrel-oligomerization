from __future__ import division
import csv
import numpy as np
from sundries import CIDict
from math import atan2
from math import pi

def direction(vector):
    return atan2(vector[1], vector[0])

# Save lists of moments to two dictionaries mapping pdbids to numpy arrays
moment_dicts = list()
for filename in ('exclusive_moment.csv',
                 'non normalized exclusive_moment.csv'):
    with open(filename, 'rb') as f:
        reader = csv.reader(f)
        moment_dict = dict()
        for line in reader:
            pdbid = line[0]
            moment = np.array([float(line[1]), float(line[2])])
            moment_dict.update([(pdbid, moment)])
        moment_dicts.append(moment_dict) 

normalized = moment_dicts[0]
non_normalized = moment_dicts[1]

# Put together the data that'll go in the csv file
differences = list()
# entries to put at the beginning, rather than the end, of the list:
priority = ['1A0S', '1QD6', '1AF6', '2O4V', '2J1N', '2POR', '1E54', '3PRN']

for pdbid in normalized.keys():
    angle_diff = (abs(direction(normalized[pdbid]) \
                    - direction(non_normalized[pdbid]))) * 180/pi
    normalized_magnitude = np.linalg.norm(normalized[pdbid])
    non_normalized_magnitude = np.linalg.norm(non_normalized[pdbid])
    magnitude_diff = normalized_magnitude - non_normalized_magnitude
    # What's going to be added to the list:
    entry = [pdbid, angle_diff, normalized_magnitude,
             non_normalized_magnitude, magnitude_diff,
             magnitude_diff/normalized_magnitude]
    if pdbid.upper() in priority:
        differences.insert(0,entry)
    else:
        differences.append(entry)

# Write data to csv file
with open('moments compare.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerow(['PDBID', 'Angle between moments',
                    'Mag with normalization', 'Mag without normalization',
                    'Difference in magnitudes',
                    'Difference over normalized magnitude'])
    writer.writerows(differences)


