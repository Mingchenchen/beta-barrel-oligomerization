#should have like 617 rows

import csv
import re

with open('hhomp cluster list.txt', 'r') as f:
    clusternames = list()
    for line in f:
        name = re.search('OMP\..+?\.\d+\.\d+', line)
        if name is None:
            name = re.search('TolC\.\d', line)
        if name is None:
            name = re.search('NodT\.\d', line)
        if name is not None:
            clusternames.append(name.group(0))

with open('clusters.csv', 'w') as f:
    csv.writer(f).writerows([[name] for name in clusternames])

print('done')
