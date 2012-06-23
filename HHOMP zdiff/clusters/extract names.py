import csv
import re

with open('hhomp cluster list.txt', 'r') as f:
    clusternames = list()
    for line in f:
        name = re.search('cluster\d+', line)
        if name is not None:
            clusternames.append(name.group(0))

with open('clusters.csv', 'wb') as f:
    csv.writer(f).writerows([[name] for name in clusternames])

print('done')
