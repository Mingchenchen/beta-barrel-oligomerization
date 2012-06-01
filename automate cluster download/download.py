import urllib
import re
import csv

class UnrecognizedClusterName(Exception):
    pass

url_base = 'http://toolkit.tuebingen.mpg.de/hhomp_ali/'

files = list()

try:
    f = open('clusters.csv', 'r')
    for line in csv.reader(f):
        if len(line) == 0 or line[0] == '':
            continue
        elif re.match('cluster\d+', line[0]) is not None:
            files.append(line[0] + '.clu')
        elif re.match('OMP\..+?\.\d+\.\d+', line[0]) is not None\
             or re.match('NodT\.\d', line[0]) is not None\
             or re.match('TolC\.\d', line[0]) is not None:
            files.append(line[0] + '.clu')            
        elif re.match('.+?\.\d+\.\d+', line[0]) is not None:
            files.append('OMP.' + line[0] + '.clu')
        else:
            raise UnrecognizedClusterName(line[0])
finally:
    f.close()

for filename in files:
    urllib.urlretrieve(url_base + filename, filename=filename)

print('done')
