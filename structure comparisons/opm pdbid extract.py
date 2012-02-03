import re
from urllib import urlretrieve

pdbids = set()

with open('opm beta barrel tm.html', 'r') as f:
    for line in f:
        search = re.search('protein.php\?pdbid=(....)', line)
        if search is not None:
            pdbids.add(search.group(1))

for pdbid in pdbids:
    urlretrieve('http://www.rcsb.org/pdb/download/' +\
                'downloadFile.do?fileFormat=FASTA&compression=NO&' +\
                'structureId=' + pdbid.upper(),
                'bbtm sequences/{}.fasta'.format(pdbid))
