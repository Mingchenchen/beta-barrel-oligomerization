import re

def extract_pdbids_from_opm(opm_page):
    '''Given an OPM page, returns all PDBIDs on that page.
    OPM organizes its entries by family, and it's a good place to get
    an up-to-date list of known transmembrane beta barrels.
    Input should be an html file, output will be a list of strings.'''
    pdbids = list()
    for line in opm_page:
        search = re.search('protein.php\?pdbid=(....)', line)
        if search is not None:
            pdbid = search.group(1)
            # To avoid redundancy:
            if pdbid not in pdbids:
                pdbids.append(pdbid)
    return pdbids

with open('opm bbtm 104.htm') as f:
    x = extract_pdbids_from_opm(f)

print('first entry: ' + x[0])
print('fifth entry: ' + x[4])
print('last entry: ' + x[-1])
print('number of entries: ' + str(len(x)))
