one_letter = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}

class CIDict(dict):
    def __init__(self, raw = dict()):
        raw = dict(raw)
        data_to_be = []
        for key, value in raw.items():
            data_to_be.append((key.lower(), value))
        dict.__init__(self, data_to_be)
        
    def __getitem__(self, key):
        return dict.__getitem__(self,key.lower())

    def update(self, iterable):
        given = dict(iterable)
        for key, value in given.items():
            dict.update(self, ((key.lower(), value),))


