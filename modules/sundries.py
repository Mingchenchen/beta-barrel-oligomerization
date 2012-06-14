import os
import re

class CIDict(dict):
    def __init__(self, raw = dict()):
        raw = dict(raw)
        data_to_be = []
        for key, value in raw.items():
            data_to_be.append((key.lower(), value))
        dict.__init__(self, data_to_be)
        
    def __getitem__(self, key):
        return dict.__getitem__(self,key.lower())

    def __setitem__(self, key, value):
        return dict.__setitem__(self, key.lower(), value)

    def __delitem__(self, key):
        return dict.__delitem__(self, key.lower())

    def update(self, iterable):
        given = dict(iterable)
        for key, value in given.items():
            dict.update(self, ((key.lower(), value),))


one_letter = CIDict({'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'})
 
three_letter = CIDict((value, key.upper()) \
                      for key, value in one_letter.items())
    
def file_dict(path, filename_templates, key_template = '(*)'):
    """Finds file paths for all files in a folder and its subdirectories whose
    names match given templates. Creates a dictionary mapping keys to
    file paths based on a template given for the keys.
    ARGUMENTS:
    path to search
    list of regular expressions to serve as filename templates
    a string with one or more (*)'s in it: this will serve as the template for
        the keys. Each (*) will be replaced with a group from the filename
        templates. Defaults to '(*)'
    RETURNS:
    A dictionary."""
    output = dict()
    for index, string_ in enumerate(filename_templates):
        filename_templates[index] = re.compile(string_)
    for tuple_ in os.walk(path):
        currentdir = tuple_[0]
        for filename in tuple_[2]:
            for rexp in filename_templates:
                match = rexp.match(filename)
                if match is not None:
                    groups = match.groups()
                    key = ''
                    pieces = key_template.split('(*)')
                    for index, segment in enumerate(pieces):
                        key += segment
                        # Possibility for error here if the number of
                        # (*)'s in key_template is more than the number of
                        # .*'s in the filename template.
                        if index != len(pieces) - 1:
                            key += groups[index]
                    output.update({key: currentdir + '/' + filename})
    return output
