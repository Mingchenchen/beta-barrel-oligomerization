import os
import sys
sys.path.append(os.getcwd()+'/calculator modules')
sys.path.append(os.getcwd()+'/pymol modules')

from useful import one_letter
from useful import CIDict
import csv
import warnings
import numpy as np
from Bio.PDB import PDBParser
from biodata import file_dict
import ezb

structure_files = file_dict('structures with 1qd5', ['aligned_(.*).pdb'])
parser = PDBParser()

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    structures = [(name, parser.get_structure(name, path)) \
                  for name, path in structure_files.items()]
    structures = CIDict(structures)

with open('cored 1 selections with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    inclusive_selections = ezb.selections_by_resi(reader)
    
with open('beta_selections.csv', 'rb') as f:
    reader = csv.reader(f)
    beta_selections = ezb.selections_by_resi(reader)

exclusive_selections = CIDict()
    
for name, selection in beta_selections.items():
   exclusive_selections.update(((name, list(set(selection).intersection( 
                                      set(inclusive_selections[name])))),))
    

with open('cored 1 centers with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    inc_centers = ezb.load_centers(reader)

with open('exc centers.csv', 'rb') as f:
    exc_centers = ezb.load_centers(csv.reader(f))
    
def moments(path, selections, centers):
    with open('published params.csv', 'rb') as f:
        reader = csv.reader(f)
        new_calc = ezb.Calculator(reader, normalize = False)

    moments = CIDict()
    for name in structures.keys():
        if name.upper() == '1QD5':
            continue

        moments.update({name:
                        ezb.moment(structures[name], selections[name],
                                   centers[name], new_calc,
                                   paramless_option = '0',
                                   old_style_gly = True) })

        with open(path, 'wb') as f:
            target = csv.writer(f)
            for name, moment in moments.items():
                row = [name.upper()] + [str(component) for component in moment]
                target.writerow(row)

moments('non normalized inclusive_moment.csv',inclusive_selections, inc_centers)
moments('non normalized exclusive_moment.csv',exclusive_selections, exc_centers)
print('done')
