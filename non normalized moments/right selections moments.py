import sys
sys.path.append(r'C:\cygwin\home\alex\beta barrels\non normalized moments\calculator modules')

from os import chdir
from sundries import one_letter
from useful import CIDict
import csv
import warnings
import numpy as np
from Bio.PDB import PDBParser
from biodata import file_dict
import ezb
reload(ezb)

#chdir('/home/davis/Desktop/work stuff/final paper/figures')
#chdir(r'C:\Users\Nanda Lab\Desktop\Alex\final paper\figures')

# Use Biopython's pdb file parser to create a dictionary called
# "structures" containing Biopython's pdb structure objects,
# created from the files in the folder "structures with 1qd5"
structure_files = file_dict('structures with 1qd5', ['aligned_(.*).pdb'])
parser = PDBParser()

# Daniel's aligned pdb files are missing b factors, and perhaps 
# some other information, which causes Biopython to give hundreds
# or thousands of warnings when loading them. The following code
# suppresses warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    structures = [(name, parser.get_structure(name, path)) \
                  for name, path in structure_files.items()]
    structures = CIDict(structures)


# Load the selections that should have been used

with open('exc centers/exc selections.csv', 'rb') as f:
    reader = csv.reader(f)
    beta_selections = ezb.selections_by_resi(reader)

right_selections = CIDict()
    
for name, selection in beta_selections.items():
   right_selections.update(((name, selection),))

with open('exc centers.csv', 'rb') as f:
    exc_centers = ezb.load_centers(csv.reader(f))
    
def moments(path, selections, centers):
    with open('published params.csv', 'rb') as f:
        reader = csv.reader(f)
        new_calc = ezb.Calculator(reader, normalize = True)

    moments = CIDict()
    for name in structures.keys():
        if name.upper() == '1QD5':
            continue

        moments.update({name:
                        ezb.moment(structures[name], selections[name],
                                   centers[name], new_calc,
                                   paramless_option = '.5',
                                   old_style_gly = True) })

        with open(path, 'wb') as f:
            target = csv.writer(f)
            for name, moment in moments.items():
                row = [name.upper()] + [str(component) for component in moment]
                target.writerow(row)

moments('right selection exc moments.csv',right_selections, exc_centers)
print('done')
