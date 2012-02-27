import sys
sys.path.append('/home/davis/Desktop/work stuff/pymol/scripts/modules')
sys.path.append('/home/davis/Desktop/work stuff/calculator/modules')

from os import chdir
from useful import one_letter
from useful import CIDict
import csv
import warnings
import numpy as np
from Bio.PDB import PDBParser
from biodata import file_dict
import ezb

chdir('/home/davis/Desktop/work stuff/final paper/figures')
#chdir(r'C:\Users\Nanda Lab\Desktop\Alex\final paper\figures')

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

moments('inclusive_moment.csv',inclusive_selections, inc_centers)
moments('exclusive_moment.csv',exclusive_selections, exc_centers)


# Now the eisenberg stuff:

eisenberg_values = CIDict()
with open('eisenberg.csv', 'rb') as f:
    freader = csv.reader(f)    
    for row in freader:
        eisenberg_values.update({row[0]: float(row[1])})
        
def eisenberg(residue, ref = eisenberg_values):
    return ref[residue.get_resname()]
        
def flexible_moments(path, selections, centers, function):
    moments = list()
    for name in structures.keys():
        
        if name.lower() == '1qd5':
            continue

        structure = structures[name]
        center = centers[name]
        selection = selections[name]
        
        sum_ = np.zeros(3)
        for residue in structure.get_residues():
            if residue.get_id()[1] not in selection:
                continue

            # My selection files just give residue numbers. If two residues
            # have the same number, if there're insertions, then the selection
            # files are ambiguous and I need to change how I do them.
            # For now, anything identified by more than just a residue number
            # is ignored. The reason this is okay is that, as far as I can tell,
            # in Dan's structures that just means it ignores water.
            # If the structures have complex residue ids, then I need different
            # selection files.
            if residue.get_id()[0] != ' ':
                if residue.get_id()[0] != 'W':
                    print('Ignored residue with id {0} in structure {1}'\
                           .format(residue.get_id(), structure))
                continue

            resn = one_letter[residue.get_resname()]

            # Vector points from center to Ca
            coordinates = residue.child_dict['CA'].get_coord()
            vector = coordinates - center
            # Take the projection on the xy plane
            vector[2] = 0
            # Normalize the vector
            normalized = vector / np.linalg.norm(vector)        
            # Give it a magnitude determined by the 'function' argument        
            complete = normalized * function(residue)
            sum_ += complete 

        moments.append((name, sum_))

        with open(path, 'wb') as f:
            target = csv.writer(f)
            for name, moment in moments:
                row = [name.upper()] + \
                      [str(component) for component in moment]
                target.writerow(row)                      
    
flexible_moments('inc_eisenmoment.csv',inclusive_selections, inc_centers, eisenberg)
flexible_moments('exc_eisenmoment.csv',exclusive_selections, exc_centers, eisenberg)    
print('done')
