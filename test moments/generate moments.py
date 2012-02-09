import numpy as np
import zenergy
import csv
import selections

def moment(structure, selection, center, function): 
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

    return sum_

with open('cored 1 selections with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    resi_lists = selections.selections_by_resi(reader)

with open('cored 1 centers with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    centers = zenergy.load_centers(reader)

with open('published params.csv', 'rb') as f:
    reader = csv.reader(f)
    calc = zenergy.Calculator(reader, normalize = True)
