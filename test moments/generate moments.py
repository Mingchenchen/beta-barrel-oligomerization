import numpy as np
import zenergy
import csv
import selections

class MissingPDBSequence(Exception):
    pass

class AlignmentOracle(object):
    '''
    Initialized with a multiple sequence alignment and a name to find a
    particular sequence of interest. Can give amino acids from other
    sequences in the alignment that correspond to amino acids in the
    sequence
    of interest, ignoring amino acids that are lined up with gaps in the
    sequence of interest.
    '''.replace('  ', '')
    def __init__(self, alignment, pdb_name = 'pdb'):
        self.data = alignment
        self.pdb_index = None
        for index, record in enumerate(self.data):
            if pdb_name in record.id:
                self.pdb_index = index
                break
        if self.pdb_index == None:
            raise MissingPDBSequence('None of the sequences had '+\
            pdb_name + ' in their title')
        self._pdb_sequence = self.data[self.pdb_index].seq

    def pdb_sequence(self, selection = None):
        return self.sequence(self.pdb_index, selection)

    def sequence(self, index, selection = None):
        pos = 0
        for j, letter in enumerate(self.data[index].seq):
            if self._pdb_sequence[j] != '-':
                if selection is None or pos in selection:
                    yield letter.upper()
                pos += 1


def moment(structure, selection, center, function, res_retrieve): 
    sum_ = np.zeros(3)
    for residue in structure.get_residues():
        if residue.get_id()[1] not in selection:
            # Should keep count of how far, for alignment oracle
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

        # This is what's gonna use the alignment oracle:
        resn = one_letter[residue.get_resname()]

        # Also need to get the z-coordinate, to use the ez-beta calculator
        # as a function

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

structure_files = file_dict('structures with 1qd5',
                            ['aligned_(.*).pdb'])

parser = PDBParser()

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    structures = [parser.get_structure(id_, path) \
                  for id_, path in structure_files.items()]

moments = [moment(structure, resi_lists[structure.get_id()],
                  centers[structure.get_id()], calc.calculate, lambda x: x) \ 
           for structure in structures]
