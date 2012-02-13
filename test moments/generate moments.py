import numpy as np
import zenergy
import csv
import selections
from functools import partial
import Bio.AlignIO
from sundries import CIDict
from sundries import file_dict
from Bio.PDB import PDBParser
import warnings
import os

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

def load_centers(iterable):  
    dict_ = CIDict()
    for row in iterable:
        if row[0] != '':
            dict_.update({row[0]:row[1]})
    for key, value in dict_.items():
        # Turns '(1,2,3)' etc, that is, textual representations of vectors,
        # into Vector objects. Will cut off the last digit of the third
        # componant, but I don't care because the third componant will
        # always be 0.0
        dict_[key] = np.array([float(y[:-1]) for y in value[1:].split()])
    return dict_

def moment(structure, selection, center, mag_function, res_retrieve): 
    sum_ = np.zeros(3)
    for residue in structure.get_residues():

        # My selection files just give residue numbers. If two residues
        # have the same number, if there're insertions, then the selection
        # files are ambiguous and I need to change how I do them.
        # For now, anything identified by more than just a residue number
        # is ignored. The reason this is okay is that, as far as I can tell,
        # in Dan's structures that just means it ignores water.
        # If the structures have complex residue ids, then I need different
        # selection files.
        if residue.get_id()[0] != ' ':
            if residue.get_id()[0] != 'W' and residue.get_id()[1] in selection:
                print('Residue with id {0} in structure {1} was in selection'\
                       .format(residue.get_id(), structure) + \
                       ' but was ignored')
            continue
        
        try:
            resn = res_retrieve.next()
        except StopIteration:
            raise Exception('Oracle ran out of letters on residue ' +\
                            repr(residue.get_id()))
        #resn = one_letter[residue.get_resname()]
        
        if residue.get_id()[1] not in selection:
            continue




        # Vector points from center to Ca
        coordinates = residue.child_dict['CA'].get_coord()
        vector = coordinates - center
        # Take the projection on the xy plane
        vector[2] = 0
        # Normalize the vector
        normalized = vector / np.linalg.norm(vector)        
        # Give it a magnitude determined by the 'function' argument        
        try:
            complete = normalized * mag_function(residue)
        except zenergy.NoParameters:
            # For now, for the purposes of replicating old data:
            complete = normalized * .5
        sum_ += complete 

    return sum_

def calculator_adapter(calc, residue):
    # The moment function uses a function of a residue
    # The ez_beta calculator is a function of a residue type and depth
    return calc.calculate(residue.get_resname(),
                          residue.child_dict['CA'].get_coord()[2])

# Load selections:
with open('cored 1 selections with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    resi_lists = selections.selections_by_resi(reader)
print('Selections retrieved... ' + repr(type(resi_lists)))

# Load centers:
with open('cored 1 centers with 1qd5.csv', 'rb') as f:
    reader = csv.reader(f)
    centers = load_centers(reader)
print('Centers retrieved... ' + repr(type(centers)))    

# Initialize a calculator:
with open('published params.csv', 'rb') as f:
    reader = csv.reader(f)
    calc = zenergy.Calculator(reader, normalize = True)
print('Calculator created... ' + repr(type(calc)))

# The slow part, open the structures:
structure_files = file_dict('structures with 1qd5',
                            ['aligned_(1QD6).pdb',
                             'aligned_(1QJP).pdb',
                             'aligned_(1A0S).pdb'])

parser = PDBParser()

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    structures = [parser.get_structure(id_, path) \
                  for id_, path in structure_files.items()]
print('Structured parsed...' + repr(structures))

# Open the multiple sequence alignments:
alignments = dict([(name, Bio.AlignIO.read(filename, 'clustal')) \
              for name, filename in \
              file_dict(os.getcwd(), ['(.*) with pdb.aln']).items()])
# Remake the alignment dictionary with the pdbids that are being used
# everywhere else in this script, instead of protein names:
pdbnames ={'1A0S': alignments['ScrY'],
           '1QD6': alignments['OMPLA'],
           '1QJP': alignments['OmpA']}
# Make an oracle for each alignment:
oracles = CIDict([(pdbid, AlignmentOracle(alignment, pdb_name = 'chaina'))\
                  for pdbid, alignment in pdbnames.items()])
print('Alignments loaded... ' + repr(oracles))

moments = CIDict([(structure.get_id(),
                   moment(structure, resi_lists[structure.get_id()],
                          centers[structure.get_id()],
                          partial(calculator_adapter, calc),
                          oracles[structure.get_id()].pdb_sequence()))
           for structure in structures])
print('Moments calculated! ' + repr(moments))
