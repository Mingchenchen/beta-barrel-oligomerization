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
import matrices
import math

class MissingPDBSequence(Exception):
    pass

class MultiplePDBSequences(Exception):
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
                if self.pdb_index is None:
                    self.pdb_index = index
                else:
                    raise MultiplePDBSequences('More than one sequence had '+\
                                               pdb_name + ' in its title')
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

    def get_alignment(self):
        return self.data

    def get_pdb_index(self):
        return self.pdb_index

    def get_pdb_seq_record(self):
        return self.data[self.pdb_index]
       

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
            if residue.get_id()[0] != 'W' and \
               residue.get_id()[1] in selection:
                print('Residue with id {0} in structure {1}'\
                      .format(residue.get_id(), structure) + \
                      'was in selection but was ignored')
            continue
        
        # Not always going to give a sequence with a moment
        if res_retrieve is not None:
            try:
                resn = res_retrieve.next()
            except StopIteration:
                raise Exception('Oracle ran out of letters on residue ' +\
                                repr(residue.get_id()))
            #resn = one_letter[residue.get_resname()]

        else:
            resn = None
        
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
            complete = normalized * mag_function(residue, resn)
        except zenergy.NoParameters:
            # To match "non normalized moments/
            #     non normalized exclusive_moment.csv":
            complete = normalized * 0
        sum_ += complete 

    return sum_

def calculator_adapter(calc, residue, resn):
    # The moment function uses a function of a residue
    # The ez_beta calculator is a function of a residue type and depth
    return calc.calculate(resn,
                          residue.child_dict['CA'].get_coord()[2])

def no_seq_calculator_adapter(calc, residue, resn):
    # Like calculator_adapter, but ignores the sequence being used
    # Instead, uses resids taken from the PDB structure
    # This is used as a double-check. If everything is functioning as it
    # should, then when calculating the moment for the pdb sequence,
    # using this should be no different than using calculator_adapter.
    return calc.calculate(residue.get_resname(),
                          residue.child_dict['CA'].get_coord()[2])

# Load selections:
selections_file = 'exc selections.csv'
with open(selections_file, 'rb') as f:
    reader = csv.reader(f)
    resi_lists = selections.selections_by_resi(reader)
print('Retrieved selections from {}... '.format(selections_file) \
      + repr(type(resi_lists)))

# Load centers:
centers_file = 'exc centers.csv'
with open(centers_file, 'rb') as f:
    reader = csv.reader(f)
    centers = load_centers(reader)
print('Retrieved centers from {}... '.format(centers_file) \
      + repr(type(centers)))    

# Initialize a calculator:
normalize = False
params_file = 'published params.csv'
with open(params_file, 'rb') as f:
    reader = csv.reader(f)
    calc = zenergy.Calculator(reader, normalize = normalize)
print('Calculator created with normalize set to {}, '
      + 'and parameters from {}... '.format(normalize, params_file) \
      + repr(type(calc)))

# Load matrix:
matrix_file = 'identity.csv'
with open(matrix_file, 'rb') as f:
    reader = csv.reader(f)
    identity = matrices.retrieve_matrix(reader)
print('Retrieved matrix from {}... '.format(matrix_file) \
      + repr(type(identity)))

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
    structure_dict = CIDict([(structure.get_id(), structure)
                             for structure in structures])
    # Yes this is stupid, next time I'll just use a dictionary, no list
print('Structured parsed... ' + repr(structures))

# Open the multiple sequence alignments:
alignments = dict([(name, Bio.AlignIO.read(filename, 'clustal')) \
              for name, filename in \
              file_dict(os.getcwd(), ['(.*) with pdb.aln']).items()])
# Remake the alignment dictionary with the pdbids that are being used
# everywhere else in this script, instead of protein names:
alignments ={'1A0S': alignments['ScrY'],
             '1QD6': alignments['OMPLA'],
             '1QJP': alignments['OmpA']}
# Make an oracle for each alignment:
oracles = CIDict([(pdbid, AlignmentOracle(alignment, pdb_name = 'chaina'))\
                  for pdbid, alignment in alignments.items()])
print('Alignments loaded... ' + repr(oracles))

# Calculate the moments for the pdb sequences
pdb_moments = CIDict([(structure.get_id(),
                   moment(structure, resi_lists[structure.get_id()],
                          centers[structure.get_id()],
                          partial(calculator_adapter, calc),
                          oracles[structure.get_id()].pdb_sequence()))
           for structure in structures])
print('pdb moments calculated! ' + repr(pdb_moments))

# Calculate the family moments, that is, the moments for all 
# sequences in the alignments
family_moments = CIDict((pdbid, list()) for pdbid in alignments.keys())

for pdbid in family_moments.keys():
    for seq_index in range(len(oracles[pdbid].get_alignment())):
        # Calculate the moment
        family_moment = moment(structure_dict[pdbid], resi_lists[pdbid],
                               centers[pdbid],
                               partial(calculator_adapter, calc),
                               oracles[pdbid].sequence(seq_index))

        # Calculate the %identity with the pdb sequence
        pdb_sequence = oracles[pdbid].get_pdb_seq_record().seq
        sequence = oracles[pdbid].get_alignment()[seq_index].seq
        normalized_distance = matrices.compare(pdb_sequence, sequence,
                                               identity)

        seq_id = oracles[pdbid].get_alignment()[seq_index].id

        family_moments[pdbid].append((seq_id, normalized_distance,
                                      family_moment))
        
        # If this is the pdb sequence, then calculating the moment using
        # only structural information, without using the MSA, should give
        # the same result. If it doesn't, something went wrong!
        if 'chaina' in seq_id:
            no_seq_moment = moment(structure_dict[pdbid], resi_lists[pdbid],
                                   centers[pdbid],
                                   partial(no_seq_calculator_adapter, calc),
                                   None)
            assert (no_seq_moment == family_moment).all(), \
                   "when calculating moment using pdb sequence of "+pdbid \
                   + ", using pdb sequence moment was " \
                   + str(family_moment)\
                   + ", but using structure, moment was " \
                   + str(no_seq_moment)

print('family moments calculated! ' + str(type(family_moments)))

for pdbid in family_moments.keys():
    with open('exc test {}.csv'.format(pdbid), 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(['name', '%id w/ pdb seq',
                         'mag', 'dir', 'x','y'])
        # Sort moments by distance from pdb sequence, descending:
        sorted_moments = sorted(family_moments[pdbid],
                                key = lambda x: x[1])[::-1]
        for seq_id, normalized_distance, family_moment in sorted_moments:
            writer.writerow([seq_id, normalized_distance,
                            np.linalg.norm(family_moment),
                            (180/math.pi) * math.atan2(family_moment[1],
                                                       family_moment[0]),
                            family_moment[0],
                            family_moment[1]])
