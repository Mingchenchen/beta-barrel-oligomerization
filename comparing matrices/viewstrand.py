from __future__ import division
from __future__ import print_function

import sys
sys.path.append("/home/davis/work stuff/beta-barrel-oligomerization/modules")

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

# Create dictionaries that map names of proteins to other dictionaries
bbtm_alignments = dict()
gonnet_alignments = dict()

for protein in ('OMPLA',):
    # Create dictionaries that map names of sequences to sequences containing
    # only the positions that are aligned to non-gap positions in the pdb
    # sequence
    abridged_dicts = list()
    for filename in (protein + ' bbtm.aln',
                     protein + ' gonnet.aln'):
        alignment = Bio.AlignIO.read(filename, 'clustal')
        oracle = AlignmentOracle(alignment, pdb_name = 'chaina')

        abridged = dict()
        for i in range(len(alignment)):
            abridged.update({oracle.data[i].id: (''.join(oracle.sequence(i)))})

        abridged_dicts.append(abridged)

    bbtm_alignments.update({protein: abridged_dicts[0]})
    gonnet_alignments.update({protein: abridged_dicts[1]})

with open('identity.csv', 'rb') as f:
    id_matrix = matrices.retrieve_matrix(csv.reader(f))

def check_range(start, end, protein):
    # Observe what positions of each protein in the gonnet alignment and the
    # bbtm alignment correspond to a given range in the pdb sequence.
    # Positions ocrresponding to gaps in the pdb sequence not shown.
    print('name then seq distance then bbtm segment then gonnet segment')
    pdb_seq_name = 'chaina_' + protein.lower()
    for name in bbtm_alignments[protein].keys():
        print(name)
        print(matrices.compare(bbtm_alignments[protein][name],
                               bbtm_alignments[protein][pdb_seq_name],
                               id_matrix))
        print(matrices.compare(gonnet_alignments[protein][name],
                               gonnet_alignments[protein][pdb_seq_name],
                               id_matrix))
        print(bbtm_alignments[protein][name][start: end])
        print(gonnet_alignments[protein][name][start:end])
        print()
    
check_range(56,68,'OMPLA')
