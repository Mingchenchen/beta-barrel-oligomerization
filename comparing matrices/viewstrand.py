from __future__ import division
from __future__ import print_function

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

abridged_dicts = list()
for filename in (r'C:\cygwin\home\alex\temp\scry with bbtm.aln',
                 r'C:\cygwin\home\alex\beta barrels\test moments\ScrY with pdb.aln'):
    alignment = Bio.AlignIO.read(filename, 'clustal')
    oracle = AlignmentOracle(alignment, pdb_name = 'chaina')

    abridged = dict()
    for i in range(len(alignment)):
        abridged.update({oracle.data[i].id: (''.join(oracle.sequence(i)))})

    abridged_dicts.append(abridged)

bbtm = abridged_dicts[0]
gonnet = abridged_dicts[1]

strand_limits = (150, 181)
print('name then bbtm segment then gonnet segment')
for name in bbtm.keys():
    print(name)
    print(bbtm[name][strand_limits[0]: strand_limits[1]])
    print(gonnet[name][strand_limits[0]:strand_limits[1]])
    print()

