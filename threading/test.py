from __future__ import division

from sundries import CIDict
import csv
from Bio import AlignIO
from Bio import SeqIO

# Uses CIDict
def retrieve_matrix(csvfile):
    # Create a dictionary mapping column numbers to resns
    colnum_to_resn = dict()
    # The first line is going to have an empty space as its first entry,
    # so start at the second entry
    for index, entry in enumerate(csvfile.next()[1:]):
        colnum_to_resn.update({index: entry})
    
    # Create a dictionary so dict['A']['V'] will give you score for
    # alanine mutating into valine
    first_to_second = CIDict([(resn, CIDict())
                            for resn in colnum_to_resn.values()])
    for line in csvfile:
        second_resn = line[0]
        for colnum, number in enumerate(line[1:]):
            first_resn = colnum_to_resn[colnum]
            # Add number as floating point so I can sum these things up
            # to get distances
            first_to_second[first_resn].update({second_resn: float(number)})
    
    return first_to_second

# Doesn't need any modules
def distance(seq1, seq2, matrix):
    '''Get a distance score using a given matrix for a pair of aligned sequences
    Or it might be a similarity score, depending on the matrix.
    '''
    if len(seq1) != len(seq2):
        raise ValueError('Sequences of different length, must compare ' + \
                         'aligned sequences')

    distance = 0
    for pos in range(len(seq1)):
        # These matrices are not necessarily symmetric, order matters
        distance += matrix[seq1[pos]][seq2[pos]]

    return distance

# Uses future division
def compare(seq1, seq2, matrix):
    '''Comparison function I made just for this script, divides a distance
    by number of positions that are not gaps in both of the sequences compared
    to return a normalized distance'''
    count = 0
    for pos in range(len(seq1)):
        if seq1[pos] != '-' or seq2[pos] != '-':
            count += 1

    return distance(seq1, seq2, matrix) / count




matrices = []    
for matrix in ('hamming.csv', 'bbTMall40.csv'):
    with open('matrices/' + matrix, 'rb') as f:
        freader = csv.reader(f)
        matrices.append(retrieve_matrix(freader))

hamming = matrices[0]
bbtm = matrices[1]

x = 'LALA'
y = 'LARG'
# Hamming distance is 2, bbTM distance is 6+6+(-11)+1=2


ompla_msa = AlignIO.read('test data/1qd6.clu', 'clustal')
ompla_pdb = ompla_msa[5]


sorted_seqs = sorted(ompla_msa[4:],
                     key = lambda seq: compare(ompla_pdb, seq, hamming))[::-1]

