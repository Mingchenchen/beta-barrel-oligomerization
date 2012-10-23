from __future__ import division
from __future__ import print_function

from sundries import CIDict
import scipy.linalg
import math
import numpy as np
import itertools
import os
import copy

# These first three functions are antiquated, kept for compatibility
# with my older code
def retrieve_matrix(csvfile):
    '''
    retrieve_matrix(csvfile)
    Returns a CIDict

    Given a csv file whose rows and columns are labeled, returns a 
    case insensitive dictionary of case insensitive dictionaries
    that will retrieve numberical elements in the csvfile as
    floating point numbers, using the syntax
    dictionary['columnlabel']['rowlable']

    The intention is to do dictionary['original']['mutated'] to get
    a score for a mutation, if the originals are in the rows then
    I did it wrong
    '''
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

def distance(seq1, seq2, matrix):
    '''
    distance(seq1, seq2, matrix)
    Returns a number (floating point if using a matrix generated
                      by retrieve_matrix)
    Get a distance score using a given matrix for a pair of aligned sequences
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

def compare(seq1, seq2, matrix):
    '''
    compare(seq1, seq2, matrix)
    Returns a number

    Just uses 'distance' function with same arguments, and then divides
    by the number of positions in the two sequences that are not both
    gaps
    '''
    count = 0
    for pos in range(len(seq1)):
        if seq1[pos] != '-' or seq2[pos] != '-':
            count += 1

    return distance(seq1, seq2, matrix) / count

# Here's the functions that use the modern way I deal with matrices
class MatrixMapping(dict):
    '''A representation of a matrix as a dictionary. Supports scalar
    multiplication and matrix exponentiation. If x is a MatrixMapping,
    try 5*x, x*5, and x**5.
    '''
    # The mathematical meaningfulness of matrix exponentiation with this
    # object seems to depend upon self.keys() being a Caresian
    # product (set of row names is same as set of column names).
    # I haven't done the linear algebra to work that out yet, though
    def __rmul__(self, other):
        output = copy.deepcopy(self)
        for key in output.keys():
            output[key] *= other
        return output
    
    def __lmul__(self, other):
        return self.__rmul__(self, other)

    def __pow__(self, other):
        '''Matrix exponentiation.'''
        # This is just a matrix_power call. However, since the order of the
        # "resns" argument of matrix_power is supposedly arbitrary,
        # why not just use self.keys()? Then, without that keyword argument,
        # there's no reason matrix_power should exist as a separate
        # function.
        # Seems clear, but needs testing before I make changes.
        if set(i[0] for i in self.keys()) != set(published_q_ordering):
            raise Exception('Exponentiation is only supported for matrices'\
                            + ' of the 20 standard amino acids at the'\
                            + ' moment, and the keys must be the'\
                            + ' one-letter codes.')

        return matrix_power(self, other)

def parse(mat_filename):
    '''Take filename of a file in the format in which matrices are presented
    in "Patterns of Amino Acid Substitutions..." Jimenez-Morales, Jie Liang
    and return a dictionary that can be used like q['A']['T'] to find
    the entry in row A, column T.'''
    with open(mat_filename, 'r') as mat_file:
        found_resns = False
        for line in mat_file:
            # Ignore comments
            if line[0] == '#':
                continue

            # Find which resn corresponds to which column number
            if not found_resns:
                found_resns = True
                col_names = line.split()
                # Check to make sure they're all there
                for resn in ['C', 'N', 'H', 'D', 'S', 'Q', 'K', 'M', 'P',
                             'T', 'F', 'A', 'G', 'I', 'L', 'R', 'W', 'E',
                             'Y', 'V']:
                    try:
                        assert resn in col_names, "missing " + resn
                    except AssertionError:
                        print(line)
                        raise

                # Make the matrix that will be returned
                output = MatrixMapping((tuple_, None)\
                                       for tuple_ in itertools.product(\
                                                       col_names, col_names))
                continue

            row = line.split()
            row_name = line[0]
            for rate, col_name in zip(row[1:], col_names):
                output[row_name,col_name] = float(rate)
        
    return output

published_q_ordering = ['S', 'T', 'N', 'Q', 'D', 'E', 'R', 'K', 'H', 'C',
                     'P', 'G', 'W', 'Y', 'F', 'V', 'I', 'L', 'A', 'M']

published_bbtm_ordering = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
                           'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                           'Y', 'V']
def to_mat(m, order):
    '''The parser returns a dictionary. This function turns one of
    those dictionaries into a matrix, with the elements in the
    specified order. "order" should be a list of one-letter
    residue names.'''

    # I know more about manipulating lists than matrices,
    # so the output is constructed as a list, then converted into a
    # matrix right before the return statement

    mat_as_list = list()
    for row_name in order:
        # Append a row:
        mat_as_list.append([m[row_name,col_name] for col_name in order])
    mat = scipy.matrix(mat_as_list)
    return mat
  

def transition_probability_matrix(q,t, resns = published_q_ordering):
    '''Take a matrix in dictionary form, whose elements are instantaneous
    transition rates, and return a transition probability matrix for a given
    time.'''

    # resns is a list of entry names (probably the 20 amino acids).
    # I think the order is arbitrary.
    # Made into a keyword argument so that I can test this function on
    # small matrices.

    # Convert q to a matrix
    mat = to_mat(q, order = resns)

    # Create P(t) = exp(q * t)
    p_as_mat = scipy.linalg.expm(t*mat)

    # Make a dictionary mapping tuples of resns to elements of P(t)
    # Create the empty dictionary
    output = MatrixMapping((respair, None) \
                  for respair in itertools.product(resns, resns))

    # Fill it up with the values from the matrix
    numrows, numcols = p_as_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i],resns[j]] = p_as_mat[i,j]

    return output

def matrix_power(q,t, resns = published_q_ordering):
    # resns is a list of entry names (probably the 20 amino acids).
    # I think the order is arbitrary.
    # Made into a keyword argument so that I can test this function on
    # small matrices.

    # Convert q to a matrix
    mat = to_mat(q, order = resns)

    # Take the power
    exponentiated_mat = np.linalg.matrix_power(mat, t)

    # Make a dictionary mapping tuples of resns to elements of P(t)
    # Create the empty dictionary
    output = MatrixMapping((respair, None) \
                  for respair in itertools.product(resns, resns))

    # Fill it up with the values from the matrix
    numrows, numcols = exponentiated_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i],resns[j]] = exponentiated_mat[i,j]

    return output


def scoring_matrix(p, pi, lam = 1/math.log(16)):
    '''p(A, G) becomes (1/lam) * log(p(A,G)/pi(G))'''

    b = MatrixMapping()
    for from_, to in p.keys():
        try:
            b.update({(from_, to): (1/lam) * math.log(p[from_,to] / pi[to])})
        except ValueError:
            raise ValueError('p[from,to] is {0}, pi[to] is {1}, '/\
                             format(p[from_,to], pi[to]) \
                             + 'can\'t take log of ratio')

    return b

def expected_changes(mat, pi):
    '''Expected number of changes after update with transition probability
    matrix mat on frequency distribution pi'''
    return sum((1 - mat[i,i])*pi[i] for i in pi.keys())

# Frequency distributions for use with expected_changes
# Frequencies of amino acids in the dataset used to derive the PAM
# matrices, from Dayhoff et. al 1978
pam_freq = {'G': .089, 'R': .041,
            'A': .087, 'N': .040,
            'L': .085, 'F': .040,
            'K': .081, 'Q': .038,
            'S': .070, 'I': .037,
            'V': .065, 'H': .034,
            'T': .058, 'C': .033,
            'P': .051, 'Y': .030,
            'E': .050, 'M': .015,
            'D': .047, 'W': .010}

# Frequency distributions for the datasets in Jimenez-Morales et. all
# 2011, unpublished:
pi_all = {'A': 0.090686, 'C': 0.000344, 'E': 0.037944, 'D': 0.029356,
          'G': 0.109694, 'F': 0.06036, 'I': 0.041904, 'H': 0.014036,
          'K': 0.02891, 'M': 0.018884, 'L': 0.104019, 'N': 0.037032,
          'Q': 0.043966, 'P': 0.012709, 'S': 0.066126, 'R': 0.04019,
          'T': 0.067868, 'W': 0.031041, 'V': 0.071151, 'Y': 0.09384}
pi_in = {'A': 0.078272, 'C': 0.000478, 'E': 0.072977, 'D': 0.047942,
         'G': 0.1498, 'F': 0.028562, 'I': 0.018918, 'H': 0.011209,
         'K': 0.049208, 'M': 0.019776, 'L': 0.037144, 'N': 0.057714,
         'Q': 0.065387, 'P': 0.006456, 'S': 0.107541, 'R': 0.069281,
         'T': 0.086192, 'W': 0.015694, 'V': 0.026265, 'Y': 0.05133}
pi_out = {'A': 0.103414, 'C': 0.000253, 'E': 0.003965, 'D': 0.010899,
          'G': 0.071558, 'F': 0.088656, 'I': 0.064497, 'H': 0.016882,
          'K': 0.009315, 'M': 0.018249, 'L': 0.168981, 'N': 0.016985,
          'Q': 0.023042, 'P': 0.018898, 'S': 0.025996, 'R': 0.012083,
          'T': 0.050352, 'W': 0.045422, 'V': 0.115135, 'Y': 0.135606}


def parse_david(path):
    '''Open the matlab format matrix files that David Jiminez-Morales
    sent me (in the "pout" folder)'''
    with open(path, 'r') as f:
        output = MatrixMapping()
        for row_resn, line in zip(published_bbtm_ordering, f):
            for col_resn, entry in zip(published_bbtm_ordering,
                                       line.split()):
                output.update({(row_resn, col_resn): float(entry)})
    return output

