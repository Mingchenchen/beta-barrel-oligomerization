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
    '''A dictionary such that if x is a MatrixMapping, (n*x)[key] == n*(x[key])
    '''
    def __rmul__(self, other):
        output = copy.deepcopy(self)
        for key in output.keys():
            output[key] *= other
        return output
    
    def __lmul__(self, other):
        return self.__rmul__(self, other)


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

def round_to_nearest_int(n):
    '''int(n)+1 if n%1>=.5, int(n) if n%1 <.5'''
    if n%1 >= .5:
        return int(n) + 1
    else:
        return int(n)

def expected_changes(mat, pi):
    return sum((1 - mat[i,i])*pi[i] for i in pi.keys())
