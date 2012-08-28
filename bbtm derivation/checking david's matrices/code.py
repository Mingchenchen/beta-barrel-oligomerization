# Excerpts from the August 25 2012 edition of "bbtm derivation/make bbtm.py":
from __future__ import division
from __future__ import print_function
import scipy.linalg
import math
import numpy as np
import itertools
import os
import copy

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


# This function is slightly modified from the original version:
# I changed the default ordering
def transition_probability_matrix(q,t, resns = published_bbtm_ordering):
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
 
qout = 10**-4 * parse('matrices/qout.txt')

pi_out = {'A': 0.103414, 'C': 0.000253, 'E': 0.003965, 'D': 0.010899,
          'G': 0.071558, 'F': 0.088656, 'I': 0.064497, 'H': 0.016882,
          'K': 0.009315, 'M': 0.018249, 'L': 0.168981, 'N': 0.016985,
          'Q': 0.023042, 'P': 0.018898, 'S': 0.025996, 'R': 0.012083,
          'T': 0.050352, 'W': 0.045422, 'V': 0.115135, 'Y': 0.135606}

# From another file, not really much tested anyway though it's irrelevant
# where I took the code from

def expected_changes(mat, pi):
    return sum((1 - mat[i,i])*pi[i] for i in pi.keys())

def avg_rate(q, pi):
    return -1 * sum(q[i,i] *pi[i] for i in pi.keys())

# From here down is new stuff, not from make_bbtm.py

        
def parse_david(path):
    '''Open the matlab format matrix files that David Jiminez-Morales
    sent me (in the "pout" folder)'''
    with open(path, 'r') as f:
        output = MatrixMapping()
        for row_resn, line in zip(published_bbtm_ordering, f):
            for col_resn, entry in zip(published_bbtm_ordering, line.split()):
                output.update({(row_resn, col_resn): float(entry)})
    return output


def log(m):
    '''Logarithm base 10 of absolute value of every entry in a matrix'''
    o = np.empty(m.shape, dtype='float')
    for i, j in itertools.product(range(m.shape[0]), range(m.shape[1])):
        o[i, j] = math.log10(abs(m[i, j]))
    return o

p = to_mat(transition_probability_matrix(qout, 1), published_bbtm_ordering)
q = to_mat(qout, published_bbtm_ordering)
I = np.matrix(np.identity(20))
