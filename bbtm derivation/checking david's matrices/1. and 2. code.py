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

def scoring_matrix(p, pi, lam = 1/math.log(16)):
    '''p(A, G) becomes (1/lam) * log(p(A,G)/pi(G))'''

    b = MatrixMapping()
    for from_, to in p.keys():
        try:
            b.update({(from_, to): (1/lam) * math.log(p[from_,to] / pi[to])})
        except ValueError:
            raise ValueError('p[from,to] is {}, pi[to] is {}, can\'t take\
                              log of ratio'.format(p[from_,to], pi[to]))

    return b

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


 
qout = 10**-4 * parse('matrices/qout.txt')

pi_out = {'A': 0.103414, 'C': 0.000253, 'E': 0.003965, 'D': 0.010899,
          'G': 0.071558, 'F': 0.088656, 'I': 0.064497, 'H': 0.016882,
          'K': 0.009315, 'M': 0.018249, 'L': 0.168981, 'N': 0.016985,
          'Q': 0.023042, 'P': 0.018898, 'S': 0.025996, 'R': 0.012083,
          'T': 0.050352, 'W': 0.045422, 'V': 0.115135, 'Y': 0.135606}

# From here down is new stuff, not from make_bbtm.py

def round_to_nearest_int(n):
    '''Round to nearest integer'''
    if n >= 0:
        if abs(int(n) + 1 - n) < abs(int(n) - n):
            return int(n) + 1
        else:
            return int(n)
    elif n < 0:
        # int(n) rounds towards zero, so rounding "up" is subtracting 1
        if abs(int(n) - 1 - n) < abs(int(n) - n):
            return int(n) - 1
        else:
            return int(n)

        
def parse_david(path):
    '''Open the matlab format matrix files that David Jiminez-Morales
    sent me (in the "pout" folder)'''
    with open(path, 'r') as f:
        output = MatrixMapping()
        for row_resn, line in zip(published_bbtm_ordering, f):
            for col_resn, entry in zip(published_bbtm_ordering, line.split()):
                output.update({(row_resn, col_resn): float(entry)})
    return output

def round_matrix(m):
    '''Return a matrix in which each entry of m is rounded to the nearest
    integer'''
    # o is for output
    o = np.empty(m.shape, dtype='int')

    # Round each entry
    for i, j in itertools.product(range(m.shape[0]), range(m.shape[1])):
        o[i,j] = round_to_nearest_int(m[i,j])

    return o
        

mtmout20 = parse_david('MTMout20.p')
computed_bbtmout20_as_dict = scoring_matrix(mtmout20, pi_out)
computed_bbtmout20_as_floatmat = to_mat(computed_bbtmout20_as_dict,
                                        published_bbtm_ordering)
computed_bbtmout20 = round_matrix(computed_bbtmout20_as_floatmat)

published_bbtmout40_as_dict = parse('bbtmout.txt')
published_bbtmout40 = to_mat(published_bbtmout40_as_dict,
                             published_bbtm_ordering)

difference = published_bbtmout40 - computed_bbtmout20


