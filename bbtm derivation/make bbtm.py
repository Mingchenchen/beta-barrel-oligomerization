from __future__ import division
from __future__ import print_function
import scipy.linalg
import math
import numpy as np
import itertools

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
                output = dict([(resn, None) for resn in col_names])
                for key, value in output.items():
                    output[key] = dict([(resn, None) for resn in col_names])
                continue

            row = line.split()
            row_name = line[0]
            for rate, col_name in zip(row[1:], col_names):
                output[row_name][col_name] = float(rate)
        
    return output

published_ordering= ['S', 'T', 'N', 'Q', 'D', 'E', 'R', 'K', 'H', 'C',
                     'P', 'G', 'W', 'Y', 'F', 'V', 'I', 'L', 'A', 'M']    
def to_mat(m, order = published_ordering):
    '''The parser returns a dictionary. This function turns one of
    those dictionaries into a matrix, with the elements in the
    specified order. "order" should be a list of one-letter
    residue names.'''
    if order is None:
        resns = m.keys()
    else:
        resns = order

    mat_as_list = list()
    for row_name in resns:
        mat_as_list.append([m[row_name][col_name] for col_name in resns])
    mat = scipy.matrix(mat_as_list)
    return mat
   
def transition_probability_matrix(q,t):
    '''Take a matrix in dictionary form, whose elements are instantaneous
    transition rates, and return a transition probability matrix for a given
    time.'''
    resns = q.keys()
    mat = to_mat(q)
    p_as_mat = scipy.linalg.expm(t*mat)
    output = dict((resn, None) for resn in resns)
    for key in output.keys():
        output[key] = dict((resn, None) for resn in resns)

    numrows, numcols = p_as_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i]][resns[j]] = p_as_mat[i,j]

    return output


def scoring_matrix(p, lam, pi):
    '''p(A, G) becomes (1/lam) * log(p(A,G)/pi(G))'''
    b = dict()
    for x in p.keys():
        b.update({x: dict((x, None) for y in p.keys())})
    for from_, to in itertools.product(p.keys(), p.keys()):
        b[from_][to] = (1/lam) * math.log(p[from_][to] / pi[to])
    return b
