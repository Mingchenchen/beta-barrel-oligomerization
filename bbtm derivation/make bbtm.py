from __future__ import division
from __future__ import print_function
import scipy.linalg
import math
import numpy as np
import itertools
import os
import copy

# The guts

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

def round_to_nearest_int(n):
    '''int(n)+1 if n%1>=.5, int(n) if n%1 <.5'''
    if n%1 >= .5:
        return int(n) + 1
    else:
        return int(n)

# Constants

# Background frequencies, e-mailed to me by David Jimenez-Morales
# on August 20 2012. I consider it a birthday present
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

# Q matrices, taken from figure S3 in Jimenez-Morales, Liang PLoS One 2011
# These are the rate matrices. They have to be multipled by 10**-4 to be
# usable though.
qall = 10**-4 * parse('matrices/qall.txt')
qin = 10**-4 * parse('matrices/qin.txt')
qout = 10**-4 * parse('matrices/qout.txt')


# The API
def bbtm(dataset, t):
    '''Return a BBTM matrix from the given dataset ("all", "out" or "in")
    and at the given number of evolutionary time units (defined as the
    time it takes for one out of a hundred residues to change, on average).
    '''
    if dataset == 'all':
        pi = pi_all
        q = qall
    elif dataset == 'out':
        pi = pi_out
        q = qout
    elif dataset == 'in':
        pi = pi_in
        q = qin
    else:
        raise ValueError('the only datasets are "all", "out" and "in"')

    p = transition_probability_matrix(q, t)
    output = scoring_matrix(p, pi)
    
    return output

def write_bbtm(dataset, t, output_filename, comments = None):
    '''Write a file containing, in what I think is BLAST format, a BBTM
    matrix from the given dataset("all", "out" or "in") and at the given
    number of evolutionary time units (defined as the time it takes for
    one out of a hundred residues to change, on average)
    String in keyword argument 'comments' will also be written'''
    with open(output_filename, 'w') as o:
        # Write the first comment line, containing the dataset and t
        o.write('# bbTM{0} at t={1}\n'.format(dataset, t))
        
        # Add comments
        if comments != None:
            for line in comments.split('\n'):
                o.write('# ' + line + '\n')
        
        # Create the bbTM matrix
        matrix = bbtm(dataset, t)

        # Write the row names
        o.write(' '*4 + (' '*3).join(published_bbtm_ordering) + '\n')

        # Write the actual entries in the matrix, rounded to the nearest
        # integer
        for row in published_bbtm_ordering:
            o.write(row)
            for col in published_bbtm_ordering:
                # Create a string that is the integer entry, rounded down,
                # padded with spaces so that it takes up four characters
                entry = str(round_to_nearest_int(matrix[row,col]))
                padded_entry = ' '*(4-len(entry)) + entry
                o.write(padded_entry)
            o.write('\n')

# My dream for the following function is for it to take as an argument
# the function that actually produces the matrices, and then to write a 
# function that retrieves it from the folder David Jimenez-Morales sent me
def write_bbtm_series(dataset, id_t_triplets, output_filename, comments = None):
    '''Call write_bbtm multiple times to produce a series of matrices,
    as well as a Clustal matrix series file'''

    # Information on the series file format can be found here:
    # http://bips.u-strasbg.fr/en/Documentation/ClustalX/
    # In case the link breaks, the Wayback Machine has an archived page
    # from July 27 2011

    # Function so that I can add stuff to the output filename
    # right before a file extension, to make the filenames of the individual
    # matrices (output_filename becomes the filename of the series file)
    def make_matrix_filename(t):
        separated = output_filename.split('.')
        separated[-2] += str(t)
        return '.'.join(separated)
        
        

    # Make the matrix files
    for id_bottom, id_top, t in id_t_triplets:
        write_bbtm(dataset, t, make_matrix_filename(t), comments)

    # Make the series file
    with open(output_filename, 'w') as o:
        o.write('CLUSTAL_SERIES\n')
        o.write('\n')
        for id_bottom, id_top, t in id_t_triplets:
            o.write(' '.join([str(id_bottom), str(id_top),
                              make_matrix_filename(t)]))
            o.write('\n')

