from __future__ import division
from __future__ import print_function
import scipy.linalg
import math
import numpy as np
import itertools
import os

class MatrixMapping(dict):
    def __rmul__(self, other):
        output = dict()
        for i in self.keys():
            inner_dict = dict()
            for j in self.keys():
                inner_dict.update({j: self[i][j] * other})
            output.update({i: inner_dict})
        return MatrixMapping(output)

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
        
    return MatrixMapping(output)

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
    resns = published_q_ordering
    mat = to_mat(q, order = resns)
    p_as_mat = scipy.linalg.expm(t*mat)
    output = dict((resn, None) for resn in resns)
    for key in output.keys():
        output[key] = dict((resn, None) for resn in resns)

    numrows, numcols = p_as_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i]][resns[j]] = p_as_mat[i,j]

    return MatrixMapping(output)


def scoring_matrix(p, lam, pi):
    '''p(A, G) becomes (1/lam) * log(p(A,G)/pi(G))'''
    b = dict()
    for x in p.keys():
        b.update({x: dict((x, None) for y in p.keys())})
    for from_, to in itertools.product(p.keys(), p.keys()):
        b[from_][to] = (1/lam) * math.log(p[from_][to] / pi[to])
    return MatrixMapping(b)

def expected_changes(mat, pi):
    return sum((1 - mat[i][i])*pi[i] for i in pi.keys())

def avg_rate(q, pi):
    return -1 * sum(q[i][i] *pi[i] for i in pi.keys())


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

os.chdir(r'C:\cygwin\home\alex\beta barrels\bbtm derivation\matrices')

qall = parse('qall.txt')
qin = parse('qin.txt')
qout = parse('qout.txt')

# Normalize by from
#for i, j in itertools.product(pi_out.keys(), pi_out.keys()):
#    qout[i][j] /= pi_out[i]
    
mult = 10**-4
pout = transition_probability_matrix(qout, 40 * mult)
bout = scoring_matrix(pout, 1/2.798, pi_out)
to_mat(bout, published_bbtm_ordering)
print(expected_changes(pout, pi_out))
print(avg_rate(mult * qout, pi_out))

bbtmout = parse('bbtmout.txt')

wrongones = dict()
for i, j in itertools.product(pi_out.keys(), pi_out.keys()):
    display = bbtmout[i][j] - bout[i][j]
    if abs(display) > 3:
        print(i + j + ': ' + str(display))
        if j not in wrongones.keys():
            wrongones.update({j: 1})
        else:
            wrongones[j] += 1
print(wrongones)

'''
for i, j in itertools.product(pi_out.keys(), pi_out.keys()):
    print(i + j + ' ' + str(bbtmout[i][j] / bout[j][i]) + ', ' + str(bout[j][i]))
'''

print(sum((bbtmout[i][j] / bout[i][j]) for i, j in \
          itertools.product(pi_out.keys(), pi_out.keys()))/400)

#print(np.array(\
#               to_mat(bout, published_bbtm_ordering),
#               dtype = 'int32'))
