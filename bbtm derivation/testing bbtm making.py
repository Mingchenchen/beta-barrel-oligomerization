from __future__ import division
from __future__ import print_function
import scipy.linalg
import math
import numpy as np
import itertools
import copy
from scipy.stats import gmean

published_ordering= ['S', 'T', 'N', 'Q', 'D', 'E', 'R', 'K', 'H', 'C',
                     'P', 'G', 'W', 'Y', 'F', 'V', 'I', 'L', 'A', 'M']

def parse(mat_filename):
    # Take filename of a file in the format in which matrices are presented
    # in "Patterns of Amino Acid Substitutions..." Jimenez-Morales, Jie Liang
    # and return a dictionary that can be used like q['A']['T'] to find
    # the entry in row A, column T.
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
    # The parser returns a dictionary. This function turns one of
    # those dictionaries into a matrix, with the elements in the
    # specified order. "order" should be a list of one-letter
    # residue names.
    if order is None:
        resns = m.keys()
    else:
        resns = order

    mat_as_list = list()
    for row_name in resns:
        mat_as_list.append([m[row_name][col_name] for col_name in resns])
    mat = scipy.matrix(mat_as_list)
    return mat
    
def rowsums(array):
    # Return a vector such that element i of the vector is the sum of the
    # elements in row i of a given matrix
    return np.matrix(array) * np.matrix([[1.] for i in range(array.shape[1])])

def get_params(p, b, two = 'T', one = 'S'):
    # "T" will be the equivalent of "2" in the specification, and
    # "S" will be the equivalent of "1".
    l = math.log(p[one][one]/p[two][one])/(b[one][one] - b[two][one])
    
    # Make an arbitrary ordering of one-letter residue names
    ordering = p.keys()

    pi = dict()

    for resn in ordering:
        pi.update({resn: p[one][resn]*math.exp(-1 * l * b[one][resn])})

    return (l, pi)

def get_avg_l(p, b, two = 'T', one = 'S'):
    # Like get_params, but calculates l in all possible ways and 
    # averages them, in an attempt to make up for the terrible
    # inaccuracy in the calculation
    ells = list()
    for t1, t2 in itertools.product(p.keys(), p.keys()):
        try: 
            ells.append(math.log(p[t1][t1]/p[t2][t1]) \
                        /(b[t1][t1] - b[t2][t1]))
        except ZeroDivisionError:
            continue

        # assert ells[-1] >= 0, str(ells[-1])

    l = sum(ells) / len(ells)

    
    # Make an arbitrary ordering of one-letter residue names
    ordering = p.keys()

    pi = dict()

    for resn in ordering:
        pi.update({resn: p[one][resn]*math.exp(-1 * l * b[one][resn])})

    return (l, pi)

def transition_probability_matrix(q,t):
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

def mathematica_mat(m):
    array = np.array(m)
    rows = list()
    for i in range(20):
        rows.append('{' + ', '.join(str(j) for j in array[i]) + '}')
    print('{\n' + ',\n'.join(rows) + '\n}')

b = parse('matrices/bbtmout.txt')
q = parse('matrices/qout.txt')

# Create a matrix whose rows sum to zero
q_fixed = copy.deepcopy(q)
for row in q_fixed.keys():
    rowsum = sum(q_fixed[row][col] for col in q_fixed.keys())
    q_fixed[row][row] -= rowsum
    y = q[row][row]

# Do some sanity checks on the derived parameters assuming b = exp(q*t)
# Returns lambda derived in two different ways (should be the same),
# and the sum of the residue frequencies derived in two different ways
# (should both be 1)
def test(t, q = q):
    p = transition_probability_matrix(q, t)

    l1, pi1 = get_params(p, q)
    l2, pi2 = get_params(p, q, two = 'L', one = 'A')
    print(l1)
    print(l2)
    print(sum(pi1.values()))
    print(sum(pi2.values()))

# Similar thing, but using get_avg_l instead of get_params. So,
# checking for consistency in l is meaningless.
def test2(t, q = q):
    p = transition_probability_matrix(q, t)

    l1, pi1 = get_avg_l(p, q)
    l2, pi2 = get_avg_l(p, q, two = 'L', one = 'A')
    print(l1)
    print(l2)
    print(sum(pi1.values()))
    
