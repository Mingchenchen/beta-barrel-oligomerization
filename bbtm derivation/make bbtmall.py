import scipy.linalg
import math
import numpy as np
import itertools

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
    
def to_mat(m, order = None):
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


b = parse('matrices/bbtmout.txt')
q = parse('matrices/qout.txt')
p = transition_probability_matrix(q, 10**-4 * 40)

l1, pi1 = get_params(p, q)
l2, pi2 = get_params(p, q, two = 'L', one = 'A')
