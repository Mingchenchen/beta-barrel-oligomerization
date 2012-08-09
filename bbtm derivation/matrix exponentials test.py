import scipy.linalg
import itertools
import math
import os
import numpy as np

os.chdir(r'C:\cygwin\home\alex\beta barrels\bbtm derivation')


def parse(mat_file):
    # Take in an open file in the format in which matrices are presented
    # in "Patterns of Amino Acid Substitutions..." Jimenez-Morales, Jie Liang
    # and return a dictionary that can be used like q['A']['T'] to find
    # the entry in row A, column T.
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

def to_mat(m, order = 'None'):
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

def p(q,t):
    mat = to_mat(q)
    p_as_mat = scipy.linalg.expm(t*mat)
    output = dict((resn, None) for resn in resns)
    for key in output.keys():
        output[key] = dict((resn, None) for resn in resns)

    numrows, numcols = p_as_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i]][resns[j]] = p_as_mat[i,j]

    return output

def rowsums(array):
    # Return a vector such that element i of the vector is the sum of the
    # elements in row i of a given matrix
    return np.matrix(array) * np.matrix([[1.] for i in range(array.shape[1])])

def test(power, subject):
    # Return a vector of the sum of each row of e^(subject*10^power)
    # Used to test how these change as the factor by which the subject
    # is multiplied increases
    return rowsums(scipy.linalg.expm(subject*10**power))

def test2(power, subject):
    # Same as test using scipy.linalg.expm2 instead of scipy.linalg.expm
    # Same function, different algorithm
    return rowsums(scipy.linalg.expm2(subject*10**power))

def test3(power, subject, q=20):
    return rowsums(scipy.linalg.expm3(subject*10**power, q=q))

with open('matrices/qout.txt', 'r') as f:
    r_dict = parse(f)

published_ordering= ['S', 'T', 'N', 'Q', 'D', 'E', 'R', 'K', 'H', 'C',
                     'P', 'G', 'W', 'Y', 'F', 'V', 'I', 'L', 'A', 'M']
r = to_mat(r_dict, order = published_ordering)

# Make a new matrix whose rows actually sum to zero like they do in the
# unrounded matrix
r_fixed = np.copy(r)
for rownum in range(20):
    rowsum = 0
    for colnum in range(20):
        rowsum += r_fixed[rownum, colnum]
    r_fixed[rownum, rownum] -= rowsum


# Make a list of the elements of e^(r * 40)
elem_list = [scipy.linalg.expm(r*40)[i, j] \
             for i, j in itertools.product(range(20), range(20))]
print(max([i - .05 for i in elem_list]))

# Same thing for r_fixed:
elem_list = [scipy.linalg.expm(r_fixed*40)[i, j] \
             for i, j in itertools.product(range(20), range(20))]
print(max([i - .05 for i in elem_list]))
