import numpy as np
import scipy
<<<<<<< HEAD
import matrices
=======
import itertools
import copy
>>>>>>> b765eb524c9407563c8aa45467e5e47635c784a5

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

published_bbtm_ordering = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
                           'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                           'Y', 'V']

# Background frequencies, e-mailed to me by David Jimenez-Morales
# on August 20 2012. I consider it a birthday present

pi_out = {'A': 0.103414, 'C': 0.000253, 'E': 0.003965, 'D': 0.010899,
          'G': 0.071558, 'F': 0.088656, 'I': 0.064497, 'H': 0.016882,
          'K': 0.009315, 'M': 0.018249, 'L': 0.168981, 'N': 0.016985,
          'Q': 0.023042, 'P': 0.018898, 'S': 0.025996, 'R': 0.012083,
          'T': 0.050352, 'W': 0.045422, 'V': 0.115135, 'Y': 0.135606}

# Amino acid frequencies from http://www.tiem.utk.edu/~gross/bioed/webmodules/aminoacid.htm retrieved August 27 2012
pi_ver = dict({'A': 7.4e-2,
               'R': 4.2e-2,
               'N': 4.4e-2,
               'D': 5.9e-2,
               'C': 3.3e-2,
               'E': 5.8e-2,
               'Q': 3.7e-2,
               'G': 7.4e-2,
               'H': 2.9e-2,
               'I': 3.8e-2,
               'L': 7.6e-2,
               'K': 7.2e-2,
               'M': 1.8e-2,
               'F': 4.0e-2,
               'P': 5.0e-2,
               'S': 8.1e-2,
               'T': 6.2e-2,
               'W': 1.3e-2,
               'Y': 3.3e-2,
               'V': 6.8e-2})
<<<<<<< HEAD
pam1 = 10**-4 * matrices.parse('pam1.txt')

=======

def expected_changes(mat, pi):
    '''Expected fraction of positions different from original, after update
    with mat, given original frequency distribution pi.
    Very different from expected number of TRANSITIONS.'''
    return sum((1 - mat[i,i])*pi[i] for i in pi.keys())

def avg_rate(q, pi):
    return -1 * sum(q[i,i] *pi[i] for i in pi.keys())

def parse(mat_filename):
    '''Take filename of a file in the format in which matrices are presented
    in "Patterns of Amino Acid Substitutions..." Jimenez-Morales, Jie Liang
    and return a dictionary that can be used like q['A','T'] to find
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

pam1 = 10**-4 * parse('pam1.txt')

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
        
>>>>>>> b765eb524c9407563c8aa45467e5e47635c784a5
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

<<<<<<< HEAD
def p_retrieve(t):
    '''Return one of David Jimenez-Morales's transition probability
    matrices, corresponding to the given time t, from the matrices
    he sent me on August 23 2012.'''
    return parse_david('TMout/pout/MTMout{0}.p'.format(t))

def id_given_time(t, mat_name='bbtm'):
    '''Given a time and the name of a matrix family (bbtm and pam
    currently available) return the expected value of %identity between
    a sequence and the same sequence after a time t has passed,
    given that the amino acid frequencies in the original sequence
    are equal to the background frequencies'''
    if mat_name == 'pam':
        pam_at_t = matrices.matrix_power(pam1, t)
        return 1 - matrices.expected_changes(pam_at_t, pi_ver)
    
    if mat == 'bbtm':
        bbtm_at_t = p_retrieve(t)
        return 1 - matrices.expected_changes(bbtm_at_t, pi_out)
        
=======
def david_changes(n):
    return expected_changes(parse_david('TMout/pout/MTMout{0}.p'\
                                        .format(n)), pi_out)

>>>>>>> b765eb524c9407563c8aa45467e5e47635c784a5
