import numpy as np
import scipy
import matrices

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
pam1 = 10**-4 * matrices.parse('pam1.txt')

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
        
