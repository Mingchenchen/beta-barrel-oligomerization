import scipy.linalg
import itertools
import math

def parse(mat_file):
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
                assert resn in col_names, "missing " + resn

            # Make the matrix that will be returned
            output = dict([(resn, None) for resn in col_names])
            for key, value in output.items():
                output[key] = dict([(resn, None) for resn in col_names])
            continue

        row = line.split()
        row_name = line[0]
        for rate, col_name in zip(row[1:], col_names):
            output[row_name][col_name] = rate
    
    return output

def p(q,t):
    resns = q.keys()
    mat_as_list = list()
    for row_name in resns:
        mat_as_list.append([q[row_name][col_name] for col_name in resns])
    mat = scipy.matrix(mat_as_list)
    p_as_mat = scipy.linalg.expm(t*mat)
    output = dict((resn, None) for resn in resns)
    for key in output.keys():
        output[key] = dict((resn, None) for resn in resns)

    numrows, numcols = p_as_mat.shape
    for i, j in itertools.product(range(numrows), range(numcols)):
        output[resns[i]][resns[j]] = p_as_mat[i,j]

    return output

