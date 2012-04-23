from __future__ import division

from sundries import CIDict

def retrieve_matrix(csvfile):
    '''
    retrieve_matrix(csvfile)
    Returns a CIDict

    Given a csv file whose rows and columns are labeled, returns a 
    case insensitive dictionary of case insensitive dictionaries
    that will retrieve numberical elements in the csvfile as
    floating point numbers, using the syntax
    dictionary['columnlabel']['rowlable']

    The intention is to do dictionary['original']['mutated'] to get
    a score for a mutation, if the originals are in the rows then
    I did it wrong
    '''
    # Create a dictionary mapping column numbers to resns
    colnum_to_resn = dict()
    # The first line is going to have an empty space as its first entry,
    # so start at the second entry
    for index, entry in enumerate(csvfile.next()[1:]):
        colnum_to_resn.update({index: entry})
    
    # Create a dictionary so dict['A']['V'] will give you score for
    # alanine mutating into valine
    first_to_second = CIDict([(resn, CIDict())
                            for resn in colnum_to_resn.values()])
    for line in csvfile:
        second_resn = line[0]
        for colnum, number in enumerate(line[1:]):
            first_resn = colnum_to_resn[colnum]
            # Add number as floating point so I can sum these things up
            # to get distances
            first_to_second[first_resn].update({second_resn: float(number)})
    
    return first_to_second

def distance(seq1, seq2, matrix):
    '''
    distance(seq1, seq2, matrix)
    Returns a number (floating point if using a matrix generated
                      by retrieve_matrix)
    Get a distance score using a given matrix for a pair of aligned sequences
    Or it might be a similarity score, depending on the matrix.
    '''
    if len(seq1) != len(seq2):
        raise ValueError('Sequences of different length, must compare ' + \
                         'aligned sequences')

    distance = 0
    for pos in range(len(seq1)):
        # These matrices are not necessarily symmetric, order matters
        distance += matrix[seq1[pos]][seq2[pos]]

    return distance

def compare(seq1, seq2, matrix):
    '''
    compare(seq1, seq2, matrix)
    Returns a number

    Just uses 'distance' function with same arguments, and then divides
    by the number of positions in the two sequences that are not both
    gaps
    '''
    count = 0
    for pos in range(len(seq1)):
        if seq1[pos] != '-' or seq2[pos] != '-':
            count += 1

    return distance(seq1, seq2, matrix) / count
