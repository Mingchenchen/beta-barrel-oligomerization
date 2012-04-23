from Bio.SubsMat import SeqMat

def read_text_matrix(data_file, average_pairs = False):
    # An modified version of the Biopython function

    # Make a table of entries from the file. Will contain gunked up
    # comments too
    matrix = dict()
    tmp = data_file.read().split('\n')
    tmp = [line.split() for line in tmp]

    # Remove comments (lines beginning with #)
    tmp = filter(lambda x: x[0] != '#', tmp)

    # Remove blank lines
    tmp = filter(lambda x: len(x) > 0, tmp)

    # Build a dictionary
    alphabet = table[0]
    for index, line in enumerate(table):
        row_resn = table[0]
        col_resn = alphabet[index]
        entries = table[1:]

        for entry in entries:
            if average_pairs:
                factor = .5
            else:
                factor = 1
            matrix[(row_resn, col_resn)] = float(entry)*factor

    # delete entries with an asterisk
    for i in matrix.keys():
        if '*' in i:
            del matrix[i]

    return SeqMat(matrix)



def score_match(pair, matrix):
    try:
        output = matrix[pair]
    except KeyError:
        output = matrix[pair[::-1]]
    return output
    
def score_pairwise(seq1, seq2, matrix,
                   gap_opening_penalty, gap_extension_penalty):
    gap = False
    score = 0

    for pair in zip(seq1, seq2):
        if not gap and '-' not in pair:
            score += score_match(pair, matrix)
        elif not gap and '-' in pair:
            gap = True
            score -= gap_opening_penalty
        elif gap and '-' not in pair:
            gap = False
            score += score_match(pair, matrix)
        elif gap and '-' in pair:
            score -= gap_extension_penalty

    return score
