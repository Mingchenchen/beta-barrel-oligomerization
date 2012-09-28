import Bio.AlignIO

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must be the same length')
    
