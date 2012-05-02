from __future__ import division

class MissingPDBSequence(Exception):
    pass

class AlignmentOracle(object):
    '''
    Initialized with a multiple sequence alignment and a name to find a
    particular sequence of interest. Can give amino acids from other
    sequences in the alignment that correspond to amino acids in the
    sequence
    of interest, ignoring amino acids that are lined up with gaps in the
    sequence of interest.
    '''.replace('  ', '')
    def __init__(self, alignment, pdb_name = 'pdb'):
        self.data = alignment
        self.pdb_index = None
        for index, record in enumerate(self.data):
            if pdb_name in record.id:
                self.pdb_index = index
                break
        if self.pdb_index == None:
            raise MissingPDBSequence('None of the sequences had '+\
            pdb_name + ' in their title')
        self._pdb_sequence = self.data[self.pdb_index].seq

    def pdb_sequence(self, selection = None):
        return self.sequence(self.pdb_index, selection)

    def sequence(self, index, selection = None):
        pos = 0
        for j, letter in enumerate(self.data[index].seq):
            if self._pdb_sequence[j] != '-':
                if selection is None or pos in selection:
                    yield letter.upper()
                pos += 1


def pairwise_identity(seq1, seq2, gap_setting='ignore conjunction'):
    '''Return % sequence id between two sequences, disregarding columns
    where both have gaps
    gap_setting='ignore conjunction' (default):
        ignore columns where both have gaps
    gap_setting='ignore disjunction':
        ignore columns where either has a gap'''

    if len(seq1) != len(seq2):
        raise ValueError('pairwise_identity function '\
                        + 'only for aligned sequences')
    total = 0
    same = 0
    for pos1, pos2 in zip(seq1, seq2):
        if gap_setting == 'ignore conjunction':
            if pos1 == '-' and pos2 == '-':
                continue
        elif gap_setting == 'ignore disjunction':
            if pos1 == '-' or pos2 == '-':
                continue
        else:
            raise ValueError('not recognized gap setting')

        total += 1
        if pos1 == pos2:
            same += 1
    return same/total


def identity(reference_seq, msa, gap_setting='ignore conjunction'):
    '''Returns % sequence id between a given sequence and each sequence
    in an MSA. The subroutine calculating pairwise identity disregards
    columns where both have gaps.
    gap_setting='ignore conjunction' (default):
        ignore columns where both have gaps
    gap_setting='ignore disjunction':
        ignore columns where either has a gap'''

    output = list()
    for aligned_seq in msa:
        output.append((aligned_seq.id,
                       pairwise_identity(aligned_seq, reference_seq,
                                         gap_setting = gap_setting)))
    output = sorted(output, key=lambda x: x[1])
    return output
