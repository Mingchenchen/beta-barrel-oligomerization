from __future__ import division
import re

def identity(seq1, seq2):
    '''Return fraction of letters identical between two sequences'''
    # Check that the sequences are the same length - this will always
    # be the case for aligned sequences
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must be the same length')

    # Count number of positions that are the same, and number of positions
    # that are different. Skip positions in which both are gaps.
    same = 0
    different = 0
    for top, bottom in zip(seq1, seq2):
        # Skip a pair of gaps
        if top == '-' and bottom == '-':
            continue
        else:
            if top == bottom:
                same += 1
            elif top != bottom:
                different += 1
    
    # Return fraction identical - number between zero and 1
    return same / (same + different)

class AmbiguousMarker(Exception):
    '''Raised by find_in_alignment when a marker is found in the names
    or identifiers of two or more sequences of an alignment'''
    pass

class MeaninglessMarker(Exception):
    '''Raised by find_in_alignment when a marker is not found in the names
    or identifiers of any sequence in an alignment'''
    pass

def find_in_alignment(alignment, marker):
    '''Given an alignment and a "marker" - a substring of the identifier
    of a sequence - return the unique sequence in the alignment that
    has the marker. Raises exceptions if 0 or >1 sequences have the
    marker.'''
    marked_seq = None
    for seq in alignment:
        # I don't really know the difference between names and id's, so
        # check both:
        if marker in seq.name or marker in seq.id:
            # There should only be one
            if marked_seq is not None:
                raise AmbiguousMarker('More than one sequence has "'\
                                        + marker + '" in its '\
                                        + 'name or id')
            else:
                marked_seq = seq
    # There should be more than zero
    if marked_seq is None:
        raise MeaninglessMarker('No sequence has {0} in its name or id'\
                         .format(marker))
    
    return marked_seq
    
def identities_with_template(alignment, template_identifier):
    '''Return a list of each sequence identifier, along with its sequence
    identity with the sequence of a template structure.  The template
    structure's sequence must be in the alignment, and must be the only
    sequence with the string given as the second argument contained within
    its name or id.
    Returned list is sorted from lowest sequence identity to highest.
    Raises a AmbiguousMarker exception if more than one sequence is
    identified as a template, and a MeaninglessMarker exception of no
    sequences are.'''
    # Find the template sequence
    template_seq = find_in_alignment(alignment, template_identifier)

    # Make a list of sequence identities with the template
    identities = list()
    for seq in alignment:
        identities.append((seq.id, identity(template_seq, seq)))

    return sorted(identities, key=lambda x: x[1])
