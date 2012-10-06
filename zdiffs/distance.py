from __future__ import division
import Bio.AlignIO
import glob
import re

def identity(seq1, seq2):
    '''Return fraction identical between two sequences'''
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
    
def identity_of_pair(alignment, marker1, marker2):
    '''Given an alignment, return the sequence identity of
    two sequences that contain, in their names or id's, the given markers.
    Raise an AmbiguousMarker exception if either marker matches two
    or more sequences, and a MeaninglessMarker exception if either marker
    matches no sequence.'''
    # Find the sequences
    seq1 = find_in_alignment(alignment, marker1)
    seq2 = find_in_alignment(alignment, marker2)

    # Return the sequence identity as a number from 0 to 1
    return identity(seq1, seq2)

# Find out how many sequences in each cluster are closer to the 
# template than the protein used in the zdiff comparison
for filename in glob.glob('gonnet_aligned/cluster*'):
    # Discover the cluster number, target, and template in the alignment
    info_re = r"cluster(\d\d), (....) as target, (....) as template"
    clusternum, target, template = re.search(info_re, filename).groups()
    
    # Parse the alignment
    alignment = Bio.AlignIO.read(filename, 'clustal')
    # Announce which alignment is being examined
    print('Examining cluster ' + str(clusternum))

    # Report the sequence identity between the target and the template
    temp_targ_id = identity_of_pair(alignment, 'target', 'template')
    print('Identity between zdiff test protein {} and template {} is {}'\
          .format(target, template, temp_targ_id))
    
    # Count how many have at least this much sequence identity with the
    # template
    # identities_with_template returns a list sorted from lowest to 
    # highest, so this reduces to the problem of finding out how far in
    # the list is the target
    id_list = identities_with_template(alignment, 'template')
    # Search through just the names of the sequences
    for index, name in enumerate(x[0] for x in id_list):
        if 'target' in name:
            target_index = index
            break
    number_above = len(id_list) - target_index
    # Fractional coverage is number above plus one over total. The plus
    # one is to include the target sequence itself, which obviously is part
    # of the coverage.
    fractional_coverage = (number_above+1) / len(alignment)
    # Report the results
    print('Number of sequences closer to the template than the protein '
          + 'for which zdiff was calculated: ' + str(number_above))
    print('Total number of sequences in the cluster: '\
          + str(len(alignment)))
    print('Fractional coverage: ' + str(fractional_coverage))

    # Print a blank line as a separator before the results for the next
    # cluster
    print('')
        
