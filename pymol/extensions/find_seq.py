from collections import deque
from sundries import one_letter

def find_seq(target, silent = False):
    '''Returns a list of all positions at which a given 5-residue
    sequence is found.'''
    
    stored.one_letter = one_letter
    
    # Turn the sequence into a list of (resi, resn) pairs
    # This way, when the target sequence of resn's is found,
    # the corresponding resi's can be returned
    stored.resi_resn = list()
    cmd.iterate('n. ca',
                'stored.resi_resn.append((resi, one_letter[resn]))')
    
    # Create a variable to hold the five residues currently being checked,
    # a sort of sliding window over the sequence
    last_five = deque(stored.resi_resn[:5])
    
    # Loop over the sequence, moving the sliding window
    # Append a resi to the output list whenever the target is found
    output = list()
    for resi, resn in stored.resi_resn[5:]:
        if ''.join([x[1] for x in last_five]) == target:
            output.append(last_five[0][0])

        # Shift the sliding window one to the right
        last_five.popleft()
        last_five.append((resi, resn))
    
    if not silent:
        for i in output:
            print(i)
    
    return output

cmd.extend("find_seq", find_seq)

class MultipleMatches(Exception):
    pass
    
class NoMatches(Exception):
    pass

def find_one(target):
    pos_list = find_seq(target, silent = True)
    if len(pos_list) > 1:
        raise MultipleMatches("Multiple matches for {}: {}" \
                              .format(target, ", ".join(pos_list)))
    if len(pos_list) == 0:
        raise NoMatches("No matches for " + target)
    return pos_list[0]
    
def select_strand(name, seq1, seq2):
    '''Given a 5-residue sequence that starts on the first residue of a
    strand, and a 5-residue sequence that starts on the last residue of a
    strand, selects the strand and colors it red'''

    start = find_one(seq1)
    end = find_one(seq2)
    cmd.select(name, "i. {}-{}".format(start, end))
    print("i. {}-{}".format(start, end))
    cmd.color("red", name)
    
cmd.extend("select_strand", select_strand)