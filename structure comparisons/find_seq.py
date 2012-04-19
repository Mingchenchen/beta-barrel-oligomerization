from collections import deque
from sundries import one_letter

stored.one_letter = one_letter

def find_seq(target):
    stored.resi_resn = list()
    cmd.iterate('n. ca', 'stored.resi_resn.append((resi, one_letter[resn]))')
    
    last_five = deque(stored.resi_resn[:5])
    for resi, resn in stored.resi_resn[5:]:
        if ''.join([x[1] for x in last_five]) == target:
            return resi
        else:
            last_five.popleft()
            last_five.append((resi, resn))
    
    return ''.join([x[1] for x in last_five]) 
