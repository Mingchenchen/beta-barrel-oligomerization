import sys
from Bio import AlignIO
import seqtools

class Report(object):
    '''Object storing the output of this script'''
    def __init__(self):
        self.conj_id = None
        self.disj_id = None

msa_path = sys.argv[1]
msa = AlignIO.read(msa_path, 'clustal')

# Find the sequence that contains the second argument
target_snippet = sys.argv[2]
possible_targets = list(filter(lambda x: target_snippet in x.id,
                               msa))
assert len(possible_targets) == 1,\
       "more than one seq with target snippet in id"
target = possible_targets[0]

# Compare sequence "target" with each sequence in "msa" (including itself)
output = Report()
output.conj_id = seqtools.identity(target, msa,
                                   gap_setting = 'ignore conjunction')
output.disj_id = seqtools.identity(target, msa,
                                   gap_setting = 'ignore disjunction')

# Display output:
print('Gap conjunction ignored identity:')
for identifier, id_percent in output.conj_id:
    print(identifier + ': ' + str(id_percent))
print('')
print('Gap disjunction ignored identity:')
for identifier, id_percent in output.disj_id:
    print(identifier + ': ' + str(id_percent))
