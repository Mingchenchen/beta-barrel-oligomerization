from __future__ import division
import os
import re


psum = 0
ssum = 0
sum_ = 0
product = 0
for filename in os.listdir('non ppi residues'):
    match = re.match('(\d...)\.csv', filename)
    if match is not None:
        p = os.stat('non ppi residues/' + filename).st_size
        s = os.stat('../pymol/structures/aligned_{}.pdb'\
                    .format(match.group(1))).st_size

        psum += p
        ssum += s

        sum_ += p + s
        product += p*s


    
print(psum)
print(ssum)
