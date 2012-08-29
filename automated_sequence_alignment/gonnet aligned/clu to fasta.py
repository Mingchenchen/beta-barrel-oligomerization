import os
from Bio import AlignIO

print os.listdir('.')
for i in os.listdir('.'):
    if i[-4:] == '.clu':
        msa = AlignIO.read(i, 'clustal')
        AlignIO.write(msa, i[:-4] + '.fasta', 'fasta')
