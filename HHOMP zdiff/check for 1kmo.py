from Bio import AlignIO
from Bio import SeqIO

def hdist(a,b):
	match = 0
	total = 0
	a = filter(lambda x: str(x) != '-', a)
	b = filter(lambda x: str(x) != '-', b)
	for i, j in zip(a, b):
		if str(i).upper() == str(j).upper():
			match += 1
		total += 1
	return float(match)/total

# Expected: 1 (if the shorter one is of length n, only considers the first n
# characters of the longer one)
print(hdist('ab','abc'))
# Expected: .5
print(hdist('ab','ac'))


x = AlignIO.read('clusters/cluster91.clu', 'clustal')
y = SeqIO.read('1KMO.fasta', 'fasta')
z = [hdist(y, seq) for seq in x]

z = sorted(z)
print(z[0])
print(z[-1])
