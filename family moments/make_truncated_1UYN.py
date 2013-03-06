import Bio.AlignIO

msa = Bio.AlignIO.read('gonnet aligned/1UYN with 12.1.6.clu', 'clustal')
template_seq = msa[-3]
first_col = template_seq.seq.find('AATVYA')
assert first_col == 1893
trunc_msa = msa[:,1893:]
Bio.AlignIO.write(trunc_msa, 'truncated 1UYN with 12.1.6.fasta', 'fasta')
