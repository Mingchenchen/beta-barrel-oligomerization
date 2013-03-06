import Bio.AlignIO

msa = Bio.AlignIO.read('better 1UYN with 12.1.6.clu', 'fasta')
template_seq = msa[-3]
first_col = template_seq.seq.find('AATVYA')
trunc_msa = msa[:,first_col:]
Bio.AlignIO.write(trunc_msa, 'better truncated 1UYN with 12.1.6.fasta', 'fasta')
