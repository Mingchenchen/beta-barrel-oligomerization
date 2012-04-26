from Bio import AlignIO

def z_diff_blist(target, alignment_path, z_diff_path,
                 alignment_format = 'clustal'):
    # Make a list of residue ids
    stored.resis = list()
    cmd.iterate(target + ' & n. CA',
                'stored.resis.append(resi)')
    assert len(stored.resis) == len(set(stored.resis)),\
           'duplicates in resi list'
    
    # Retrieve sequences as MSA object
    with open(alignment_path, 'r') as f:
        alignment = AlignIO.read(f, alignment_format)
    # Check that first is 1AF6 and the second is 1A0S
    

    # Retrieve z_diffs as list of floats
    with open(z_diff_path, 'r') as f:
        z_diffs = [float(i) for i in f]
   
    blist = list()
    z_differator = iter(z_diffs)
    resiterator = iter(resis)

    for letterA, letterB in zip(alignment[0], alignment[1]):
        if letterA != '-' and letterB != '-':
            blist.append((resiterator.next(), z_differator.next()))

    return blist

def color_b(target, blist):
    cmd.select('on_blist', 'none')
    for resi, b in blist:
        cmd.alter(target + ' & resi ' + resi,
                  'b = ' + str(b))
        cmd.select('on_blist', 'on_blist | resi ' + resi)

def color_by_z_diff(target, alignment_path, z_diff_path,
                    alignment_format = 'clustal'):
    blist = z_diff_blist(target, alignment_path, z_diff_path,
                         alignment_format = alignment_format)
    color_b(target, blist)
    return len(blist)

print(color_by_z_diff('aligned_1AF6',
                      'cluster73 bbtmall 1a0s 1af6.aln',
                      'known 1AF6, unknown 1A0S z_diff, msa, bbtmall.csv'))
