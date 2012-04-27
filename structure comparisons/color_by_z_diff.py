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
    maltoporin = str(alignment[0].seq)
    assert str(filter(lambda x: x!='-', maltoporin))[:5]=='VDFHG',\
           'first seq not maltoporin'
    scry = str(alignment[1].seq)
    assert str(filter(lambda x: x!='-', scry))[:5]=='SGFEF',\
           'second seq not sucrose porin'

    # Retrieve z_diffs as list of floats
    with open(z_diff_path, 'r') as f:
        z_diffs = [abs(float(i)) for i in f]
   
    blist = list()
    z_differator = iter(z_diffs)
    resiterator = iter(stored.resis)


    print(maltoporin)


    for malt_resn, scry_resn in zip(maltoporin, scry):
        if malt_resn != '-' and scry_resn != '-':
            blist.append((resiterator.next(), z_differator.next()))
        elif malt_resn == '-' and scry_resn != '-':
            # No maltoporin residue corresponding to this scry position
            resiterator.next()
        elif malt_resn != '-' and scry_resn == '-':
            # A maltoporin residue that doesn't map onto anywhere
            pass

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

print(color_by_z_diff('aligned_1A0S',
                      'cluster73 bbtmall 1a0s 1af6.aln',
                      'known 1A0S, unknown 1AF6 z_diff, msa, bbtmall.csv'))
