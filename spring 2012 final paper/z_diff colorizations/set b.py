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
    scry = str(alignment[1].seq)
    try:
        assert str(filter(lambda x: x!='-', maltoporin))[:5]=='VDFHG'
        assert str(filter(lambda x: x!='-', scry))[:5]=='SGFEF'
    except AssertionError:
        # Switch the sequences and try again
        scry, maltoporin = maltoporin, scry
        assert str(filter(lambda x: x!='-', maltoporin))[:5]=='VDFHG'
        assert str(filter(lambda x: x!='-', scry))[:5]=='SGFEF'

    # Retrieve z_diffs as list of floats
    with open(z_diff_path, 'r') as f:
        z_diffs = [abs(float(i)) for i in f]
   
    blist = list()
    z_differator = iter(z_diffs)
    resiterator = iter(stored.resis)

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

def auto(name, rainbow='cyan_white_magenta',excluded_color='black'):
    output = color_by_z_diff('aligned_1A0S', auto_list[name][0], auto_list[name][1], alignment_format = auto_list[name][2])
    cmd.spectrum('b',rainbow,'on_blist',0,7)
    cmd.color(excluded_color,'!on_blist')
    return output

sc_dir = r'C:\cygwin\home\alex\beta barrels\structure comparisons' + '\\'
ga_dir = r'C:\cygwin\home\alex\beta barrels\spring 2012 final paper\gonnet alignment' + '\\'
auto_list = {'structural':\
             (sc_dir+'Swiss-PDB structural alignment.fasta',
              sc_dir+'known 1A0S, unknown 1AF6 z_diff, swiss pdb'\
              +' structure alignment.csv',
              'fasta'),
             'bbtmout':\
             (sc_dir+'extracted 1a0s 1af6 from cluster73 bbtm.aln',
              sc_dir+'known 1A0S, unknown 1AF6 z_diff, msa, bbtmout.csv',
              'clustal'),
             'bbtmall':\
             (sc_dir+'extracted 1a0s 1af6 from cluster73 bbtmall.aln',
              sc_dir+'known 1A0S, unknown 1AF6 z_diff, msa, bbtmall.csv',
              'clustal'),
             'gonnet':
              (ga_dir+'1a0s 1af6 gonnet.aln',
               ga_dir+'known 1A0S, unknown 1AF6 z_diff,'\
               +' msa gonnet alignment.csv',
               'clustal')}

