import csv
from Bio import AlignIO

def z_diff_blist(target, alignment_path, z_diff_path,
                 alignment_format = 'clustal'):
    # Make a list of residue ids
    stored.resis = list()
    cmd.iterate(target + 'n. CA',
                'stored.resis.append(resi)')
    assert len(stored.resis) == len(set(stored.resis)), \ 
           'duplicates in resi list'
    
    # Retrieve sequences as MSA object
    with open(alignment_path, 'r') as f:
        alignment = AlignIO.read(f, alignment_format)

    # Retrieve z_diffs as list of floats
    with open(z_diff_path, 'r') as f:
        z_diffs = [float(i) for i in csv.reader(f)]
   
    blist = list()
    z_differator = iter(z_diffs)

    for resi, letterA, letterB \
        in zip(stored.resis, alignment[0], alignment[1]):
        if alignment[0] != '-' and alignment[1] != '-':
            blist.append(resi, z_differator.next())

    return blist

def color_b(target, blist):
    for resi, b in blist:
        cmd.alter(target + ' & resi ' + resi,
                  'b = ' + b)

def color_by_z_diff(target, alignment_path, z_diff_path,
                    alignment_format = 'clustal'):
    
