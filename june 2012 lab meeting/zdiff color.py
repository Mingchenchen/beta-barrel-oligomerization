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
