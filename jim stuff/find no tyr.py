from __future__ import division

def counts():
    for pdbid, group in groupdict.items():
        count = cmd.count_atoms(pdbid + '.below_18' \
                                + pdbid + ' & {}.non_ppi ' \
                                + pdbid +  '& {}.not_removed' \
                                + ' & resn tyr & n. ca')
        print((pdbid, count, group.ni_str_count, count/group.ni_str_count))

cmd.extend('counts', counts)

def tyr_color(pdbid):
    cmd.color('forest', pdbid)
    cmd.color('green', pdbid + ' & resn tyr')
    counted = '{0}.non_ppi & {0}.not_removed & {0}.below_18'.format(pdbid)
    cmd.color('deeppurple', counted)
    cmd.color('lightmagenta', counted + ' & resn tyr')
