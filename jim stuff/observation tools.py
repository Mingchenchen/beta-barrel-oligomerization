from __future__ import division

def counts(resn = 'tyr'):
    for pdbid, group in groupdict.items():
        count = cmd.count_atoms(pdbid + '.below_18 & ' \
                                + pdbid + '.non_ppi & ' \
                                + pdbid +  '.not_removed & ' \
                                + 'resn {} & n. ca'.format(resn))
        print((pdbid, count, group.ni_str_count, count/group.ni_str_count))
    print('(pdbid, count, strands, count/strands)')

cmd.extend('counts', counts)

def aro_color(pdbid, resn = 'tyr'):
    cmd.hide('everything', pdbid)
    cmd.show('cartoon', pdbid)

    cmd.color('forest', pdbid)
    cmd.color('limon', pdbid + ' & resn ' + resn)

    counted = '{0}.non_ppi & {0}.not_removed & {0}.below_18'.format(pdbid)
    cmd.color('deeppurple', counted)
    cmd.color('violet', counted + ' & resn ' + resn)

    cmd.color('grey30', '{0} & !{0}.below_18'.format(pdbid))
    cmd.show('sticks', pdbid + ' & resn ' + resn)
    

cmd.extend('arocolor', aro_color)
