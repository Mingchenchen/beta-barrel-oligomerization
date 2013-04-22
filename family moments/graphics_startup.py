fams = families()
rand_param_fams = families('normal random params.csv')

def beta_exposed(pdbid):
    return Selection(fams[pdbid], lambda y: y.ss == 'E' and y.rel_acc > .2)

def rand_param_beta_exposed(pdbid):
    return Selection(rand_param_fams[pdbid], lambda y: y.ss == 'E' and y.rel_acc > .2)
def beta(pdbid):
    return Selection(fams[pdbid], lambda y: y.ss == 'E')

def pretty():
    cmd.bg_color('white')
    cmd.set('opaque_background', 'off')
    cmd.hide('everything')
    cmd.show('cartoon', 'polymer & ss s')
    cmd.color('yellow', 'polymer')

    cmd.show('everything', 'moment*')
    cmd.hide('labels')
    cmd.color('blue', 'moment*')
