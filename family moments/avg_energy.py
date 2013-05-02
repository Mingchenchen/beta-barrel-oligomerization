from __future__ import division
import zenergy

execfile('family_moments.py')

def average_ezb(params):
    fams = families(params)
    total_ezb = 0
    total_residues = 0
    for family in fams.values():
        sel = Selection(family, lambda x: x.rel_acc > .2 and x.ss == 'E')
        for res_dos in sel:
            for seq_id in family.msa.keys():
                try:
                    total_ezb += res_dos.ez_b(seq_id)
                    added = True
                except zenergy.NoParameters:
                    pass
                if added:
                    total_residues += 1
    return total_ezb/total_residues

normal = average_ezb('published params.csv')
random = average_ezb('normal random params.csv')

print('Average ez beta with published parameters:')
print(normal)
print('')

print('Average ez beta with the randomly generated parameters:')
print(random)
