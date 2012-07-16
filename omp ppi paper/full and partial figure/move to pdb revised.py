import numpy as np

# Store the coordinates of each residue in the tobmodel model
stored.tob_coords = dict()
cmd.iterate_state(1, 'topranked*',
                'stored.tob_coords.update({resi: np.array([x,y,z])})')

# Make a mapping between pdb structure residue numbers and
# tobmodel residue numbers
stored.count = 1
tob_to_pdb = dict()
cmd.iterate('pdb & n. ca',
            'tob_to_pdb.update({str(stored.count): resi}); stored.count+=1')
pdb_to_tob = dict([(value, key) for key, value in tob_to_pdb.items()])


for tob_resi, new in stored.tob_coords.items():
    # Find the traslation from the old to the new position, apply
    # it to each atom in the residue
    pdb_resi = tob_to_pdb[tob_resi]
    cmd.iterate_state(1, 'pdb & n. ca & i. ' + pdb_resi,
                      'stored.old = np.array([x,y,z])')
    stored.trnslt = new - stored.old 
    cmd.alter_state(1, 'pdb & i. ' + pdb_resi,
                    'x = float(x + stored.trnslt[0]);'\
                    + 'y = float(y + stored.trnslt[1]);'\
                    + 'z = float(z + stored.trnslt[2])')

stored.count = 0
def delete_me_maybe(pdb_resi):
    if pdb_to_tob[pdb_resi] not in stored.tob_coords.keys():
        cmd.remove('pdb & i. ' + pdb_resi)
        stored.count += 1

cmd.iterate('pdb & n. ca', 'delete_me_maybe(resi)')
print('deleted ' + str(stored.count) + ' residues')
