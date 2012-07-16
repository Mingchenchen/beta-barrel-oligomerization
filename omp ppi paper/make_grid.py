cmd.disable('*.pdb')
cmd.set('grid_mode', '1')
for index, protein in enumerate(oligomers.keys()):
    cmd.set('grid_slot', str(index+1), protein + '.molecule')
