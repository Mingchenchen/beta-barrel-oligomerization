stored.coords = dict()
cmd.iterate_state(1, 'topranked*',
                  'stored.coords.update({int(resi): (x,y,z)})')
# Okay, that's the coordinates numbered by position in the
# sequence (starting with 0)

stored.resns = dict()
cmd.iterate('topranked*',
            'stored.resns.update({int(resi): resn})')

def get_if_there(key, dictionary, default):
    if key in dictionary.keys():
        return dictionary[key]
    else:
        return default

stored.count = 0
cmd.alter_state(1, 'pdb & n. ca',
          'x, y, z = get_if_there(stored.count, stored.coords, (x,y,z));' \
          + 'stored.count+=1')
# Now check if this process actually lines them up right:
stored.count = 0
stored.check = list()
cmd.alter('pdb & n. ca',
          'stored.check.append(resn == ' \
          + 'get_if_there(stored.count, stored.resns, True));' \
          + 'stored.count+=1')
if all(stored.check) & len(stored.check)>0:
    print('probably worked')
else:
    print('probably didn\'t work')
'''
# As an experiment. Move all atoms to the position of their corresponding Calpha.
stored.ca_coord = dict()
cmd.iterate_state(1,'pdb &  n. ca', 'stored.ca_coord.update({resi: (x,y,z)})')
cmd.alter_state(1, 'pdb & ! n. ca', 'x, y, z = stored.ca_coord[resi]')
'''