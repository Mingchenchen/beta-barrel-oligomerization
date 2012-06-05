def ss(group, target_resi):
    '''ss(group,resi)
    retrieves secondary structure using non ppi dataset spreadsheet'''
    for resi, ss in group.non_ppi.get_columns(['resi','ss']):
        if resi != '' and int(resi)==int(target_resi):
            return ss

def select_from_list(resi_lists, names, superset):
    for resi_list, name in zip(resi_lists, names): 
        cmd.select(pdbid + name, 'none')
        for resi in resi_list:
            cmd.select(name, '{} | (i. {} & {})'\
                             .format(name, resi, superset))
            
# Define half selections
for pdbid, group in groupdict.items():
    # Make lists of the residues to be included in each selection
    stored.ec_half = set()
    stored.pp_half = set()
    
    def add(resi, z):
        if z>0:
            stored.ec_half.add(resi)
        elif z<0:
            stored.pp_half.add(resi)
       
    stored.add = add
    iterate_state(1, pdbid + '.molecule', 
                  'stored.add(resi, z)')
    
    # Make the selection
    names = ('ec_half', 'pp_half')
    names = [pdbid + '.' + x for x in names]
    resi_lists = [ec_half, pp_half]
    select_from_list(resi_lists, names, pdbid + '.molecule')

# Define the region selections
lower_limit = 6
upper_limit = 18
for pdbid, group in groupdict.items():
    # Make lists of the resiues to be included in each selection
    stored.upper_region = set()
    stored.lower_region = set()
    def add(resi, z):
        if abs(z) > lower_limit and abs(z) < upper_limit:
            if z > 0:
                stored.upper_region.add(resi)
            elif z < 0:
                stored.lower_region.add(resi)
    stored.add = add
    iterate_state(1, pdbid + '.molecule',\
                  'stored.add(resi, z')
    
    # Make the selection
    names = ('ec_{}-{}'.format(lower_limit, upper_limit),
             'pp_{}-{}'.format(lower_limit, upper_limit))
    names = [pdbid + '.' + name for name in names]
    resi_lists = (stored.upper_region, stored.lower_region)
    select_from_list(resi_lists, names, pdbid + '.molecule')
