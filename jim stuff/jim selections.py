# Open the asymmetric ezbeta spreadsheets used for retrieving DSSP results
# and the residue numbers of the residues included in the dataset
phrasebooks = biodata.phrasebooks('weights phrasebook.csv')
weights_phrasebook = phrasebooks['weights']
for pdbid, group in groupdict.items():
    csv_name = pdb.upper() + '.csv'
    if csv_name in os.listdir('non_ppi'):
        group.non_ppi_data = biodata.Spreadsheet('non_ppi/' + csv_name,
                                            phrasebook = weights_phrasebook)
    if csv_name in os.listdir('ppi'):
        group.ppi_data = biodata.Spreadsheet('ppi/' + csv_name,
                                        phrasebook = weights_phrasebook)

def select_by_predicate(name, function, arg_list, superset):

    def conditional_add(resi, *args):
        if function(*args):
            cmd.select(name, '{} | ({} & i. {})'\
                             .format(name, superset, resi))
    
    stored.conditional_add = conditional_add

    cmd.select(name, 'none')


    cmd.iterate_state(1, superset,
                      'conditional_add(resi, '\
                      + ', '.join(arg_list) \
                      + ')')



# Select top/bottom halves

def extracellular_half(z):
    return z>0
def periplasmic_half(z):
    return z<0

for pdbid in groupdict.keys():
    select_by_predicate(pdbid + '.ec_half', extracellular_half,
                       ['z'], pdbid + '.molecule')
    select_by_predicate(pdbid + '.pp_half', periplasmic_half,
                       ['z'], pdbid + '.molecule')

# Select a range

def from_6_to_18(z):
    return abs(z) > 6 and abs(z) < 18

for pdbid in groupdict.keys():
    select_by_predicate(pdbid + '.6to18', from_6_to_18,
                       ['z'], pdbid + '.molecule')


# Select residues in the ppi and non_ppi datasets
def in_non_ppi_dataset(resi, group):
    return resi in group.non_ppi_data.get_column('resi')

def in_ppi_dataset(resi, group):
    return resi in group.ppi_data.get_column('resi')

for pdbid, group in groupdict.items():
    stored.group = group
    select_by_predicate(pdbid + '.non_ppi', in_non_ppi_dataset,
                        ['resi', 'stored.group'], pdbid + '.molecule')
    select_by_predicate(pdbid + '.ppi', in_ppi_dataset,
                        ['resi', 'stored.group'], pdbid + '.molecule')
    



# I only have DSSP results for residues included in the asymmetric ezB
# dataset. 
# Select beta sheet residues WITHIN the non_ppi dataset

def ss(group, target_resi):
    '''ss(group,resi)
    retrieves secondary structure using non ppi dataset spreadsheet'''
    for resi, ss in group.non_ppi_data.get_columns(['resi','ss']):
        if resi != '' and int(resi)==int(target_resi):
            return ss

# Select beta sheets
def beta_sheet(group, resi):
    return ss(group, target_resi) == 'E'

for pdbid, group in groupdict.items():
    stored.group = group
    select_by_predicate(pdbid + '.beta_sheet', beta_sheet,
                       ['stored.group', 'resi'], pdbid + '.non_ppi')
