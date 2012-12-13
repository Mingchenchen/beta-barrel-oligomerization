from __future__ import division

# Go to the correct directory so that zdiff_report will work
cmd.cd(r'C:\cygwin\home\alex\beta-barrel-oligomerization\zdiffs')

def zdiff_blist(zdiff_path):
    '''Return a list of pairs (resi, abs(zdiff)); that is, identification
    numbers of residues in the protein whose structure was predicted,
    paired with the error in the z value of the prediction of that
    residue's Calpha'''
    blist = list()
    with open(zdiff_path, 'r') as f:
        for line in f:
            # Ignore comments
            if line[0] == '#':
                continue
            # Extract data
            resi_str, zdiff_str = line.split(',')
            resi = int(resi_str)
            b = abs(float(zdiff_str))
            blist.append((resi, b))
    return blist

def set_b(target, blist):
    cmd.select('has_data', 'none')
    for resi, b in blist:
        cmd.alter('{} & resi {}'.format(target, resi),
                  'b = ' + str(b))
        cmd.select('has_data', 
                   'has_data | ({} & resi {})'.format(target, resi))

def color_by_zdiff(target, zdiff_path, max_,
                   spectrum='cyan_white_magenta', no_data_color = 'black'):
    '''Color a given PyMOL selection "target" by the zdiffs at "zdiff_path".
    A selection called "has_data" will be created, containing all residues
    for which there was a value in the zdiff file.

    Those for which their is data will be colored by the given spectrum
    with the given maximum. Those for which there was no data will
    be colored the given no_data_color.

    Will fail if not every resi in the zdiff file is within the target.'''
    blist = zdiff_blist(zdiff_path)
    set_b(target, blist)

    cmd.spectrum('b', spectrum, '{} & has_data'.format(target),
                 0, max_)
    cmd.color(no_data_color, '{} & ! has_data'.format(target))

    return len(blist)

def zdiff_report(target_pdbid, max_=3):
    '''Create, in PyMOL, a visual representation of the Gonnet matrix's
    success in predicting the z coordinates of the protein with the given
    PDBID. Produces the structure of the protein with those that had
    an error in z coordinate of 0 colored cyan, those that had an error
    in z coordinate of 3 or greater colored magenta, and those with 
    errors between 0 and 3 colored on a spectrum from cyan to white to 
    magenta. Residues for which no prediction was made are colored black.

    Also opens the template structure that was used.

    The max_ keyword argument (default 3) sets the error to which the color
    magenta is assigned.

    pdbid must be in (1BY5, 1T16, 2MPR, 1PHO, 2OMF). See "zdiffs/zdiff
    notes.docx" September 4th entry for the reason these were chosen.

    Must be run with "zdiffs" folder as working directory.'''

    # If run as a PyMOL command with a value of max specified, this will
    # be passed as a string
    # Convert it to a number:
    max_ = float(max_)

    # Open the real structure of the protein whose structure was predicted
    target_id = target_pdbid.upper()
    cluster_map = dict({'1BY5': 'cluster18',
                        '1T16': 'cluster71',
                        '2MPR': 'cluster73',
                        '2OMF': 'cluster99',
                        '1PHO': 'cluster99'})

    template_map = dict({'1BY5': '3EFM',
                         '1T16': '3BS0',
                         '2MPR': '1A0S',
                         '2OMF': '2J1N',
                         '1PHO': '2J1N'})

    # Check that the given PDBID is one of those for which zdiff was
    # calculated
    if target_id not in cluster_map.keys():
        raise ValueError('zdiffs not calculated for that structure')
        
    cluster = cluster_map[target_id]
    template_id = template_map[target_id]

    target_path = "comparison structures/{}/{} aligned to daniel's {}.pdb"\
                  .format(cluster, target_id, template_id)
    cmd.load(target_path)

    # Color the structure of the prediction target
    # If PyMOL is ever updated to allow spaces or apostrophes
    # in object names, this code will probably stop working
    target_name = '{}_aligned_to_daniel_s_{}'.format(target_id, template_id)
    zdiff_path = 'gonnet_zdiff/{}, {} as target, {} as template.zdiff'\
                 .format(cluster, target_id, template_id)
    color_by_zdiff(target_name, zdiff_path, max_)

    # Open the template structure
    template_path = 'comparison structures/{}/aligned_{}.pdb'\
                    .format(cluster, template_id)
    cmd.load(template_path) 

cmd.extend('zdiff_report', zdiff_report)

def fraction_below(f, condition = lambda z, ss: z<=15 and ss=='S'):
    '''Prints what fraction of predicted z values were within f of the
    real z value, and for what fraction of residues there was a prediction.

    By default, only checks beta strand residues within fifteen Angstroms
    of the membrane center. You can change this by changing the "condition"
    keyword argument. The "condition" argument is a function of z and
    ss (the properties of atoms to which cmd.iterate_state has access).
    For example, the default setting is lambda z,ss: z<=15 and ss='S'.'''
    cmd.select('to_be_checked', 'none')
    def check(z, ss, resi):
        if condition(z, ss):
            cmd.select('to_be_checked', 'to_be_checked | i. ' + resi)
    stored.check = check
    cmd.iterate_state(1, '*_aligned_to_*', 'stored.check(z, ss, resi)')

    if cmd.count_atoms('to_be_checked') == 0:
        raise ValueError('No c alphas match condition')

    # Helper function for calculating conditional frequency of b<f
    def p(property_, background):
        return cmd.count_atoms(property_ + '&' + background)\
               /cmd.count_atoms(background)
    print('Fraction closer than {} Angstroms:'.format(f))
    print(p('b<'+str(f), '*_aligned_to_* & n. ca & has_data & to_be_checked'))
    print('Fraction for which no prediction was made:')
    print(p('!has_data', '*_aligned_to_* & to_be_checked'))

cmd.extend('fraction_below', fraction_below)


def check_all_for_all():
    for pdbid in ('1BY5', '1T16', '2MPR', '1PHO', '2OMF'):
        cmd.reinitialize()
        zdiff_report(pdbid)
        for distance in (3, 6):
            print((pdbid, distance))
            # Check fraction below, with no condition (checking each calpha)
            fraction_below(distance, condition = lambda z, ss: True)
            print('')

