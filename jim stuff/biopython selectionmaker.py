'''Make histograms where the data points are protein structures and the
numbers are counts of aromatic residues.
Decent code - the predicate-based selection is clever - but it's really
annoying that the attributes of the group objects are created as the
script goes along. What if you want to know what type they are?
Or get a list of them?
I don't know exactly what to change but that's definitely something
to avoid in code derived from this script. Histogram class is a keeper.
The lack of a function for saving the selections is ridiculous.'''

from __future__ import division
import re
from Bio.PDB import PDBParser
from sundries import CIDict
from sundries import three_letter
import functools
import csv
import os
import warnings
import biodata
import time
import itertools

start = time.time()

class Group(object):
    def __init__(self, name):
        self.name = name

    def selection(self, function):
        return set(filter(functools.partial(function, self),
                          self.structure.get_residues()))
    
# Retrieve all the pdbids that are in the asymmetric dataset
filenames = os.listdir('non ppi residues')
matches = [re.match('(\d...)\.csv', x) for x in filenames]
real_matches = filter(lambda x: x is not None, matches)
pdbids = [x.group(1) for x in real_matches]

# Uncomment this line to speed things up for testing:
pdbids = ['1QD6']

# Make group objects. I'm going to be associating a lot of stuff with
# each protein, it's easiest to just group them together.
groupdict = CIDict([(pdbid, Group(pdbid)) for pdbid in pdbids])

# The slow part: load structures.
structure_dir = '../pymol/structures'
def filename(pdbid):
    return structure_dir + '/aligned_{}.pdb'.format(pdbid)
# Daniel's aligned structures give "invalid/missing occupancy" and
# "invalid/missing B factor" warnings - thousands of them! Have to filter
# warnings or the structures won't get loaded
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for group in groupdict.values():
        group.structure = PDBParser().get_structure(group.name,
                                                    filename(group.name))

print('structures loaded after ' + str(time.time() - start))
afterstructures = time.time()

# Open the asymmetric ezbeta spreadsheets used for retrieving DSSP results
# and the residue numbers of the residues included in the dataset
phrasebooks = biodata.phrasebooks('weights phrasebook.csv')
weights_phrasebook = phrasebooks['weights']
for group in groupdict.values():
    csv_name = group.name.upper() + '.csv'
    if csv_name in os.listdir('non ppi residues'):
        group.non_ppi_data=biodata.Spreadsheet('non ppi residues/'\
                                                 + csv_name,
                                        phrasebook = weights_phrasebook)
    if csv_name in os.listdir('ppi residues'):
        group.ppi_data = biodata.Spreadsheet('ppi residues/'\
                                             +csv_name,
                                        phrasebook = weights_phrasebook)

print('loading spreadsheets took ' + str(time.time() - afterstructures))

# Instead of select_by_predicate, selections will be created by
# using filter on a list of residues
# I want to use them in PyMOL, so:
def selection_file_line(group, selection):
    '''Return a list to be a line in a selection file readable by my
    PyMOL scripts'''
    output = list()
    output.append(group.name)
    output += [residue.get_id()[1] for residue in selection]
    return output

# Each PDB structure has a method get_residues(), which is a generator
# returning residue objects.
# The attributes of residue objects I'm gonna need are:
# residue.get_id()[1] is resi, as an integer
# residue.get_resname() is resn, as a 3-letter string
# residue.child_dict['CA'].get_coord()[2] is the z coordinate
 
# Other information will come from the spreadsheet


# Selection of residues in the non_ppi portion of the asymmetric dataset:
# note three phases. Define predicate. Make selection. Write to file.
def in_non_ppi_dataset(group, residue):
    return residue.get_id()[1] in group.non_ppi_resis

making_selections_start= time.time()
for group in groupdict.values():
    group.non_ppi_resis = set()
    for entry in group.non_ppi_data.get_column('resi'):
        if entry != '':
            group.non_ppi_resis.add(int(entry))
    
    group.non_ppi_res = group.selection(in_non_ppi_dataset)

print('making non_ppi selections took '\
      + str(time.time() - making_selections_start))


with open('selections/non ppi.csv', 'wb') as f:
    csv.writer(f).writerows([selection_file_line(group, group.non_ppi_res)\
                             for group in groupdict.values()])


# Selection of residues above and below z=0
def on_ec_side(group, residue):
    try:
        output = residue.child_dict['CA'].get_coord()[2] > 0
    except KeyError:
        output = False
    return output


def on_pp_side(group, residue):
    try:
        output = residue.child_dict['CA'].get_coord()[2] < 0
    except KeyError:
        output = False
    return output

for group in groupdict.values():
    group.ec_half_res = group.selection(on_ec_side)
    group.pp_half_res = group.selection(on_pp_side)

with open('selections/ec half.csv', 'wb') as ec,\
     open('selections/pp half.csv', 'wb') as pp:
    for group in groupdict.values():
        csv.writer(ec).writerow([group.name] + [residue.get_id()[1]\
                                 for residue in group.ec_half_res])
        csv.writer(pp).writerow([group.name] + [residue.get_id()[1]\
                                 for residue in group.pp_half_res])


# Selection of residues such that |z|<18, to exclude the tall extracellular
# parts of tall proteins

def below_18(group, residue):
    try:
        z = residue.child_dict['CA'].get_coord()[2]
    except KeyError:
        return False
    return abs(z) < 18

for group in groupdict.values():
    group.below_18_res = group.selection(below_18)

with open('selections/below 18.csv', 'wb') as f:
    
    for group in groupdict.values():
        row_to_write = list()
        row_to_write.append(group.name)
        row_to_write += [residue.get_id()[1] \
                         for residue in group.below_18_res]
        csv.writer(f).writerow(row_to_write)


# Selection of residues that were not removed manually. Manually 
# removed residues look like they are on the inside of the protein, or
# on the water-facing rather than lipid-facing part of their surface.
# Very subjective, but a lot of stuff that definitely needs to be removed
# was.

manually_removed = dict()
with open('removed.csv', 'rb') as f:
    for row in csv.reader(f):
        pdbid = row[0]
        resis = (int(i) for i in row[1:])
    
for group in groupdict.values():
    group.manually_removed = set(resis)
    
def not_removed(group, residue):
    return residue.get_id()[1] not in group.manually_removed

for group in groupdict.values():
    group.not_removed_res = group.selection(not_removed)

with open('selections/not removed.csv', 'wb') as f:
    for group in groupdict.values():
        row_to_write = list()
        row_to_write.append(group.name)
        row_to_write += [residue.get_id()[1] \
                         for residue in group.not_removed_res]
        csv.writer(f).writerow(row_to_write)


# Load non-interface strand counts. These will be used to normalize
# the the residue counts.
with open('strand counts olig minus.csv', 'rb') as f:
    ni_str_counts = CIDict((pdbid, int(count)) \
                            for pdbid, count in csv.reader(f))

for group in groupdict.values():
    group.ni_str_count = ni_str_counts[group.name]




class Histogram(list):
    '''A class that inherits from lists. Attempts to set or retrieve
    from indices beyond the length of the list will result in it extending
    with 0 entries until it is long enough to contain the index.

    Initialized with an integer "binsize" argument.

    Method "h.add_data(value)" places a data point in a bin, depending on
    the binsize. Rule is that the entry in index value // binsize is
    incremented by 1.

    Method "h.save(filename)" saves the histogram to a csv file with
    two columns. Left column is first number of bin, right column is count
    in that bin.'''

    def __init__(self, binsize):
        self.binsize = binsize

    def extend_to(self, i):
        '''x.extend(i): extend so that highest index is i
        All entries added are zeros'''
        if len(self) <=i:
            # Its highest index should be i, so its length should be i+1
            self += [0 for n in range(i+1 - len(self))]

    def __getitem__(self, i):
        '''Retrieve item i, or make list bigger until there is an item i'''
        self.extend_to(i)
        return list.__getitem__(self, i)

        
    def __setitem__(self, i, y):
        '''Set item i, extending list so that item i exists'''
        self.extend_to(i)
        return list.__setitem__(self, i, y)
    
    def add_data(self, value):
        self[int(value / self.binsize)] += 1

    def save(self, filename):
        with open(filename, 'wb') as f:
            csv.writer(f).writerows((i * self.binsize, count) \
                                    for i, count in enumerate(self))
           
def histogram(groupdict, binsize, resns, sele_names):
    '''Function to produce Histogram objects from a groupdict
    (a dictionary mapping strings to groups). Will return a histogram
    with frequency of a specified group of residues within a specified
    selection on the x axis, and frequency of proteins in the groupdict
    that have that frequency of residues on the y axis.

    Automatically divides its counts by the number of non-interface
    strands, group.ni_str_count

    groupdict: the groupdict containing groups with attributes that
        are sets of residues ('selections')
    binsize: size of the bins
    resns: List of names of residues to count. Can be uppercase or
        lowercase, one-letter or three-letter convention.
    sele_names: names of the selections counted residues must be a member
        of. Entering three selections, for example, the function will
        count the specified residues in the INTERSECTION of the three
        selections given.

    Remember, a residue is counted if it is ONE of the types in the 
    "resns" container, and in ALL of the selections in the sele_names
    container.'''

    # Change resns to standard format biopython uses
    # (three letters, uppercase)
    # First, change to list since resns must be a mutable object for this
    
    resns = list(resns)

    for i, resn in enumerate(resns):
        if len(resn) == 1:
            resns[i] = three_letter[resn]
    
    for i, resn in enumerate(resns):
        resns[i] = resn.upper()

    # Create the histogram object that will be returned
    output = Histogram(binsize)

    # For each protein, count the number of residues in the intersection
    # of the given selections that are of one of the given types.
    for group in groupdict.values():
        protein_count = 0
        
        # Create a set of residue objects corresponding to the set
        # of given selections
        selections = [group.__getattribute__(x) for x in sele_names]
        intersection = set.intersection(*selections)
        
        # Check if it is one of the given types
        for residue in intersection:
            if residue.get_resname() in resns:
                protein_count += 1
        
        # Divide by non interface strand count
        normalized_count = protein_count/group.ni_str_count
        # Increment the appropriate bin of the histogram
        output.add_data(protein_count)

    return output

for binsize in (.1, .2, .3):
    s_hist = functools.partial(histogram, groupdict, binsize)
    resn_combos = (['y'], ['y', 'w'], ['y', 'f'], ['y', 'w', 'f'])
    sele_combos = (['not_removed_res', 'non_ppi_res'],
                   ['not_removed_res', 'non_ppi_res', 'ec_half_res'],
                   ['not_removed_res', 'non_ppi_res', 'pp_half_res'],
                   ['not_removed_res', 'non_ppi_res', 'ec_half_res',
                    'below_18_res'],
                   ['not_removed_res', 'non_ppi_res', 'pp_half_res',
                    'below_18_res'])

    # itertools.product is Cartesian product. The effect is the same as that
    # of nested for loops, the outer over the members of resn_combos and the
    # inner over the members of sele_combos: an operation gets carried out
    # on every pair of a member of resn_combos and a member of sele_combos.
    for resn_combo, sele_combo \
        in itertools.product(resn_combos, sele_combos):
        # Produces a filenamen of the same kind as
        # "yw non_ppi_res pp_half_res aro_peak_res.csv"
        filename = 'structure histograms binsize ' + str(binsize) + '/' \
                   + ''.join(resn_combo) + ' ' \
                   + ' '.join(sele_combo) + '.csv'
        # Create a histogram, then call its save method to save it to that
        # filename:
        s_hist(resn_combo, sele_combo).save(filename)

print('done')
