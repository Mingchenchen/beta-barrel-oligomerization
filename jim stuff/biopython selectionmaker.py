import re
from Bio.PDB import PDBParser
from sundries import CIDict
import functools
import csv
import os
import warnings
import biodata
import time

start = time.time()

class Group(object):
    def __init__(self, name):
        self.name = name

    def selection(self, function):
        return filter(functools.partial(function, self),
                      self.structure.get_residues())
    
# Retrieve all the pdbids that are in the asymmetric dataset
filenames = os.listdir('non ppi residues')
matches = [re.match('(\d...)\.csv', x) for x in filenames]
real_matches = filter(lambda x: x is not None, matches)
pdbids = [x.group(1) for x in real_matches]

# Uncomment this line to speed things up for testing:
# pdbids = ['1QD6']

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


# Selection of residues in the |z| in [6,15] range, where Daniel's
# asymmetric EzB charts show peaks of aromatic residue frequency

def in_aromatic_peak(group, residue):
    try:
        z = residue.child_dict['CA'].get_coord()[2]
    except KeyError:
        return False
    return abs(z) > 6 and abs(z) < 18

for group in groupdict.values():
    group.aro_peak_res = group.selection(in_aromatic_peak)

with open('selections/aromatic peak.csv', 'wb') as f:
    
    for group in groupdict.values():
        row_to_write = list()
        row_to_write.append(group.name)
        row_to_write += [residue.get_id()[1] \
                         for residue in group.aro_peak_res]
        csv.writer(f).writerow(row_to_write)
