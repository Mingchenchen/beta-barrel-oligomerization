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
        group.non_ppi_data = biodata.Spreadsheet('non ppi residues/' \
                                                 + csv_name,
                                        phrasebook = weights_phrasebook)
    if csv_name in os.listdir('ppi residues'):
        group.ppi_data = biodata.Spreadsheet('ppi residues/' + csv_name,
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
# Other information will come from the spreadsheet


# Selection of residues in the non_ppi portion of the asymmetric dataset:
# This is SUPERSLOW.
def in_non_ppi_dataset(group, residue):
    ppi_resis = set()
    for entry in group.non_ppi_data.get_column('resi'):
        if entry != '':
            ppi_resis.add(int(entry))
    return residue.get_id()[1] in ppi_resis

making_selections_start= time.time()
for group in groupdict.values():
    group.non_ppi_res = group.selection(in_non_ppi_dataset)

print('making selections took ' + str(time.time() - making_selections_start))


with open('selections/non ppi.csv', 'wb') as f:
    csv.writer(f).writerows([selection_file_line(group, group.non_ppi_res)\
                             for group in groupdict.values()])

