'''Rough script to load Daniel's selections. If I use this a lot, I'll refine
it and extend it ot load PPI selections as well.
Maybe make a context manager instead of that awkard try/finally?'''

import re
import os

import biodata
from sundries import CIDict

try:
    # Save the current working directory so you can restore it at end:
    previous_wd = os.getcwd()
    os.chdir(r'C:\cygwin\home\alex\beta barrels\jim stuff')

    # Get the phrasebook for interacting with Daniel's "weights" files
    phrasebooks = biodata.phrasebooks('weights phrasebook.csv')

    # Retrieve names of csv files containing the resis
    # is there something I can use that's simpler than os.walk,
    # and just returns all the filenames in a given directory?
    folder = 'non ppi residues'
    filename_superset = os.listdir(folder)

    def filename_match(filename):
        return re.match('\d...\.csv', filename) is not None

    weight_file_filenames = filter(filename_match, filename_superset)
    weight_file_paths = [folder + '/' + filename \
                         for filename in weight_file_filenames]


    # Make spreadsheets
    spreadsheets = [biodata.Spreadsheet(filename,
                                        phrasebook = phrasebooks['weights'])\
                    for filename in weight_file_paths]

    # weights maps pdbids to spreadsheets
    weights = CIDict()
    for spreadsheet in spreadsheets:
        pdbid_list = spreadsheet.get_column('pdbid')
        pdbid_filtered = filter(lambda x: x != '', pdbid_list)
        pdbid_set = set(pdbid_filtered)
        # I expect there to only be one pdbid
        assert len(pdbid_set) == 1, 'more than 1 pdbid in one spreadsheet'
        pdbid = list(pdbid_set)[0]
        weights.update({pdbid: spreadsheet})

    # selections maps pdbids to sets of resis
    selections = CIDict()
    for pdbid, spreadsheet in weights.items():
        def not_blank(string):
            return string != ''
        resis = filter(not_blank, spreadsheet.get_column('resi'))
        selections.update({pdbid: set(resis)})

    # A new global variable for looping over these proteins
    asymmetric_dataset = CIDict([(pdbid, groupdict[pdbid]) \
                                 for pdbid in weights.keys()])

    # Make the spreadsheets available through groupdict
    for pdbid, group in asymmetric_dataset.items():
        group.non_ppi = weights[pdbid]


    # This colorizing script is REALLY SLOW
    for pdbid, selection in selections.items():
        cmd.select(pdbid.upper() + '.non_ppi', 'none')
        for resi in selection:
            cmd.select(pdbid.upper() + '.non_ppi',
                       '{0}.molecule & i. {1} | {0}.non_ppi' \
                       .format(pdbid.upper(),resi))

finally:
    os.chdir(previous_wd)

def ss(group, target_resi):
    '''ss(group,resi)
    retrieves secondary structure using non ppi dataset spreadsheet'''
    for resi, ss in group.non_ppi.get_columns(['resi','ss']):
        if resi != '' and int(resi)==int(target_resi):
            return ss
