"""percent_sasa_and_ss(minsasa,targetss)
percent_sasa_and_ss(.2,'E') and percent_sasa_and_ss(.2,'[^E]?$') and
percent_sasa_and_ss(.2,'.*') are what I usually use"""

import functools
import re
import csv
from pymol import cmd
from pymol import stored

# SELECTORS

def sel_by_percent_sasa_and_ss_generator(spreadsheet, minsasa, targetss):
    for resi, sasa,ss in spreadsheet.get_columns(['resi','%sasa','ss']):
        if sasa != '' and float(sasa) >= float(minsasa) and re.match(targetss,ss) is not None:
            yield 'i. ' + resi

def percent_sasa_and_ss(spreadsheet, minsasa, targetss):
    """percent_sasa_and_ss(minsasa,targetss)
    percent_sasa_and_ss(.2,'E') and percent_sasa_and_ss(.2,'[^E]?$') are
    what I usually use
    Returns a selection string"""
    return '(' + '|'.join(sel_by_percent_sasa_and_ss_generator(spreadsheet,
                                                               minsasa,
                                                               targetss)) + ')'
def from_file(selectionname, moleculename, csvfile):
    '''Create selections from a csv file, formatted like
    1A0S,31,32,35,...
    1AF6,1,2,3...
    
    ARGUMENTS:
    The name you want for your selection (if you want '1AF6.blar', this should
        be 'blar'
    The name of the molecule part of the group from which you're selecting (if
        it's '1AF6.blar', this should be 'blar')
    A csv reader

    RETURNS:
    Nothing'''

    for row in csvfile:
        groupname = row[0]
    
        try:
            end = row[1:].index('')
        except ValueError:
            end = len(row[1:])
        selection = row[1:end]

        selectionstring = '(' +\
                          '|'.join(['i. ' + resi for resi in selection]) \
                          + ')'
        
        cmd.select(groupname + '.' + selectionname,
                   groupname + '.' + moleculename + ' & ' + selectionstring)
                   
    

# OTHER


def append_selection(name, list_, csvfile):
    csvfile.writerow([name] + list(list_))

def write_selection_file(groupdict, target_file,
                         selection_name = '.exposed & n. CA'):
    '''Writes a selection to a csv file, like
    1A0S,31,32,35,...
    1AF6,1,2,3...
    
    ARGUMENTS:
    a groupdict
    a csv writer
    selection_name = '.exposed & n. CA', the selection to write
    
    RETURNS:
    Nothing'''
    file_ = open(target_file, 'wb')
    csvwriter = csv.writer(file_)
    for name, group in groupdict.iteritems():
        stored.list = []
        cmd.iterate(name + selection_name, 'stored.list.append(resi)')
        append_selection(name, stored.list, csvwriter)
    file_.close()       



    
class SequenceMismatch(Exception):
    pass
    
def add_if_included(resi):
    if resi in stored.included_resis:
        stored.included_seq_positions.append(stored.count)
        
def get_seq_positions(object, selection):
    
    '''
    get_seq_positions(object, selection)
    Returns a list of the indices of the full sequence that are included
    in a given selection
    '''.replace('  ', '')
    selection += ' & ' + object
    refinement = ' & polymer & n. CA'
    selection += refinement
    object += refinement
    
    # Make list of included resids
    stored.included_resis = set()
    cmd.iterate(selection,
                'stored.included_resis.add(resi)')
    
    # Make list of included seq positions,
    # where a seq position is an index beginning at 0
    stored.count = 0
    stored.included_seq_positions = list()
    stored.function = add_if_included
    cmd.iterate(object,
                'stored.function(resi); stored.count += 1')
    
    # Check that this can be used to generate the sequence of the selection
    stored.full_sequence = ''
    cmd.iterate(object,
                'stored.full_sequence += one_letter[resn]')
    
    stored.sequence_obtained_directly = ''
    cmd.iterate(selection,
                'stored.sequence_obtained_directly += one_letter[resn]')
                
    sequence_by_positions = ''
    for i in stored.included_seq_positions:
        sequence_by_positions += stored.full_sequence[i]
        
    if sequence_by_positions != stored.sequence_obtained_directly:
        print('Sequence by positions: ')
        print(sequence_by_positions)
        print('Sequence obtained directly: ')
        print(stored.sequence_obtained_directly)
        raise SequenceMismatch('Sequence obtained using included seq ' + \
                               'positions does not match sequence of ' + \
                               'selection')
                               
    return stored.included_seq_positions

def write_seq_positions(groupdict, selection, csvfile, object_name = 'molecule'):
    for name in groupdict:
        # Get the selection:
        seq_positions = get_seq_positions(name + '.' + object_name,
                                          name + '.' + selection)
        
        # The first position of each line holds the pdbid
        to_be_written = [name]

        # The second position holds the sequence of the selection
        # Programs loading these selections will see if they can replicate
        # the sequence. If they can't, they may be using different indices or
        # something else bad.
        stored.sequence = ''
        cmd.iterate(name + '.' + object_name + ' & polymer & n. CA',
                    'stored.sequence += one_letter[resn]')
        filtered_sequence = ''
        for n in seq_positions:
            filtered_sequence += stored.sequence[n]
        to_be_written.append(filtered_sequence)
        
        # All the positions after the second hold indices.
        to_be_written += seq_positions

        csvfile.writerow(to_be_written)
