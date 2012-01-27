workingdir = r'C:\cygwin\home\alex\beta barrels\pymol'

import biodata
import selectors
from sundries import CIDict

import os
import re
import functools
import csv

import numpy as np

# This way, I can load objects as G.N and a group called G will be created,
# containing an object called G.N. Also, I can refer to *.N to talk about
# equivalent members of every group (to use commands like disable *.pdb)
cmd.set('group_auto_mode',2)

# To make the moments look clearer:
cmd.set('dash_gap',  0)
cmd.set('dash_width', 6)

one_letter = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
    'GLY':'G', 'PRO':'P', 'CYS':'C'}

class PymolObject(object):    
    def _get_name(self):
        return self._name
    def _set_name(self, name):
        if hasattr(self, '_name'):
            cmd.set_name(self.name,name)
        self._name = name
    def _del_name(self):
        del self._name
    name = property(_get_name, _set_name, _del_name)
    

class Group(PymolObject):
    """Used to associate Python objects with a group. For example, PyMOL can
    associate a graphical representation of a moment with a group, but if you
    also want to keep track of the numbers, you should create a group object
    for that group, and give it an attribute representing the moment
    numerically.
    Initialize with the name of the group. If there is no PyMOL object with
    the given name, the initialization method will create a group with that
    name.
    The Group object has a "name" attribute. Changing the "name" attribute
    should also change the name of the corresponding PyMOL object.
    It doesn't, though. I should fix that."""
    def __init__(self, groupname):
        self.name = groupname
        if groupname not in cmd.get_names('objects'):
            cmd.group(groupname)
    

default_colors = {'selected': 'green', 'unselected': 'purpleblue'}

def color_by_selection(target, source, colors = default_colors):
	cmd.color(colors['selected'], '(' + target + ') & (' + source + ')')
	cmd.color(colors['unselected'], '! (' + target + ') & (' + source + ')')

def groups_from_folder(path, filename_templates, groupname_template = "(*)", load = True):
    """Load a bunch of pdb files. Put them in groups with names taken from
    the filename, name each object groupname.molecule. Return a dictionary
    mapping group names to group objects.
    ARGUMENTS:
    Path in which the pdb files are stored
    List of regular expressions that the filenames must match.
    Template for the names of the groups. When creating a group name for a
        file, the group name will be this template with each instance of '(*)'
        replaced by a group from the regular expression, the first instance
        replaced with the first group, etc. Defaults to '(*)'"""
    dictionary = biodata.file_dict(path, filename_templates,
                                   groupname_template)
    for key, value in dictionary.iteritems():
        dictionary[key] = Group(key)
        if load:
            cmd.load(value, key + '.molecule')
            cmd.disable(key)
    return dictionary

def load_spreadsheets(groupdict, path, filename_templates,
                      groupname_template = '(*)', phrasebook = None,
                      attribute_name = 'spreadsheet'):
    """Loads all csv files in a given path as Spreadsheet objects,
    and makes them into attributes of corresponding group objects.
    
    ARGUMENTS:
    Dictionary mapping names to Group objects
    Path in which to find the spreadsheets
    List of regular expressions to match filenames in the given path.
        Make sure there's a group in the regular
        expression that matches the name, so the corresponding Group object
        in the dictionary can be found. Like, if '1AF6' is a key in the
        dictionary, and you have a file '1AF6_data.csv', you could use
        '(.*)_data.csv'.
    groupname_template = '(*)': template for turning
        group in the regular expressions into the actual
        groupname. Instances of '(*)' in this template will be replaced, in
        order, by groups in the regular expression used to match a filename.
        So, if you need to get '1AF6_group' from '1AF6.csv', this would be
        '(*)._group'. Defaults to '(*)'
    phrasebook = None: phrasebook for the spreadsheets (dictionary mapping the
        names you actually use for columnsto the names for columns in the csv
        file).
    attribute_name = 'spreadsheet': name of the attribute to add to the Group
        objects."""
        
    workbook = biodata.workbook(path, filename_templates, groupname_template)
    for name, spreadsheet in workbook.iteritems():
        try:
            spreadsheet.replace_titles(phrasebook)
        except:
            print 'fucked up on ' + name
            raise
    for group in groupdict.itervalues():
        group.__dict__.update({attribute_name: workbook[group.name]})
        
def load_centers(groupdict, path, attribute_name = 'center'):
    """Loads centers from a csv file. The csv file must have its first column
    full of pdb ids and its second column full of vectors rendered
    (c1, c2, c3).
    
    ARGUMENTS:
    groupdict, path, attribute_name = 'center'
    
    RETURNS:
    Nothing"""
    f = open(path, 'rb')
    freader = csv.reader(f)
    dict_ = {}
    for row in freader:
        if row[0] != '':
            dict_.update({row[0]:row[1]})
    for key, value in dict_.iteritems():
        # Turns '(1,2,3)' etc, that is, textual representations of vectors,
        # into Vector objects. Will cut off the last digit of the third
        # componant, but I don't care because the third componant will
        # always be 0.0
        dict_[key] = np.array([float(y[:-1]) for y in value[1:].split()])
    for group_name, group in groupdict.iteritems():
        groupdict[group_name].__dict__.update({attribute_name: dict_[group_name]})


# CREATE_SESSION HELPER FUNCTIONS:

def cs_load_structures(workingdir, load):
    groupdict = groups_from_folder(workingdir + '/structures',
                                   ['aligned_(.*).pdb'], load = load)
    oligomers = {}
    monomers = {}
    olig_names = set(['1A0S','1AF6','1QD6','2J1N','2O4V', '2POR','1E54','3PRN'])
    for name, group in groupdict.iteritems():
        if name in olig_names:
            cmd.group('oligomers',name)
            oligomers.update({name: group})
        else:
            cmd.group('monomers',name)
            monomers.update({name: group})
    return groupdict, oligomers, monomers

def cs_load_centers(workingdir, groupdict):
    for name, path in [('sheets_center',
                        '/centers/beta sheets and %sasa over .2.csv'),
                       ('cored_center','/centers/cored 1.csv')]:
        load_centers(groupdict, workingdir + path, attribute_name = name)

def cs_load_spreadsheets(workingdir, groupdict):
    phrasebooks = biodata.phrasebooks(workingdir + '/misc/phrasebooks.csv')
    ezbeta_phrasebook = phrasebooks['dan_ezbeta']
    load_spreadsheets(groupdict, workingdir + \
                      '/ez beta data',
                      ['(.*)_DeltaBurial.csv$'],phrasebook = ezbeta_phrasebook,
                      attribute_name = 'ez_data')

def cs_make_selections(groupdict):
    for name, group in groupdict.iteritems():
        cmd.select(name + '.exposed', name + ' & ' + \
                   selectors.percent_sasa_and_ss(group.ez_data,.2,'.*'))
        cmd.select(name + '.exp_sheets', name + ' & ' + \
                   selectors.percent_sasa_and_ss(group.ez_data,.2,'E'))
        
    cored_selections = open(workingdir + '/selections/cored 1.csv', 'rb')
    cored_csv = csv.reader(cored_selections)
    selectors.from_file('cored','molecule',cored_csv)
    cored_selections.close()

def cs_fetch_oligomers(oligomers, load):
    if load:
        for name, group in oligomers.iteritems():
            cmd.fetch(name, name+'.pdb')
            cmd.align(name+'.pdb', name+'.molecule')
    
def create_session(workingdir, load = True):
    # Load structures:
    groupdict, oligomers, monomers = cs_load_structures(workingdir, load)

    # In case something goes wrong, so you can look at the work in progress:
    stored.groupdict = groupdict
    
    # Load centers:
    cs_load_centers(workingdir, groupdict)
    
    # Load spreadsheets:
    cs_load_spreadsheets(workingdir, groupdict)
    
    # Make selections:
    cs_make_selections(groupdict)
    
    # Fetch pdbs for oligomers:
    cs_fetch_oligomers(oligomers, load)
    
    # Change from line to cartoon representation
    cmd.hide('lines','*')
    cmd.show('cartoon','*')
    
    return (groupdict, oligomers, monomers)

try:
    groupdict
    print 'already loaded, do groupdict, oligomers, monomers = '+\
          'create_session(workingdir, load=True) to reload'
except NameError:
    groupdict, oligomers, monomers = create_session(workingdir, load=True)
    print 'done. groupdict, monomers and oligomers defined.'
