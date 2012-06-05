workingdir = r'C:\cygwin\home\alex\beta barrels\pymol'

import biodata
import selectors
from sundries import CIDict

import os
import re
import functools
import csv

import numpy as np

os.chdir(r'C:\cygwin\home\alex\beta barrels\jim stuff')

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
    
# CREATE_SESSION HELPER FUNCTIONS:
    
def create_session(workingdir, load = True):
    # Load structures:
    groupdict = CIDict(groups_from_folder(workingdir + '/structures',
                                   ['aligned_(.*).pdb'], load = load))
                                   
    # Remove structures not in the datset
    included_proteins = set()
    for filename in os.listdir('non ppi residues'):
        match = re.match('(\d...)\.csv', filename) 
        if match is not None:
            included_proteins.add(match.group(1))
    
    for pdbid in groupdict.keys():
        if pdbid.upper() not in included_proteins:
            del groupdict[pdbid]
            cmd.delete(pdbid)

    # In case something goes wrong, so you can look at the work in progress:
    stored.groupdict = groupdict
    
    
    
    # Change from line to cartoon representation
    cmd.hide('lines','*')
    cmd.show('cartoon','*')
    
    
    
    return groupdict

try:
    groupdict
    print 'already loaded, do groupdict=create_session(workdingdir) to reload'
except NameError:
    groupdict = create_session(workingdir, load=True)
    print 'done. groupdict, monomers and oligomers defined.'
