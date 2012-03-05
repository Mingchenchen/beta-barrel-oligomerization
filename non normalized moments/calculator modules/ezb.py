import math
import numpy as np
from Bio.PDB import PDBParser
from useful import CIDict
from useful import one_letter

class NoParameters(Exception):
    pass

class BeCareful(Exception):
    pass

class Calculator(object):
    '''
    Carries out ez-beta calculations using a set of parameters given to
    it at initialization.
    
    Initializing:
    The 'normalize' option must be set to True or False. Be careful!

    The set of parameters must be a spreadsheet represented as
    a list of lists, with the inner lists representing rows. The first row
    must contain the one-letter codes of each amino acid for which
    parameters
    are to be given. Underneath each letter is a column containing its
    parameters in this order:
    Curve type ('gaussian' or 'sigmoidal')
    E0/Emin
    Zmid/Zmin
    n/sigma

    Calculating pseudo-energies:
    calculate(self, resn, z): gives pseudoenergy given a one-letter or
    three-letter code for an amino acid, and a z coordinate
    '''.replace('  ', '')
    def __init__(self, iterable, normalize = None):

        if normalize is None:
            raise BeCareful("Should this calculator normalize results?")
        self.normalize = normalize

        self.ref = dict()
        colmap = dict()
        for column, letter in enumerate(iterable.next()):
            if letter != '':
                self.ref.update({letter: dict()})
                colmap.update({letter: column})
        curvetypes = iterable.next()
        for letter, column in colmap.items():
            self.ref[letter].update({'curve': curvetypes[column]})
        for parameter in [{'sigmoidal': 'e0', 'gaussian': 'emin'},
                          {'sigmoidal': 'zmid', 'gaussian': 'zmin'},
                          {'sigmoidal': 'n', 'gaussian': 'sigma'}]:
            paramrow = iterable.next()
            for letter in self.ref.keys():
                curvetype = self.ref[letter]['curve']
                self.ref[letter].update({parameter[curvetype]: \
                                         float(paramrow[colmap[letter]])})

    def calculate(self, resn, z):
        '''
        calculate(self, resn, z, normalize = None):
        gives pseudoenergy given a one-letter or
        three-letter code for an amino acid, and a z coordinate
        
        Raises a NoParameters exception when you use an amino acid that
        it doesn't have parameters for; always have some way of handling
        this when you call this method!
        '''.replace('  ', '')

        if len(resn) == 3:
            resn = one_letter[resn]        

        if resn not in self.ref:
            raise NoParameters('No parameters for resn ' + str(resn))

        params = self.ref[resn]
        
        if params['curve'] == 'gaussian':
            output = params['emin'] * \
                     math.exp(-1*(abs(z)-params['zmin'])**2 \
                                 /(2*params['sigma']**2))
        elif params['curve'] == 'sigmoidal':
            output = params['e0']/(1+(abs(z)/params['zmid'])**params['n'])
        
        if self.normalize:
            if params['curve'] == 'gaussian':
                output /= params['emin']

            if params['curve'] == 'sigmoidal':
                output /= params['e0']

            # Normalized trends are high energy in middle of the membrane
            # for sigmoidal, high energy in the head-group region for
            # gaussian. For aromatics and small hydrophobics (anything with
            # negative E0 or Emin) these trends should be reversed

            if resn in ('F','Y','W','A','L','I','V'):
                output = 1 - output


        return output

class MissingPDBSequence(Exception):
    pass

class AlignmentOracle(object):
    '''
    Initialized with a multiple sequence alignment and a name to find a
    particular sequence of interest. Can give amino acids from other
    sequences in the alignment that correspond to amino acids in the
    sequence
    of interest, ignoring amino acids that are lined up with gaps in the
    sequence of interest.
    '''.replace('  ', '')
    def __init__(self, alignment, pdb_name = 'pdb'):
        self.data = alignment
        self.pdb_index = None
        for index, record in enumerate(self.data):
            if pdb_name in record.id:
                self.pdb_index = index
                break
        if self.pdb_index == None:
            raise MissingPDBSequence('None of the sequences had '+\
            pdb_name + ' in their title')
        self._pdb_sequence = self.data[self.pdb_index].seq

    def pdb_sequence(self, selection = None):
        return self.sequence(self.pdb_index, selection)

    def sequence(self, index, selection = None):
        pos = 0
        for j, letter in enumerate(self.data[index].seq):
            if self._pdb_sequence[j] != '-':
                if selection is None or pos in selection:
                    yield letter.upper()
                pos += 1

def selections_by_index(iterable):
    selections = CIDict()
    test_sequences = CIDict()
    for line in iterable:
        selections.update({line[0]: [int(x) for x in line[2:]]})
        test_sequences.update({line[0]: line[1]})
    return selections, test_sequences

def selections_by_resi(iterable):
    selections = CIDict()
    for line in iterable:
        selections.update({line[0]: [int(x) for x in line[1:]]})
    return selections


def load_centers(iterable):  
    dict_ = CIDict()
    for row in iterable:
        if row[0] != '':
            dict_.update({row[0]:row[1]})
    for key, value in dict_.items():
        # Turns '(1,2,3)' etc, that is, textual representations of vectors,
        # into Vector objects. Will cut off the last digit of the third
        # componant, but I don't care because the third componant will
        # always be 0.0
        dict_[key] = np.array([float(y[:-1]) for y in value[1:].split()])
    return dict_


def moment(structure, selection, center, calculator,
           paramless_option = None, old_style_gly = False):
    # Add an optional 'sequence' parameter
    # Easy enough, instead of residue.get_resname(), you get the
    # next thing in the sequence

    sum_ = np.zeros(3)
    for residue in structure.get_residues():
        if residue.get_id()[1] not in selection:
            continue

        # My selection files just give residue numbers. If two residues
        # have the same number, if there're insertions, then the selection
        # files are ambiguous and I need to change how I do them.
        # For now, anything identified by more than just a residue number
        # is ignored. The reason this is okay is that, as far as I can tell,
        # in Dan's structures that just means it ignores water.
        # If the structures have complex residue ids, then I need different
        # selection files.
        if residue.get_id()[0] != ' ':
            if residue.get_id()[0] != 'W':
                print('Ignored residue with id {0} in structure {1}'\
                       .format(residue.get_id(), structure))
            continue

        resn = one_letter[residue.get_resname()]

        # Vector points from center to Ca
        coordinates = residue.child_dict['CA'].get_coord()
        vector = coordinates - center
        # Take the projection on the xy plane
        vector[2] = 0
        normalized = vector / np.linalg.norm(vector)        

        # Z energy uses coordinate of Cb
        if resn != 'G':
            z = residue.child_dict['CB'].get_coord()[2]
        else:
            # MasterTMBlist used Ca positions of glycines instead of Cb
            if old_style_gly:
                z = residue.child_dict['CA'].get_coord()[2]
            else:
                # Pseudo cbeta, as if it was alanine
                cb_hydrogen = find_glycine_cb_hydrogen(residue,
                                                       bondlength = 1.5)
                z = cb_hydrogen[2]

        try:
            ezb = calculator.calculate(resn,z)
        except NoParameters:
            if paramless_option == '.5':
                ezb = .5
            elif paramless_option == '0':
                ezb = 0
            else:
                raise
        
        sum_ += normalized * ezb

    return sum_

def find_glycine_cb_hydrogen(glycine, bondlength = 1.1):
    '''
    Uses complete_tetrahedron function to estimate the position of
    the cbeta hydrogen of a glycine.
    Arguments:
    A glycine residue
    bondlength = 1.1, optional parameter for when
    {}the length of the Ca-H bond is not 1.1 angstroms
    '''.replace('  ','').format('    ')
    
    ca_abs = glycine.child_dict['CA'].get_coord()
    n_relative = glycine.child_dict['N'].get_coord() - ca_abs
    c_relative = glycine.child_dict['C'].get_coord() - ca_abs
    cb_h_relative = bondlength * complete_tetrahedron(n_relative,
                                                      c_relative)[0]
    cb_h_abs = cb_h_relative + ca_abs

    return cb_h_abs

def complete_tetrahedron(n, c, angles = (2.186, -2.186)):
    '''
    Given two vectors of a regular tetrahedron centered at the origin,
    will return the other two. The lengths will not be correct - only 
    the directions. The lengths will be one.
    Returns a pair of arrays. To determine which of the pair is the vector
    you're looking for:

    If (n, c) are given as parameters, imagine those two vectors pointing
    toward you, from the origin, with n to your left and c to your right.
    The first of the pair returned will be above the plane of n and c,
    and the second of the pair returned will be below it.
    If n and c are the amino nitrogen and the carboxy carbon of an L amino
    acid, the first returned will be the Cb, and the second returned will
    be the hydrogen.

    Optional parameter 'angles' if you don't want a regular tetrahedron - 
    specify the angles from basis vector m1 that the corners you want are.
    '''.replace('  ', '')

    # Create basis vectors for the plane that contains the other two
    # vectors of the tetrahedron
    m1 = (n + c) / 2
    m2 = np.cross(n,c)
    m1 /= np.linalg.norm(m1)
    m2 /= np.linalg.norm(m2)

    output = [math.cos(angle)*m1 + math.sin(angle)*m2 \
              for angle in angles]
    output = [vec / np.linalg.norm(vec) for vec in output]
    
    return tuple(output)
    
