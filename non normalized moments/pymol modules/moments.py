""" MOMENT MODULE
THE VECTOR CLASS: initialize with a tuple of three numbers.
    I couldn't figure out how to install numpy on PyMOL so I made a
    vector class.
MOMENT FUNCTIONS:
moment(sele, pairs, state=1, only_c_alphas=True)
consurf_moment(sele, grades, state=1, only_c_alphas=True)"""

import csv
import math
from pymol import cmd
from pymol import stored

# THE VECTOR CLASS

class VectorLengthError(Exception):
    # Called when someone forgets that my vector class is only for vectors of length 3.
    # Example of usage: __init__ method of the vector class.
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class Vector(object):
    """ A tuple with some new methods. Immutable and iterable, like tuples.
    
    a * b: returns a dot product if a and b are vectors, scalar product
        if one is a scalar
    a + b: vector addition
    a[n]: get component, if n is 0, 1 or 2 and a is a vector
    norm: returns length of the vector
    normalize: returns normalized vector
    getTuple: returns the tuple containing the vector's components
    """
    
    def __init__(self,data):
        if len(data) > 3:
            raise VectorLengthError('Vectors cannot have more than three components')
        else:
            componentlist = list(data)
            for n in xrange(3 - len(data)):
                componentlist.append(0)
            self.data = tuple(componentlist)
            
    def __iter__(self):
        return iter(self.get_tuple())
        
    def __len__(self):
        return len(self.data)
            
    def __repr__(self):
        return 'Vector(' + repr(self.data) + ')'
    
    def __str__(self):
        return str(self.data)
        
    def __add__(self,other):
        # Vector addition:
        if type(other) != type(self):
            raise TypeError('Vectors can only be added to other vectors')
        return Vector((self.data[0] + other.data[0], self.data[1] + other.data[1], self.data[2] + other.data[2]))
        
    def __mul__(self,other):
        # Dot product
        if type(other) == type(self):
            return (self.data[0] * other.data[0]) + (self.data[1] * other.data[1]) + (self.data[2] * other.data[2])
        # Scalar multiplication
        else:
            return Vector( (self.data[0] * other, self.data[1] * other, self.data[2] * other) )
    
    def __rmul__(self,other):
        # Adds support for right side multiplication        
        return self.__mul__(other)
        
    def __sub__(self,other):
        return self + -1 * other
        
    def __getitem__(self,n):
        # So you can talk about vector[1] etc.
        return self.data[n]
    
    def norm(self):
        return math.sqrt(self.data[0]**2 + self.data[1]**2 + self.data[2]**2)
    
    def normalize(self):
        norm = self.norm()
        if norm == 0:
            return self
        else:
            return (1 / self.norm()) * self
            
    def get_tuple(self):
        """tuple(vector object) also works, but this is more straightforward.
        If it ever actually matters I guess I'll time them and use whatever's
        faster."""
        return self.data

        
    
        
# MOMENT FUNCTIONS

def moment(sele, pairs, center, state=1, only_c_alphas=True):
    """Returns a moment calculated from the resis and values in a given
    dictionary, or list of (resi, value) tuples, with resi as an int
    ARGUMENTS:
    A selection string
    A list of (resi,value) tuples or a dictionary.
        resi must be an int, value can be an int or float.
    A position vector or tuple of the center from which the moment
        will be calculated
    KEYWORD ARGUMENTS:
    state=1, the state from which to retrieve atomic coordinates
    only_c_alphas=True, if true adds " & n. CA" to the selection
    RETURNS:
    A vector"""
    
    if only_c_alphas:
        sele += ' & n. CA'
    
    # Convert resis to strings and scores to floats, and make a
    # dictionary mapping resis to scores
    stringpairs = []
    for resi, score in pairs:
        # Some csv files have empty rows
        if str(resi) != '':
            stringpairs.append((str(resi), float(score)))
    stored.pairsdict = dict(stringpairs)    
    
    # Moments are calculated from projections on the z = 0 plane, which, when
    # I use these functions, is parallel to the membrane
    stored.center = Vector(center)
    
    stored.sum = Vector((0,0,0))    
    cmd.iterate_state(state, sele,
        'stored.sum += stored.pairsdict[resi] * (Vector((x,y,0)) + \
        -1 *stored.center).normalize()')
    return stored.sum
    
    
def consurf_moment(sele, grades, center, state=1, only_c_alphas=True):
    """Returns a consurf moment.
    ARGUMENTS:
    A selection string, a gfile, a center
    KEYWORD ARGUMENTS:  
    state=1, the state from which to retrieve atomic coordinates
    only_c_alphas=True, if true adds " & n. CA" to the selection
    RETURNS:
    A vector"""
    
    # Multiplied by -1 because low ConSurf scores indicate high conservation
    return -1 * moment(sele, grades.certain_score_list(), center,
                       state = state, only_c_alphas = only_c_alphas)
                       
def basis_moment(sele, center, only_c_alphas=True):
    """What you would get if you ran any other moment giving every residue a
    score of one.
    ARGUMENTS:
    A selection string, a center
    KEYWORD ARGUMENTS:
    only_c_alphas=True, if true adds "& n. CA" to the selection
    RETURNS:
    A vector"""
    
    if only_c_alphas:
        sele += " & n. CA"
    
    stored.list = []
    
    cmd.iterate(sele, 'stored.list.append((resi,1))')
    
    return moment(sele, stored.list, center, only_c_alphas = only_c_alphas)

def draw_vector(name, vector_, center):
    """Create a distance with the length and direction of a given vector,
    and then hide its label.
    The center should have the same xy projection as the center you gave to
    the moment function, but you might want to change the z value to make it
    more visible.
    
    ARGUMENTS:
    (name, vector_, center)
    RETURNS:
    Nothing."""
    
    # The awkward tuple(Vector(center)) thing turns (a,b) into (a,b,0)
    cmd.pseudoatom('ps1', pos = tuple(Vector(center)))
    cmd.pseudoatom('ps2', pos = tuple((Vector(center) + Vector(vector_))))
    cmd.distance(name, 'ps1', 'ps2')
    cmd.hide('labels', name)
    cmd.delete('ps1 ps2')