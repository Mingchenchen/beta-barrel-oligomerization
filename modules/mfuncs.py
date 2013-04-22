'''Some functions that map residue dossiers and sequence identifiers to
numbers. These functions can then be used in the calculation of moments.

Check out blank_function if you want to write your own. The idea is that it's
easier to write a function of a residue dossier and a residue type; usually all
you want with the sequence identifier is to get you the residue type. So write
your function f(dos, resn). Then functools.partial(blank_function, f) is your
function of a residue dossier and a sequence identifier.'''

import functools
import random

def ez_b(res_dos, seq_id):
    return res_dos.ez_b(seq_id)

def random(res_dos, seq_id):
    return random.random()

def blank_function(spec_func, res_dos, seq_id):
    column = self.family.res_to_pos[self.residue]
    resn = self.family.msa[seq_id][column]
    
    return spec_function(res_dos, resn)

