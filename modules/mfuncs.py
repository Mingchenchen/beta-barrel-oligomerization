'''Some functions that map residue dossiers and sequence identifiers to
numbers. These functions can then be used in the calculation of moments.

Check out blank_function if you want to write your own. The idea is that it's
easier to write a function of a residue dossier and a residue type; usually all
you want with the sequence identifier is to get you the residue type. So write
your function f(dos, resn). Then functools.partial(blank_function, f) is your
function of a residue dossier and a sequence identifier. Derivatives of
blank_function will also raise a GapException when given a gap position,
causing such positions to be skipped in the moment calculation.'''

import functools
import random

class GapException(Exception):
    pass

def ez_b(res_dos, seq_id):
    return res_dos.ez_b(seq_id)

def blank_function(spec_func, res_dos, seq_id):
    column = res_dos.family.res_to_pos[res_dos.residue]
    resn = res_dos.family.msa[seq_id][column]
    if resn == '-':
        raise GapException()
    
    return spec_func(res_dos, resn)

def random_prefunction(res_dos, resn):
    return random.random()

random_magnitude = functools.partial(blank_function, random_prefunction)
