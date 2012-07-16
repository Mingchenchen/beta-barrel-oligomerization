# So I can run on my stupid home computer:
#import sys
#sys.path.append('/home/davis/work stuff/beta-barrel-oligomerization/modules')

import matplotlib.pyplot as plt
import math
import csv
from sundries import one_letter
from sundries import CIDict
import re
import os

class NoParameters(Exception):
    pass


class ResTrend(object):
    def __init__(self, curvetype):
        self.curve = curvetype
        self.params = dict()

    def param(self, name, sign):
        return self.params[name + sign]

    def add_param(self, name, value):
        self.params.update({name: value})

class Calculator(object):
    def __init__(self, param_file_path):
        # The reference dictionary, which maps resns to a dictionary
        # of parameters:
        self.ref = CIDict()

        # A mapping from resn's to column numbers:
        colmap = CIDict()
        
        with open(param_file_path, 'rb') as f:
            param_file = csv.reader(f)
            # The first line contains the resn's
            resn_line = param_file.next()
            for col_num, letter in enumerate(resn_line):
                if letter != '':
                    colmap.update({letter: col_num})
    
            # The next line contains the types of curve
            # Create a ResTrend object corresponding to each residue
            curvetypes = param_file.next()
            for resn, column in colmap.items():
                self.ref.update({resn: ResTrend(curvetypes[column])})
    
            # The rest of the lines contain the parameters. The first
            # column contains the name of the parameter for a sigmoidal
            # curve, and the second column contains the name of the 
            # parameter for a gaussian curve. 
            # It's completely stupid but all of the parameter names
            # are lowercase by convention except for "Paq". That should
            # be changed someday.
            for line in param_file:
                sig_name = line[0]
                gauss_name = line[1]
                for resn, column in colmap.items():
                    value = float(line[column])
                    if self.ref[resn].curve == 'gaussian':
                        self.ref[resn].add_param(gauss_name, value)
                    if self.ref[resn].curve == 'sigmoidal':
                        self.ref[resn].add_param(sig_name, value)
    
    def dez(self, resn, z):
        if len(resn) == 3:
            resn = one_letter[resn]   

        try:
            trend = self.ref[resn]
        except KeyError:
            raise NoParameters('No parameters for resn ' + str(resn))

        if z>0:
            sign = '+'
        else:
            sign = '-'

        if trend.curve == 'gaussian':
            emin = trend.param('emin', sign)
            zmin = trend.param('zmin', sign)
            sigma = trend.param('sigma', sign)
            return emin * math.exp(-1*(abs(z)-abs(zmin))**2 \
                                   / (2*sigma**2))

        if trend.curve == 'sigmoidal':
            e0 = trend.param('e0', sign)
            zmid = trend.param('zmid', sign)
            n = trend.param('n', sign)
            return e0 \
                   / (1 + (abs(z)/abs(zmid))**n)

    def negative_log_propensity(self, resn, z):
        if len(resn) == 3:
            resn = one_letter[resn]

        if z > 0:
            sign = '+'
        else:
            sign = '-'

        E_z = self.dez(resn, z)
        Paq = self.ref[resn].param('Paq', sign)

        # It's NEGATIVE log that's like energy
        return -1 * math.log(Paq * math.exp(-1 * E_z / (0.00198717*298)))

    def normalized_dez(self, resn, z):
        if len(resn) == 3:
            resn = one_letter[resn]

        if z > 0:
            sign = '+'
        else:
            sign = '-'
        trend = self.ref[resn]

        if trend.curve == 'gaussian':
            amplitude = trend.param('emin',sign)
        if trend.curve == 'sigmoidal':
            amplitude = trend.param('e0', sign)
        
        deamplituded = self.dez(resn, z) / amplitude

        if amplitude > 0:
            normalized = deamplituded
        else:
            normalized = 1 - deamplituded

        return normalized
           
    def energy_plot(self, resn):
        '''Displays a plot of ez-beta for a given residue type (one letter)
        from -24 to 24 angstroms'''
        plt.plot(range(-24,24),
                 [self.dez(resn, z) for z in range(-24,24)])
        plt.xlabel('distance from membrane center (anstroms)')
        plt.ylabel('ez-beta')
        plt.ylim(-3,3)
        plt.show()

    def test(self, resn):
        '''Displays a plot of propensity for a given residue type (one letter)
        from -24 to 24 angstroms, for comparison against the plots that
        Daniel makes'''
        def propensity(z, E_z):
            if z > 0:
                sign = '+'
            else:
                sign = '-'
            Paq = self.ref[resn].param('Paq', sign)
            return Paq * math.exp(-1 * E_z / (0.00198717*298))
                
        plt.plot(range(-18,18),
                 [propensity(z, self.dez(resn, z)) \
                  for z in range(-18,18)])
        plt.xlabel('distance from membrane center (anstroms)')
        plt.ylabel('ez-beta')
        plt.ylim(0,3)
        plt.show()

def color(function, selection, params = 'bbtm_out', zero_is_white = True):
    '''Do not try to color more than one protein at once. It'll
    come out wrong.
    Beware: when using normalized dEz, "zero_is_white" actually means
    ".5 is white". In a future version maybe I'll change this to
    "neutral is white" or something.
    '''

    if params == 'bbtm_out':
        param_path = 'asymmetric params, bbtmout.csv'
        # Residues not included in these params, to be excluded:
        excl = ['cys']
        # Threonine's flat distribution leads to arbitrary energies if
        # you're normalizing
        if function == 'normalized_dez':
            excl.append('thr')

    if params == 'symmetric':
        param_path = "published params, asymmetric format.csv"
        excl = ('cys', 'met', 'thr')

    stored.calc = Calculator(param_path)
    stored.resi_z = dict()
    cmd.iterate_state(1, selection + ' & n. ca & !resn ' + '+'.join(excl),
                      'stored.resi_z.update({resi: float(z)})')
    stored.assigned_vals = list()
    cmd.alter(selection + ' & polymer & !resn ' + '+'.join(excl),
              'b = stored.calc.{0}(resn, stored.resi_z[resi])' \
              .format(function) + ';stored.assigned_vals.append(b)')

    # I want the spectrum to be of a range inclusive of whatever
    # values are there, but I want 0 to be in the middle.
    # Unless it's normalized, then I want .5 to be in the middle.
    if function == "normalized_dez":
        min_color = 0
        max_color = 1
    else:
        min_b = min(stored.assigned_vals)
        max_b = max(stored.assigned_vals)
        limit = max((abs(min_b), abs(max_b)))
        min_color = -1 * limit
        max_color = 1 * limit
    
    if not zero_is_white:
        cmd.spectrum('b','blue_white_red',
                     selection + ' & !resn ' + '+'.join(excl))
    # The way to test whether this is performing as expected is to
    # run without zero_is_white, and then with. Either the deepest red
    # or the deepest blue should change, but not both. So there is always
    # one residue whose color does not change: the one whose b factor
    # has the highest absolute value. I have not performed this
    # test.
    if zero_is_white:
         cmd.spectrum('b','blue_white_red',
                     selection + ' & !resn ' + '+'.join(excl),
                     minimum = min_color, maximum = max_color)
    
    cmd.color('black', selection + ' & resn ' + '+'.join(excl))

cmd.extend('ezb_color', color)
