import math
from sundries import one_letter
from sundries import CIDict

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
    def __init__(self, param_file):
        # The reference dictionary, which maps resns to a dictionary
        # of parameters:
        self.ref = CIDict()

        # A mapping from resn's to column numbers:
        colmap = CIDict()

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
        # curve, and the second column contains the name of the parameter
        # for a gaussian curve. 
        for line in param_file:
            sig_name = line[0]
            gauss_name = line[1]
            for resn, column in colmap.items():
                value = float(line[column])
                if self.ref[resn].curve == 'gaussian':
                    self.ref[resn].add_param(gauss_name, value)
                if self.ref[resn].curve == 'sigmoidal':
                    self.ref[resn].add_param(sig_name, value)

    def calculate(self, resn, z):
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
            return math.exp(-1*(abs(z)-abs(zmin))**2 \
                            / (2*sigma**2))

        if trend.curve == 'sigmoidal':
            e0 = trend.param('e0', sign)
            zmid = trend.param('zmid', sign)
            n = trend.param('n', sign)
            return e0 \
                   / (1 + (abs(z)/abs(zmid))**n)

