import math
from sundries import one_letter
from sundries import CIDict

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
    '''
    
    def __init__(self, iterable, normalize = None):

        if normalize is None:
            raise BeCareful("Should this calculator normalize results?")
        self.normalize = normalize

        self.ref = CIDict()
        colmap = CIDict()
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
        gives pseudoenergy given a one-letter or
        three-letter code for an amino acid, and a z coordinate
        
        Raises a NoParameters exception when you use an amino acid that
        it doesn't have parameters for; always have some way of handling
        this when you call this method!
        '''

        if len(resn) == 3:
            resn = one_letter[resn]        

        try:
            params = self.ref[resn]
        except KeyError:
            raise NoParameters('No parameters for resn ' + str(resn))
        
        if params['curve'] == 'gaussian':
            output = params['emin'] * \
                     math.exp(-1*(abs(z)-params['zmin'])**2 \
                                 /(2*params['sigma']**2))
        elif params['curve'] == 'sigmoidal':
            output = params['e0']/(1+(abs(z)/params['zmid'])**params['n'])
        
        if self.normalize:
            if params['curve'] == 'gaussian':
                output /= params['emin']
                # Normalized trends are high energy in middle of the membrane
                # for sigmoidal, high energy in the head-group region for
                # gaussian. For aromatics and small hydrophobics (anything with
                # negative E0 or Emin) these trends should be reversed
                if params['emin'] < 0:
                    output = 1 - output

            if params['curve'] == 'sigmoidal':
                output /= params['e0']
                if params['e0'] < 0:
                    output = 1 - output

        return output

