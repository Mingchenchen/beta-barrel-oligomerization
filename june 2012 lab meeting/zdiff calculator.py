from sundries import CIDict

def z(residue):
    '''Returns the z coordinate of a residue object's Calpha.'''
    return residue.child_dict['CA'].get_coord()[2]

class Gap(object)
    '''Represents a gap in a Position object.'''
    pass

class Position(object):
    def __init__(self, *pairs):
        self.residues = CIDict(pairs)

    def zdiff(self, template_id, unknown_id):
        template_res = self.residues[template_id]
        unknown_res = self.residues[unknown_id]
        if type(template_res) is Gap or type(unknown_res) is Gap:
            return None

        else:
            return z(template_res) - z(unknown_res)
    
    def resi(self, stru_id):
        

class Zdiff(object):
    def __init__(self, alignment, 
