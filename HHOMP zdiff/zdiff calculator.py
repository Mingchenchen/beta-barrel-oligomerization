from sundries import CIDict
from Bio.PDB import PDBParser
import warnings
from Bio import AlignIO

def z(residue):
    '''Returns the z coordinate of a residue object's Calpha.'''
    return residue.child_dict['CA'].get_coord()[2]

class Gap(object):
    '''Represents a gap in a Position object.'''
    pass

class Position(object):
    def __init__(self, pairs):
        self.residues = CIDict(pairs)

    def zdiff(self, template_id, unknown_id):
        template_res = self.residues[template_id]
        unknown_res = self.residues[unknown_id]
        if type(template_res) is Gap or type(unknown_res) is Gap:
            return None

        else:
            return z(template_res) - z(unknown_res)
    
    def resi(self, stru_id):
        return self.residues[stru_id].get_id()[1]
        
class NotFoundError(Exception):
    pass

class Zdiff(object):
    def __init__(self, *stru_seq_pairlist):
        pos_inputs = list()
        for structure, sequence in stru_seq_pairlist:
            residues = structure.get_residues()
            pairs_for_pos = list()

            for letter in sequence:
                if letter == '-':
                    seq_unit = Gap()
                else:
                    seq_unit = residues.next()

                pairs_for_pos.append((structure.get_id(), seq_unit))

            pos_inputs.append(pairs_for_pos)
       
        self.positions = [Position(id_res_pairlist) \
                          for id_res_pairlist in zip(*pos_inputs)]
        
    def get(self, template_id, unknown_id, resi, start=0):
        for pos in self.positions[start:]:
            if pos.resi(unknown_id) == resi:
                return pos.zdiff(template_id, unknown_id)
        
        raise NotFoundError('resi {} of {} not found' \
                            .format(resi, unknown_id))
    
    def report(self, template_id, unknown_id):
        output = (pos.zdiff(template_id, unknown_id) \
                  for pos in self.positions)
        return filter(lambda x: x is not None, output)

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    scry_stru = PDBParser().get_structure('1A0S', 'aligned_1A0S.pdb')
    maltoporin_stru = PDBParser().get_structure('1AF6', 'aligned_1AF6.pdb')

alignment = AlignIO.read('Swiss-PDB structural alignment.aln', 'clustal')
for seq_record in alignment:
    if seq_record.id == 'aligned_1A0S':
        scry_seq = seq_record
    if seq_record.id == 'aligned_1AF6':
        maltoporin_seq = seq_record

output = Zdiff((scry_stru, scry_seq), (maltoporin_stru, maltoporin_seq))

with open('NEW 1a0s as template.txt', 'w') as f:
    f.writelines(str(zdiff) + '\n' \
                 for zdiff in output.report('1A0S', '1AF6'))

with open('NEW 1af6 as template.txt', 'w') as f:
    f.writelines(str(zdiff) + '\n' \
                 for zdiff in output.report('1af6', '1a0s'))
