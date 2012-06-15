from sundries import CIDict
from sundries import one_letter
from Bio.PDB import PDBParser
import warnings
import re
import os
import csv

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
        if type(self.residues[stru_id]) is Gap:
            return None
        else:
            return self.residues[stru_id].get_id()[1]
        
class NotFoundError(Exception):
    pass

class Zdiff(object):
    def __init__(self, *stru_seq_pairlist):
        pos_inputs = list()
        for structure, sequence, beginning in stru_seq_pairlist:
            residues = structure.get_residues()
            pairs_for_pos = list()

            for letter in sequence:
                if letter == '-':
                    seq_unit = Gap()
                else:
                    seq_unit = residues.next()
                    while 'CA' not in seq_unit.child_dict.keys():
                        seq_unit = residues.next()
                    while int(seq_unit.get_id()[1]) < beginning:
                        seq_unit = residues.next()
                    assert one_letter[seq_unit.get_resname()] == letter, \
                           str(seq_unit.get_id())+ str(seq_unit.child_dict)\
                           + structure.id

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


    def resi_report(self, template_id, unknown_id):
        output = ((pos.resi(template_id),
                   pos.zdiff(template_id, unknown_id)) \
                  for pos in self.positions)
        return filter(lambda x: x[1] is not None, output)


linerec = re.compile('Chain (\d): +(\d+) ([A-Z\-]+)')
for alignment_filename in os.listdir('pairwise alignments'):
    ids = re.match('(\d...)_\.(\d...)_\.txt', alignment_filename)
    id1 = ids.group(1)
    id2 = ids.group(2)
    seq1 = ''
    seq2 = ''
    beginning1 = -1
    beginning2 = -1
    with open('pairwise alignments/' + alignment_filename, 'r') as f:
        for line in f:
            x = re.match(linerec, line)
            if x is not None:
                if x.group(1) == '1':
                    if beginning1 == -1:
                        beginning1 = int(x.group(2))
                    seq1 += x.group(3)
                if x.group(1) == '2':
                    if beginning2 == -1:
                        beginning2 = int(x.group(2))
                    seq2 += x.group(3)
    
    for filename in os.listdir('aligned structures'):
        stru_match = re.match("(\d...) to daniel's \d...\.pdb", filename)
        if stru_match.group(1) == id1:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                structure1 = PDBParser() \
                             .get_structure(id1,
                                            'aligned structures/' \
                                            + filename)
        if stru_match.group(1) == id2:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                structure2 = PDBParser() \
                             .get_structure(id2,
                                           'aligned structures/' \
                                           + filename)
    ''' 
    for structure, beginning in ((structure1, beginning1),
                                 (structure2, beginning2)):
        for model in structure.get_iterator():
            for chain in model.get_iterator():
                kill_list = list()
                print('beginning is ' + str(beginning))
                for full_id, residue in chain.child_dict.items():
                    if int(residue.get_id()[1]) < beginning:
                        kill_list.append(full_id)
                    elif 'CA' not in residue.child_dict.keys():
                        kill_list.append(full_id)
                for key in kill_list:
                    del chain.child_dict[key]
                    print('eliminated ' + str(key))
                    assert key not in chain.child_dict.keys()
    '''
    output = Zdiff((structure1, seq1, beginning1), (structure2, seq2, beginning2))
    write_filename = id1 + ' mapped on ' + id2 + '.csv'
    with open('zdiffs/' + write_filename, 'wb') as f:
        csv.writer(f).writerows(output.resi_report(id2, id1))
