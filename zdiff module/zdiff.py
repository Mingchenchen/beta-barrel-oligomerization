from sundries import CIDict
from Bio.PDB import PDBParser
import warnings
from Bio import AlignIO

# The guts that it runs on
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
        '''Return resi of position in specified structure'''

        # Get the residue of the specified structure that is in this
        # position in the alignment.
        target = self.residues[stru_id]

        # It might actually be a gap; return None of it is
        if type(target) is Gap:
            return None

        # Otherwise, return the id of the residue. It's a Biopython
        # Residue object, so this is done with its get_id() method.
        else:
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

    def resi_report(self, template_id, unknown_id):
        output = ((pos.resi(template_id),
                   pos.zdiff(template_id, unknown_id)) \
                  for pos in self.positions)
        return filter(lambda x: x[1] is not None, output)
        

# The API for making zdiff files

def calc(template_name, target_name, template_structure_filename,
         target_structure_filename, alignment_filename,
         write_to, comment, format_='clustal'):

    # Open relevant files
    with open(template_structure_filename, 'r') as template_structure_file,\
         open(target_structure_filename, 'r') as target_structure_file,\
         open(alignment_filename, 'r') as alignment_file:
        # Load structures with Biopython's PDB file parser
        # Daniel's aligned structures are missing some inessential
        # information, and as a consequence the parser gives thousands
        # of warnings. Gotta ignore these.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            template_structure = PDBParser().\
                                 get_structure(template_name,
                                               template_structure_file)
            target_structure = PDBParser().\
                               get_structure(target_name,
                                             target_structure_file)
        
        # Open alignment using Biopython's parser
        alignment = AlignIO.read(alignment_filename, format_)

    # Find the template and target sequences in the alignment
    for seq_record in alignment:
        if seq_record.id == template_name:
            templ_seq = seq_record
        if seq_record.id == target_name:
            targ_seq = seq_record
    
    # Calculate zdiff
    results = Zdiff((template_structure, templ_seq),
                    (target_structure, targ_seq))
                    
    # Write results to a file
    with open(write_to, 'w') as o:
        # Write some coments so I know which zdiff file this is
        o.write('# Template: ')
        o.write(template_name + ' (' + template_structure_filename + ')\n')
        o.write('# Target: ')
        o.write(target_name + ' (' + target_structure_filename + ')\n')
        o.write('# Alignment: ' + alignment_filename + '\n')
        o.write('# ' + comment + '\n')

        # Write the actual data
        # Weird quirk of the zdiff objects - before giving a report, it
        # requires the id's of the structures. I don't know if I wrote
        # it to support more than two, or if I was planning to, or what.
        for resi, zdiff in results.resi_report(template_name, target_name):
            o.write(str(resi) + ', ' + str(zdiff) + '\n')


