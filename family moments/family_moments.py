from sundries import one_letter
import warnings
from Bio.PDB import PDBParser
import Bio.AlignIO
import DSSP_win
import zenergy
import csv

def map_res_to_pos(residues, sequence):
    '''Given a list of residues from a Biopython structure and a sequence
    of the same protein from a sequence alignment,
    return a mapping from the residues to their positions in the sequence.

    Unfortunately this function depends critically on a flaw in my
    method of deriving moments. PDB structures often have residues missing,
    but in its current implementation this function will fail unless
    those residues are missing from the sequence as well. Due to a mistake
    generating the alignments, the resiudes ARE, in fact, missing from the
    sequences, so it works. But the alignments may be of lower quality
    due to that mistake.'''
    # From the sequence, make an ordered list of pairs (index, resname)
    # excluding gapped positions
    numbered_resnames = [pair for pair in enumerate(sequence) \
                         if pair[1] != '-']
    
    # Make a mapping from residues to indices of the sequence
    res_to_index = dict((res, pair[0]) \
                        for res, pair in zip(residues, numbered_resnames))

    # Check that the residues actually match
    for res, index in res_to_index.items():
        if one_letter[res.get_resname()] != sequence[index]:
            raise ValueError("List of residues doesn't match sequence")

    return res_to_index

class Structure(object):
    def __init__(self, stru_name, stru_path, msa_path, template_name,
                 param_path):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.stru = PDBParser().get_structure(stru_name, stru_path)
        msa = Bio.AlignIO.read(msa_path, 'clustal')
        self.msa = dict((seq.id, seq) for seq in msa)
        self.template_seq = self.msa[template_name]

        self.res_to_pos = map_res_to_pos(self.stru.get_residues(),
                                         self.template_seq)
        
        self.dssp = DSSP_win.DSSP(self.stru.child_dict[0], stru_path)

        params = csv.reader(open(param_path, 'rb'))
        self.calc = zenergy.Calculator(params)

# Tests
x = Structure('1A0S', 'structures/aligned_1A0S.pdb',
              'gonnet aligned/1A0S with cluster73.clu', 'template_1A0S',
              'published params.csv')

