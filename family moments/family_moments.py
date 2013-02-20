from sundries import one_letter
import warnings
from Bio.PDB import PDBParser
import Bio.AlignIO
import DSSP_win
import zenergy
import csv
#from pymol import cmd

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

class Family(object):
    '''A class that compiles information about a protein structure and 
    related sequences. This information is meant to be sufficient to filter
    the residues and calculate an ez-beta moment from those that remain.

    Attributes:
    stru_name: Name of the structure
    stru_path: Path of the PDB format structure file
    stru: the structure, as a Biopython entity
    msa: a dictionary mapping sequence identifiers to rows in a multiple
        sequence alignment
    template_seq: the row of the MSA containing the sequence of the
        structure
    res_to_pos: a dictionary mapping residues from structures to their
        column number in the MSA (asssuming the first column is numbered 0)
    dssp: a Biopython DSSP object with a DSSP for the structure
    calc: an Ez-beta calculator'''

    def __init__(self, stru_name, stru_path, msa_path, template_name,
                 param_path):
        '''Requires a name for the structure (your choice), a path to a
        PDB format structure file, a path to a multiple sequence alignment
        containing a row with exactly the same sequence as the structure,
        the sequence identifier of this row, and a path to a CSV file
        of Ez-beta parameters (see zenergy.Calculator for how to make
        these files)'''
        self.stru_name = stru_name
        self.stru_path = stru_path

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.stru = PDBParser().get_structure(stru_name, stru_path)

        # When Daniel created the aligned structures, he removed heteroatoms
        # (though the procedure he used seems to have removed anything
        # without the residue identifer of one of the 20 standard amino
        # acids, leading to main chain selenomethionines being removed from
        # the 1FEP structure). However, he also added a sort of box of
        # water atoms (perhaps as a visual aid, so you can tell how the
        # coordinate system is defined?)
        # Therefore, remove all waters
        waters = [i for i in self.stru.get_residues() \
                  if i.get_resname() == 'HOH']
        for chain in self.stru.get_chains():
            for water in waters:
                try:
                    chain.detach_child(water.get_id())
                # Maybe it's not in this chain
                except KeyError:
                    pass
            

        msa = Bio.AlignIO.read(open(msa_path), 'clustal')
        self.msa = dict((seq.id, seq) for seq in msa)
        self.template_seq = self.msa[template_name]

        self.res_to_pos = map_res_to_pos(self.stru.get_residues(),
                                         self.template_seq)
        
        self.dssp = DSSP_win.DSSP(self.stru.child_dict[0], stru_path)

        params = csv.reader(open(param_path, 'rb'))
        self.calc = zenergy.Calculator(params)

class ResidueDossier(object):
    '''A class that compiles the information about a residue that I expect
    to use directly in moment calculations, so I don't have to remember
    how to get each piece of information from the Structure object.'''
    def __init__(self, residue, family):
        
        self.residue = residue
        self.family = family

        self.ca_coord = residue.child_dict['CA'].get_coord()

        chain_id = residue.get_full_id()[2]
        extended_resi = residue.get_full_id()[3]
        self.rel_acc = family.dssp[(chain_id, extended_resi)][3]

        self.resi = residue.get_id()[1]

    def ez_b(self, seq_id):
        '''Calulate the Ez-Beta insertion energy of this residue, assuming
        the residue type that it has in the sequence with the given id.
        
        Will raise a NoParameters exception for a gap, or for a residue
        for which there are no parameters (either because they are too
        rare to determine their depth-dependent frequency trends, like
        cysteine, or because they apparently do not have depth-dependent
        frequency trends, like threonine).'''
        column = self.family.res_to_pos[self.residue]
        resn = self.family.msa[seq_id][column]

        z = self.ca_coord[2]
        
        return self.family.calc.calculate(resn, z)

class Selection(object):
    def __init__(self, family, inclusion_condition):
        self.family = family

        filing_cabinet = [ResidueDossier(res, family)\
                          for res in family.stru.get_residues()]
        
        # It seems very wrong to have a list as an object. A selection
        # is a container object, by nature. This code fails to represent
        # that and must be repaired.
        self.dossiers = [dos for dos in filing_cabinet \
                         if inclusion_condition(dos)]
    
    def show(self):
        cmd.reinitialize()
        cmd.load(self.family.stru_path, self.family.stru_name)
        cmd.color('purple', '*')
        for dos in self.dossiers:
            cmd.color('green', 'resi ' + str(dos.resi))
        
    



        



# Tests
test_family = Family('1A0S', 'structures/aligned_1A0S.pdb',
              'gonnet aligned/1A0S with cluster73.clu', 'template_1A0S',
              'published params.csv')
full_sele = Selection(test_family, lambda y: True)
