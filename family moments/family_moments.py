from __future__ import division
from sundries import one_letter
from sundries import CIDict
import warnings
from Bio.PDB import PDBParser
import Bio.AlignIO
import DSSP_win
import zenergy
import csv
import numpy as np
import re
import glob
import math
#from pymol import cmd

def draw_vector(name, vector, center):
    """Create a distance with the length and direction of a given vector,
    and then hide its label.
    The center should have the same xy projection as the center you gave to
    the moment function, but you might want to change the z value to make it
    more visible.
    
    ARGUMENTS:
    (name, vector_, center)
    RETURNS:
    Nothing."""
    
    cmd.pseudoatom('ps1', pos = tuple(center))
    cmd.pseudoatom('ps2', pos = tuple(center + vector))
    cmd.distance(name, 'ps1', 'ps2')
    cmd.hide('labels', name)
    cmd.delete('ps1 ps2')


def geomed(starlist, quiet=True):
    '''Calculate the geometric median aka Fermat point using Weiszfeld's
    algorithm. You can learn more about this from the Wikipedia article
    "Geometric Median"'''
    # Throughout this function's comments, the metaphor is maintained
    # that this function moves a spaceship from the origin to the
    # geometric median of the surrounding stars by a series of smaller steps

    # Start with the spaceship at the origin
    ship = np.zeros(len(starlist[0]))

    for i in range(1000):
        # Create a matrix with the star positions as columns
        pos_mat = np.array(starlist).transpose()
    
        # Create a vector where the nth entry is the reciprocal of the
        # distance between the ship and the nth star
        inverse_distances = np.array([np.linalg.norm(star - ship)**-1 \
                                      for star in starlist])

        # Find the next spot dictated by Weiszfeld's algorithm
        destination = inverse_distances.sum()**-1 \
                      * np.dot(pos_mat, inverse_distances)
    
        # Break the loop if we're moving a billionth of an angstrom at a
        # time
        if np.linalg.norm(ship - destination) < 1e-12:
            break

        # Move the ship to the next spot
        ship = destination

    if not quiet:
        print(str(i) + ' iterations')

    return ship

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
    how to get each piece of information from the Structure object.

    Attributes:
        ca_coord: c alpha coordinate in three dimensions
        rel_acc: relative accessibility from DSSP

    Method:
        ez_b(self, seq_id): calculate Ez-beta using a particular sequence
                            in the family

    Two more attributes for if you need other information for some reason:
        residue: the Biopython residue object
        family: the Family object that this residue came from'''

    def __init__(self, residue, family):
        
        self.residue = residue
        self.family = family

        self.ca_coord = residue.child_dict['CA'].get_coord()

        # Get relative accessibility and secondary structure from
        # DSSP
        chain_id = residue.get_full_id()[2]
        extended_resi = residue.get_full_id()[3]
        self.rel_acc = family.dssp[(chain_id, extended_resi)][3]
        self.ss = family.dssp[(chain_id, extended_resi)][1]

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

class Selection(list):
    def __init__(self, family, inclusion_condition):

        filing_cabinet = [ResidueDossier(res, family)\
                          for res in family.stru.get_residues()]
        
        list.__init__(self, (dos for dos in filing_cabinet \
                             if inclusion_condition(dos)))
        
        self.family = family

        # Calculate the geometric median and save it right now, 
        # so that it doesn't have to be recalculated
        # every time the moment is calculated with a different sequence
        xy_coords = [dos.ca_coord[:2] for dos in self]
        self.geomed = geomed(xy_coords)
    
    def show(self):
        cmd.reinitialize()
        cmd.load(self.family.stru_path, self.family.stru_name)
        cmd.color('purple', '*')
        for dos in self:
            cmd.color('green', 'resi ' + str(dos.resi))

    def moment(self, seq_id):
        running_total = np.zeros(2)

        for dos in self:
            # Calculate ez-beta; skip this residue if this position is a gap
            # or if there is no information on this kind of residue
            try:
                ez_b = dos.ez_b(seq_id)
            except zenergy.NoParameters:
                continue
            
            # Calculate vector that points from the geometric median of the
            # structure, to this residue's c-alpha, with magnitude
            # equal to it's ez-beta
            # Then, add it to the running total
            relative_position = dos.ca_coord[:2] - self.geomed
            relative_direction = relative_position\
                                 / np.linalg.norm(relative_position)
            running_total += ez_b * relative_direction

        return running_total
    
    def show_moment(self, seq_id, normalize=False, name = None, z=0):
        if name is None:
            name = seq_id

        mx, my = self.moment(seq_id)
        if not normalize:
            moment = np.array([mx, my, 0])
        if normalize:
            moment = (mx**2 + my**2 + z**2)**(-.5) * np.array([mx,my,z])

        gx, gy = self.geomed
        geomed = np.array([gx, gy, z])
        draw_vector(name, moment, geomed)

    def spreadsheet(self, path):
        with open(path, 'wb') as f:
            csv.writer(f).writerow(['Seq id', 'X', 'Y', 'Direction',
                                    'Magnitude'])
            for seq_id in self.family.msa.keys():
                moment = self.moment(seq_id)
                direction = math.atan2(moment[1], moment[0]) * (180/math.pi)
                magnitude = np.linalg.norm(moment)
                csv.writer(f).writerow([seq_id, moment[0], moment[1],
                                        direction, magnitude])

def families(params='published params.csv', sanity_file=None):
    # Map PDBIDs to paths of structures
    stru_path_list = glob.glob('structures/aligned_*.pdb')
    match_str = r'structures[/\\]aligned_(....)\.pdb'
    pdbids = [re.match(match_str, path).group(1) \
              for path in stru_path_list]
    stru_path = CIDict(zip(pdbids, stru_path_list))
    
    # Map PDBIDS to paths of multiple sequence alignments
    msa_path_list = glob.glob('gonnet aligned/* with *.clu')
    match_str = r'gonnet aligned[/\\](....) with .*\.clu'
    pdbids = [re.match(match_str, path).group(1) for path in msa_path_list]
    msa_path = CIDict(zip(pdbids, msa_path_list))
    
    # Map PDBIDs to names of sequences of the structure
    template_id = CIDict((pdbid, 'template_' + pdbid.upper())
                          for pdbid in msa_path.keys())
    
    if sanity_file is not None:
        with open(sanity_file, 'w') as f:
            f.write(repr(stru_path))
            f.write('\n')
            f.write(repr(msa_path))
            f.write('\n')
            f.write(repr(template_id))
    # Create the families
    return CIDict((pdbid, Family(pdbid, stru_path[pdbid], msa_path[pdbid],
                          template_id[pdbid], params))
                   for pdbid in msa_path.keys())

fams = families()
scry_exposed = Selection(fams['1a0s'], lambda x: x.rel_acc > .2)
scry_exposed.show()
scry_exposed.spreadsheet('test.csv')
