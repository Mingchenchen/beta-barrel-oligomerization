from Bio.PDB import PDBParser
import warnings
from Bio import AlignIO
import os
import re

import alignments
from sundries import CIDict

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
        '''Return the zdiff between the residues at this position for the
        two structures requested. Returns None if one or both of the 
        positions has a gap.'''
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
        '''Initialize with any number of pairs (structure, sequence)
        where "structure" is a Biopython.PDB structure and
        "sequence" is the sequence of that structure, with gaps,
        taken from a multiple sequence alignment'''

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
        # Raise exception if the target was not found, and the
        # return statement never reached:
        raise NotFoundError('resi {0} of {1} not found' \
                            .format(resi, unknown_id))
    
    def report(self, template_id, unknown_id):
        '''Return a list of zdiffs; use resi_report for a list with
        residue numbers included. Sign of zdiff is positive if predicted
        is too high'''
        output = (pos.zdiff(template_id, unknown_id) \
                  for pos in self.positions)
        return filter(lambda x: x is not None, output)

    def resi_report(self, template_id, unknown_id):
        '''Return a list of pairs, (resi, zdiff), where "resi" is
        a residue number in the sequence whose structure is being
        predicted. Residues that are gapped in either sequence are
        not included. Sign of zdiff is positive if predicted is too high'''
        output = ((pos.resi(unknown_id),
                   pos.zdiff(template_id, unknown_id)) \
                  for pos in self.positions)
        return filter(lambda x: x[1] is not None, output)
        


def zdiff_align(matrix, output_dir):
    walk = list(os.walk('comparison structures'))
    # Get the names of the clusters:
    # [0][1] is the list of foldernames in the rootmost directory
    clusters = walk[0][1]

    for path, folder_list, file_list in walk[1:]:
        # Establish what cluster we're talking about
        for name in clusters:
            if name in path:
                cluster = name

        for filename in file_list:
            # Align each pair of an HHOMP structure, and the structure in
            # our dataset that was matched to the same cluster, with
            # all the sequences of the cluster
            match = re.match("(....) aligned to daniel's (....)\.pdb",
                             filename)
            if match is None:
                continue
            
            # Retrieve from the regex match the pdbid of the protein in
            # HHOMP:
            hhomp_pdbid = match.group(1)

            # Retrieve from the regex match the pdbid of the protein in 
            # our dataset that was mapped to this cluster:
            our_pdbid = match.group(2)
            
            # Retrieve the path names of the PDB files representing both
            hhomp_path = path + '/' + match.group(0)
            our_path = path + '/aligned_{0}.pdb'.format(our_pdbid)

            # Make the alignment
            # Create the output directory if it does not already exist
            try:
                # This will not successfully make the directory if
                # nested directories would have to be created
                os.mkdir(output_dir)
            except OSError:
                # This is the error you'd get if the directory already
                # exists
                # So, our work is done, don't worry about it:
                pass

            output_path = output_dir + '/{}, {} as target, {} as template'\
                                        .format(cluster, hhomp_pdbid,
                                               our_pdbid)
            alignments.align(output_path, cluster, matrix,
                             hhomp_path, our_path)

                             

            

def calc(template_name, target_name, template_structure_filename,
         target_structure_filename, alignment_filename,
         write_to, comment, format_='clustal'):
    '''Calculate zdiffs from an alignment and write them to a file.
    template_name and target_name: the id's of the template and target
        sequences in the alignment
    template_structure_filename and target_structure_filename: where to
        find the PDB-format structure files of the structures of the
        template and target
    alignment_filename: where to find an alignment containing both the
        template and the target
    write_to: where to put the resulting file of zdiffs
    comment: a note to put in the comments of the zdiff files. Record
        the scientific meaning of these zdiffs so they can still be useful
        as your memories about them fade.
    '''
    # Open relevant files
    try:
        template_structure_file = open(template_structure_filename, 'r')
        target_structure_file = open(target_structure_filename, 'r')
        alignment_file = open(alignment_filename, 'r') 
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
    except Exception:
        print('problems opening template structure file ' \
              + template_structure_filename)
        print('and target structure file ' + target_structure_filename)
        print('and alignment file ' + alignment_filename)
        raise

    finally:
        for i in (template_structure_file, target_structure_file,
                  alignment_file):
            i.close()
        
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
        for line in comment.split('\n'):
            o.write('# ' + line + '\n')

        # Write the actual data
        # Weird quirk of the zdiff objects - before giving a report, it
        # requires the id's of the structures. I was planning on making
        # a zdiff object support more than two structures, so you give
        # it the id's to tell it which pair you want, out of all those
        # in the object. I may even have implemented this.
        for resi, zdiff in results.resi_report(template_name, target_name):
            o.write(str(resi) + ', ' + str(zdiff) + '\n')

def calc_all(input_dir, output_dir, comment, format_='clustal'):
    '''Calculate zdiffs for all alignments in a given folder. Add to 
    each zdiff the given comment (multiline comments will still work).
    
    There are some harsh requirements on the names and contents
    of the alignments. The names must be of the form
    '{clustername}, {pdbid} as target, {pdbid} as template', followed by
    whatever file extension. The format must be whatever is given in the
    "format_" keyword argument (clustal by default). The sequence of the
    target protein (the one of which a model is being made) must be
    in the alignment with the id "target_{pdbid}", and the sequence
    of the protein whose structure is to be used as a template must
    be present with the id "template_{pdbid}".'''

    # Make the directory output_dir
    # Does nothing if it already exists
    # Does nothing if more than one directory would have to be created,
    # leading to errors downstream
    # This is bad coding. You never want to catch an exception in a 
    # way such taht two different errors get handled the same way,
    # when there are two different appropriate responses. I may fix
    # this later.
    try:
        os.mkdir(output_dir)
    except OSError:
        pass

    for alignment_filename in os.listdir(input_dir):
        # Check if it's a zdiff alignment, and if it is, retrieve the name
        # of the cluster, the PDBID of the target to be modelled, and the
        # PDBID of the template
        match = re.match(r"(.*?), (....) as target, (....) as template",
                         alignment_filename)
        if match is None:
            continue
        else:
            cluster = match.group(1)
            target = match.group(2)
            template = match.group(3)

        # Find the pdb files. Need z values to calculate zdiffs.
        target_path = "comparison structures/" \
                      + "{}/{} aligned to daniel's {}.pdb".format(cluster,
                                                                  target,
                                                                  template)
        template_path = "comparison structures/" \
                        + "{}/aligned_{}.pdb".format(cluster, template)

        # Construct the path and filename of the zdiff file to be created
        zdiff_path = output_dir + '/' + match.group(0) + '.zdiff'

        # Create the zdiffs
        # Remember: "template" and "target" are PDBID's
        calc('template_'+template, 'target_'+target, template_path,
             target_path, input_dir + '/' + alignment_filename, zdiff_path,
             comment, format_=format_)

