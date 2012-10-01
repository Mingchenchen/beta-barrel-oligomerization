from Bio.Align.Applications import ClustalwCommandline
import Bio.AlignIO
import Bio.SeqIO
import Bio.PDB
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import warnings
import os
import copy
import csv

from sundries import one_letter

# The guts
def call_clustalw(seq_file, mat_file, output_file):
    '''The function that actually produces the sequence alignment.'''
    cline = ClustalwCommandline('clustalw', infile=seq_file,
                                matrix=mat_file, outfile = output_file)
    
    return cline()

def cluster_retriever(cluster_name):
    '''Returns a multiple sequence alignment corresponding to the given
    HHOMP cluster name.'''
    if cluster_name[:3].upper() == 'OMP' or 'cluster' in cluster_name:
        path = os.path.dirname(__file__) \
               + '/clusters/{0}.clu'.format(cluster_name)
    
    else:
        path = os.path.dirname(__file__) \
               + '/clusters/OMP.{0}.clu'.format(cluster_name)

    return Bio.AlignIO.read(path, 'clustal')

def sequence_retriever(pdb_paths):
    '''Given an iterable of paths of PDB files, return a list of
    SeqRecord objects containing their sequences.'''

    output = list()

    for path in pdb_paths:
        # Come up with a name to identify the sequence by
        # Has special routines for names of the format "aligned_{PDBID}.pdb"
        # (one of Daniel's aligned structures, used as a template), and
        # "{PDBID} aligned to daniel's {PDBID}.pdb", a structure to be
        # modeled when deriving zdiff's
        name = None
        template_pattern = re.compile(r"aligned_(....)\.pdb")
        target_pattern = re.compile(r"(....) aligned to daniel's ....\.pdb")
        filename = os.path.basename(path)
        template_match = re.match(template_pattern, filename)
        target_match = re.match(target_pattern, filename)
        if template_match is not None:
            name = 'template_' + template_match.group(1)
        elif target_match is not None:
            name = 'target_' + target_match.group(1)
        else:
            name = filename

        # I think ClustalX ignores everything after a space which is
        # annoying, so replace the spaces with underscores
        name = name.replace(" ", "_")
            
        # Load the structure. See what happens if you try to load one
        # of Daniel's aligned structures without hiding warnings. It's not
        # even about anything significant, for this kind of thing, 
        # just missing b factors and occupancies
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            stru = Bio.PDB.PDBParser().get_structure(name, path)
            
        # Get the sequence from the structure
        seq_as_str = str('')
        for r in stru.get_residues():
            # one_letter is a case-insensitive dictionary mapping
            # 3-letter residue names to 1-letter residue names
            # Biopython PDB structures use 3-letter residue names
            standard_20 = [n.upper() for n in one_letter.keys()]
            resn = r.get_resname().upper()

            # Only take stuff that's actually a residue, and not like
            # some bound lipid or something
            if resn in standard_20:
                seq_as_str += one_letter[resn]

            # In the PDB entry for 1FEP, the FASTA file just has "M"
            # where there's selenomethionine in the structure, 
            # which makes sense because Vik
            # says that putting in the selenium is something that's done
            # to make structure determination easier, a bacterium's FepA
            # is not really going to be full of selenium.
            # So, treat selenomethionine as methionine:
            elif resn == 'MSE':
                seq_as_str += 'M'
        
        seq_as_Seq = Seq(seq_as_str, alphabet = IUPACProtein())

        # Append a SeqRecord to the output
        output.append(SeqRecord(seq_as_Seq, name=name, id=name,
                      description='sequence taken from a structure file'))
        
    return output

def combine(output_filename, msa, seqs):
    '''Write a file containing all the sequences in a given MSA and given
    iterable object of sequences'''
    output_seqs = list(seqs)

    for record in msa:
        # An HHOMP alignment has lots of useful stuff besides sequences,
        # which I want to ignore here
        # I don't know what "gi" means but all the real sequences seem to
        # have it at the beginning of their id, so,
        # Filter for only those sequences whose id's begin with "gi"
        if record.id[:2] == 'gi':
            # Get rid of gaps
            gapless_as_str = ''.join(filter(lambda x: x!= '-', record))
            gapless_seq = Seq(gapless_as_str)
            gapless_record = copy.deepcopy(record)
            gapless_record.seq = gapless_seq
            output_seqs.append(gapless_record)
    
    Bio.SeqIO.write(output_seqs, output_filename, 'fasta')

def bbtm_filename(t):
    # Make string of t padded with 0's:
    t_string = (3-len(str(t)))*'0' + str(t)

    # Create filename of one of the matrices David Jimenez-Morales
    # sent me
    filename = os.getcwd() + r'\TMout\bbtmout\bbTMout'+t_string

    return filename


def write_series(id_t_triplets, output_filename):
    # Information on the series file format can be found here:
    # http://bips.u-strasbg.fr/en/Documentation/ClustalX/
    # In case the link breaks, the Wayback Machine has an archived page
    # from July 27 2011
    with open(output_filename, 'w') as o:
        o.write('CLUSTAL_SERIES\n')
        o.write('\n')
        for id_bottom, id_top, t in id_t_triplets:
            o.write(' '.join(["MATRIX", str(id_bottom), str(id_top),
                              bbtm_filename(t)]))
            o.write('\n')


# The API
class UnrecognizedMatrix(Exception):
    pass

def align(output_name, cluster, matrix, *pdbpaths):
    '''Make a multiple sequence alignment using ClustalW. Takes as input
    the name of a cluster, a matrix, and the paths of PDB files. The
    advantage of using this function over a Biopython ClustalwCommandline
    object (which this function will eventually initialize) is that it's
    specialized for making alignments between PDB files and HHOMP clusters.
    '''
    # Find the directory to write files to 
    # (for example, the file of unaligned sequences that is created
    # as an intermediate step. These files are not deleted, so the use of
    # this function leaves some trash behind. I'll fix it someday maybe)
    target_dir = os.path.dirname(output_name)

    # Set infile parameter
    if target_dir == '':
        infile_param = 'seqs.fasta'
    else:
        # If you did this and target_dir was '', you'd end up trying
        # to write to the root directory '/'
        infile_param = target_dir + '/seqs.fasta'
    cluster_msa = cluster_retriever(cluster)
    pdb_sequences = sequence_retriever(pdbpaths)
    combine(infile_param, cluster_msa, pdb_sequences)
    
    # Set matrix parameter
    # It could be a matrix built into clustal:
    if matrix.lower() in ('blosum', 'pam', 'id', 'gonnet'):
        matrix_param = matrix
    # Or it could be BBTM:
    elif matrix.lower() == 'bbtm':
        # These pairings of sequence identities and evolutionary distances
        # are taken from a 1994 Clustal paper. They're the same as the
        # pairings in that paper 
        # of sequence identities and PAM numbers for the PAM matrices
        series_path = target_dir + r'\bbtm_series.txt'
        write_series([( 0, 40 , 350),
                      (41, 60 , 120),
                      (61, 80 , 60 ),
                      (81, 100, 20 )], series_path)
        matrix_param = os.getcwd() + '\\' + series_path
    # But otherwise I don't know what it is
    else:
        raise UnrecognizedMatrix('Matrix "{0}" not '.format(matrix) \
                                 + 'recognized. Try BBTM,  PAM, '\
                                 + 'GONNET, BLOSUM or ID.')

    return call_clustalw(infile_param, matrix_param, output_name)
